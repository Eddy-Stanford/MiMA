! (module added by mj)
!> The model time step `dt_atmos`, available throughout the code.
!>
!> `dt_atmos` is set by `coupler_nml`, which is read by the program `coupler_main`.
module coupler_mod
  integer :: dt_atmos = 500  !! [s] model time step
end module coupler_mod
!jm

!> Main program: couples the atmosphere to the surface and controls the time integration.
!>
!> Reads `coupler_nml` and sets the calendar and the start date: from `INPUT/coupler.res`
!> if it exists (unless `force_date_from_namelist`), otherwise from `current_date`. The
!> initial time of the experiment is the base date of the `diag_table` (or the start date
!> if that is zero); it must not be later than the start date. The run length is the sum of
!> `months`, `days`, `hours`, `minutes` and `seconds`, and should be a whole multiple of
!> `dt_atmos`.
!>
!> The atmosphere and the surface use the same time step and are coupled for the implicit
!> vertical diffusion of heat and moisture, with the atmosphere split into down and up
!> calls that correspond to the down and up sweeps of the tridiagonal elimination. Each
!> time step calls
!>
!> 1. `sfc_boundary_layer`: surface fluxes at the current surface temperature (without a
!>    coupled surface only if a diagnostic or a tracer needs them);
!> 2. `update_atmos_model_down`: radiation and vertical diffusion down to the surface;
!> 3. `flux_down_from_atmos`, `update_ocean_model_fast`, `update_land_model_fast` and
!>    `flux_up_to_atmos`: new ocean and land surface temperatures and the updated fluxes
!>    (skipped without a boundary layer, see `surface_is_coupled`);
!> 4. `update_atmos_model_up`: vertical diffusion up from the surface, moist processes and
!>    dynamics.
!>
!> The start and end dates are written to `time_stamp.out`; at the end the component
!> restart files and `RESTART/coupler.res` (calendar, initial and current date) are
!> written. The `RESTART` directory must exist.
!>
!> Namelist: `coupler_nml`
!> ([namelist reference](https://eddy-stanford.github.io/MiMA/Parameters/#coupler_nml)).
!>
!> Original authors: Bruce Wyman, V. Balaji.
program coupler_main
!-----------------------------------------------------------------------
!                   GNU General Public License
!
! This program is free software; you can redistribute it and/or modify it and
! are expected to follow the terms of the GNU General Public License
! as published by the Free Software Foundation; either version 2 of
! the License, or (at your option) any later version.
!
! MOM is distributed in the hope that it will be useful, but WITHOUT
! ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
! or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public
! License for more details.
!
! For the full text of the GNU General Public License,
! write to: Free Software Foundation, Inc.,
!           675 Mass Ave, Cambridge, MA 02139, USA.
! or see:   http://www.gnu.org/licenses/gpl.html
!-----------------------------------------------------------------------
  use constants_mod, only: constants_init
  use time_manager_mod, only: time_type, set_calendar_type, set_time, &
                              set_date, get_date, days_in_month, month_name, &
                              operator(+), operator(-), operator(<), &
                              operator(>), operator(/=), operator(/), &
                              operator(*), THIRTY_DAY_MONTHS, JULIAN, &
                              NOLEAP, NO_CALENDAR
  use fms_mod, only: input_nml_file, check_nml_error, &
                     uppercase, error_mesg, write_version_number, &
                     fms_init, fms_end
  use fms2_io_mod, only: file_exists

  use field_manager_mod, only: field_manager_init
  use diag_manager_mod, only: diag_manager_init, diag_manager_end, &
                              DIAG_OTHER, DIAG_ALL, get_base_date
!
! model interfaces: the atmosphere and the surface
!
  use atmos_model_mod, only: atmos_model_init, atmos_model_end, &
                             update_atmos_model_down, &
                             update_atmos_model_up, &
                             atmos_data_type, &
                             land_ice_atmos_boundary_type

  use surface_types_mod, only: surface_exchange_type
  use flux_exchange_mod, only: flux_exchange_init, sfc_boundary_layer, &
                               flux_down_from_atmos, flux_up_to_atmos, &
                               flux_exchange_end
  use ocean_model_mod, only: ocean_data_type, update_ocean_model_fast, &
                             update_ocean_model_slow
  use land_model_mod, only: land_data_type, update_land_model_fast

  use physics_driver_mod, only: surface_is_coupled

  use mpp_mod, only: mpp_clock_id, mpp_clock_begin, mpp_clock_end
  use mpp_mod, only: mpp_init, mpp_pe, mpp_npes, mpp_root_pe, &
                     stdlog, mpp_error, NOTE, FATAL, WARNING, &
                     mpp_set_current_pelist, mpp_declare_pelist

  use memutils_mod, only: print_memuse_stats

  use coupler_mod, only: dt_atmos

  implicit none

!-----------------------------------------------------------------------

  character(len=128) :: version = '$Id: coupler_main.f90,v 11.0.10.4 2005/05/25 15:47:50 pjp Exp $'
  character(len=128) :: tag = '$Name:  $'

!-----------------------------------------------------------------------
!---- model defined-types ----

  type(atmos_data_type) :: Atm
  type(land_ice_atmos_boundary_type) :: Land_ice_atmos_boundary
  type(ocean_data_type) :: Ocean
  type(land_data_type) :: Land
  type(surface_exchange_type) :: Exch

!-----------------------------------------------------------------------
! ----- coupled model time -----

  type(time_type) :: Time, Time_init, Time_end, Time_step_atmos
  type(time_type) :: Time_atmos
  integer :: num_atmos_calls, na

! ----- coupled model initial date -----

  integer :: date_init(6)
  integer :: calendar_type = -99

!-----------------------------------------------------------------------
!------ namelist interface -------

  integer, dimension(6) :: current_date = (/1, 1, 1, 0, 0, 0/)
  !! start date (year, month, day, hour, minute, second); used on a cold start, or on a restart if
  !! `force_date_from_namelist = .true.`
  character(len=17) :: calendar = 'thirty_day       '
  !! `'thirty_day'` (12 months of 30 days), `'julian'`, `'noleap'` or `'no_calendar'`
  logical :: force_date_from_namelist = .false.  !! take the date from `current_date` even when `INPUT/coupler.res` exists
  integer :: months = 0, days = 360, hours = 0, minutes = 0, seconds = 0
  !! `months`, `days`, `hours`, `minutes`, `seconds`: run length; it should be a whole multiple of
  !! `dt_atmos`
!mj exported dt_atmos into coupler_mod

  integer ::atmos_npes = 0  !! number of processes for the atmosphere (0: all)
  logical :: do_atmos = .true.  !! run the atmosphere
  namelist /coupler_nml/ current_date, calendar, force_date_from_namelist, months, days, hours, &
    minutes, seconds, dt_atmos, do_atmos, atmos_npes

  integer :: initClock, mainClock, termClock

  integer :: nfields
  character(len=80) :: text

!#######################################################################

  call mpp_init()
!these clocks are on the global pelist
  initClock = mpp_clock_id('Initialization')
  mainClock = mpp_clock_id('Main loop')
  termClock = mpp_clock_id('Termination')
  call mpp_clock_begin(initClock)

  call fms_init
  call constants_init

  call coupler_init
  call mpp_set_current_pelist()

  call mpp_clock_end(initClock) !end initialization

  call mpp_clock_begin(mainClock) !begin main loop

!-----------------------------------------------------------------------
!------ integration loop ------

  if (Atm%pe) then
    call mpp_set_current_pelist(Atm%pelist)
  end if
  call mpp_set_current_pelist()

  if (Atm%pe) then
    call mpp_set_current_pelist(Atm%pelist)

!-----------------------------------------------------------------------
!   ------ atmosphere/surface time step loop -------

    do na = 1, num_atmos_calls

      Time_atmos = Time_atmos + Time_step_atmos

      call sfc_boundary_layer(float(dt_atmos), Time_atmos, Atm, Ocean, Land, &
                              Land_ice_atmos_boundary)

!      ---- atmosphere down ----

      if (do_atmos) &
        call update_atmos_model_down(Land_ice_atmos_boundary, Atm)

      ! without a boundary layer (e.g. dry Held-Suarez) the surface
      ! state is left unchanged
      if (surface_is_coupled()) then
        call flux_down_from_atmos(Time_atmos, Atm, Exch)
        call update_ocean_model_fast(Exch, Ocean)
        call update_land_model_fast(Exch, Land)
        call flux_up_to_atmos(Time_atmos, Atm, Ocean, Land, Exch, Land_ice_atmos_boundary)
      end if

!      ---- atmosphere up ----

      if (do_atmos) &
        call update_atmos_model_up(Land_ice_atmos_boundary, Atm)

      call update_ocean_model_slow(Time_atmos, Ocean)

    end do

!   ------ end of atmospheric time step loop -----
!-----------------------------------------------------------------------

    Time = Time_atmos
  end if                     !Atm%pe block

  call mpp_set_current_pelist()
!-----------------------------------------------------------------------
  call mpp_clock_end(mainClock)
  call mpp_clock_begin(termClock)
  call coupler_end

  call diag_manager_end(Time)
  call mpp_clock_end(termClock)

  call print_memuse_stats('Memory HiWaterMark', always=.true.)
  call fms_end

!-----------------------------------------------------------------------

contains

!#######################################################################

  !> Reads `coupler_nml` and the restart date, sets the calendar, the PE list, the
  !> diagnostics manager and the run length, and initializes the atmosphere and the surface.
  subroutine coupler_init

    integer :: unit, ierr, io, m, i
    integer :: date(6)
    type(time_type) :: Run_length
    character(len=9) :: month
    integer :: pe, npes
    integer :: atmos_pe_start = 0, atmos_pe_end = 0
    integer :: diag_model_subset = DIAG_ALL
!-----------------------------------------------------------------------
!----- read namelist -------

    read (input_nml_file, nml=coupler_nml, iostat=io)
    ierr = check_nml_error(io, 'coupler_nml')

!----- write namelist to logfile -----
    call write_version_number(version, tag)
    if (mpp_pe() == mpp_root_pe()) write (stdlog(), nml=coupler_nml)

!----- read date and calendar type from restart file -----

    if (file_exists('INPUT/coupler.res')) then
!Balaji: currently written in binary, needs form=MPP_NATIVE
      open (newunit=unit, file='INPUT/coupler.res', form='formatted', action='read', status='old')
      read (unit, *, err=999) calendar_type
      read (unit, *) date_init
      read (unit, *) date
      goto 998 !back to fortran-4
!read old-style coupler.res
999   close (unit)
      open (newunit=unit, file='INPUT/coupler.res', form='unformatted', action='read', status='old')
      read (unit) calendar_type
      read (unit) date
998   close (unit)
    else
      force_date_from_namelist = .true.
    end if

!----- use namelist value (either no restart or override flag on) ---

    if (force_date_from_namelist) then

      if (sum(current_date) <= 0) then
        call error_mesg('program coupler', &
                        'no namelist value for base_date or current_date', FATAL)
      else
        date = current_date
      end if

!----- override calendar type with namelist value -----

      select case (uppercase(trim(calendar)))
      case ('JULIAN')
        calendar_type = JULIAN
      case ('NOLEAP')
        calendar_type = NOLEAP
      case ('THIRTY_DAY')
        calendar_type = THIRTY_DAY_MONTHS
      case ('NO_CALENDAR')
        calendar_type = NO_CALENDAR
      case default
        call mpp_error(FATAL, 'COUPLER_MAIN: coupler_nml entry calendar must be one of JULIAN|NOLEAP|THIRTY_DAY|NO_CALENDAR.')
      end select

    end if

    call set_calendar_type(calendar_type)

!----- write current/initial date actually used to logfile file -----

    if (mpp_pe() .eq. mpp_root_pe()) &
      write (stdlog(), 16) date(1), trim(month_name(date(2))), date(3:6)
16  format('  current date used = ', i4, 1x, a, 2i3, 2(':', i2.2), ' gmt')

!-----------------------------------------------------------------------
!------ initialize PEset management ------

!pe information
    pe = mpp_pe()
    npes = mpp_npes()

    if (atmos_npes .eq. 0) atmos_npes = npes
    atmos_pe_start = 0
    atmos_pe_end = atmos_npes - 1
!messages
    write (text, '(a,2i6)') 'Atmos PE range: ', atmos_pe_start, atmos_pe_end
    call mpp_error(NOTE, 'coupler_init: '//trim(text))
    allocate (Atm%pelist(atmos_npes))
    Atm%pelist = (/(i, i=atmos_pe_start, atmos_pe_end)/)
    Atm%pe = atmos_pe_start .le. pe .and. pe .le. atmos_pe_end
    call mpp_declare_pelist(Atm%pelist, '_atm')

!-----------------------------------------------------------------------
!------ initialize diagnostics manager ------

!jwd Fork here is somewhat dangerous. It relies on "no side effects" from
!    diag_manager_init. diag_manager_init or this section should be
!    re-architected to guarantee this or remove this assumption.
!    For instance, what follows assumes that get_base_date has the same
!    time for every pe.

    if (Atm%pe) then
      call mpp_set_current_pelist(Atm%pelist)
      if (atmos_npes /= npes) diag_model_subset = DIAG_OTHER  ! change diag_model_subset from DIAG_ALL
    end if
    call field_manager_init(nfields)
    call diag_manager_init(DIAG_MODEL_SUBSET=diag_model_subset)   ! initialize diag_manager for processor subset output
    call print_memuse_stats('diag_manager_init')
!-----------------------------------------------------------------------
!------ reset pelist to "full group" ------

    call mpp_set_current_pelist()
!----- always override initial/base date with diag_manager value -----

    call get_base_date(date_init(1), date_init(2), date_init(3), &
                       date_init(4), date_init(5), date_init(6))

!----- use current date if no base date ------

    if (date_init(1) == 0) date_init = date

!----- set initial and current time types ------

    Time_init = set_date(date_init(1), date_init(2), date_init(3), &
                         date_init(4), date_init(5), date_init(6))

    Time = set_date(date(1), date(2), date(3), &
                    date(4), date(5), date(6))

!----- compute the ending time -----

    Time_end = Time
    do m = 1, months
      Time_end = Time_end + set_time(0, days_in_month(Time_end))
    end do
    Time_end = Time_end + set_time(hours*3600 + minutes*60 + seconds, days)
    Run_length = Time_end - Time

!-----------------------------------------------------------------------
!----- write time stamps (for start time and end time) ------

    if (mpp_pe() .eq. mpp_root_pe()) &
      open (newunit=unit, file='time_stamp.out', form='formatted', action='write', status='replace')

    month = month_name(date(2))
    if (mpp_pe() .eq. mpp_root_pe()) write (unit, 20) date, month(1:3)

    call get_date(Time_end, date(1), date(2), date(3), &
                  date(4), date(5), date(6))
    month = month_name(date(2))
    if (mpp_pe() .eq. mpp_root_pe()) write (unit, 20) date, month(1:3)

    if (mpp_pe() .eq. mpp_root_pe()) close (unit)

20  format(6i4, 2x, a3)

!-----------------------------------------------------------------------
!----- compute the time steps ------

    Time_step_atmos = set_time(dt_atmos, 0)

!----- determine maximum number of iterations per loop ------

    num_atmos_calls = Run_length/Time_step_atmos

!-----------------------------------------------------------------------
!------------------- some error checks ---------------------------------

!----- initial time cannot be greater than current time -------

    if (Time_init > Time) call error_mesg('program coupler', &
                                          'initial time is greater than current time', FATAL)

!-----------------------------------------------------------------------
!------ initialize component models ------
!------ grid info now comes from grid_spec file

    if (Atm%pe) then
      call mpp_set_current_pelist(Atm%pelist)
!---- atmosphere ----
      call atmos_model_init(Atm, Time_init, Time, Time_step_atmos)
      call print_memuse_stats('atmos_model_init')

!---- surface: flux exchange, ocean and land ----
      call flux_exchange_init(Time, Atm, Ocean, Land, Land_ice_atmos_boundary, Exch)

    end if
    call mpp_set_current_pelist()

    Time_atmos = Time

    call mpp_set_current_pelist()

!-----------------------------------------------------------------------
!---- open and close dummy file in restart dir to check if dir exists --

    if (mpp_pe() .eq. mpp_root_pe()) then
      open (newunit=unit, file='RESTART/file', action='write', status='replace', iostat=io)
      if (io /= 0) call error_mesg('program coupler', &
                                   'cannot write to the RESTART directory; does it exist?', FATAL)
      close (unit, status='delete')
    end if

!-----------------------------------------------------------------------
    call print_memuse_stats('coupler_init')
  end subroutine coupler_init

!#######################################################################

  !> Terminates the surface and the atmosphere (writing their restart files) and writes
  !> `RESTART/coupler.res`.
  subroutine coupler_end

    integer :: unit, date(6)
!-----------------------------------------------------------------------

    call mpp_set_current_pelist()

!----- compute current date ------

    call get_date(Time, date(1), date(2), date(3), &
                  date(4), date(5), date(6))

!----- check time versus expected ending time ----

    if (Time /= Time_end) call error_mesg('program coupler', &
                                          'final time does not match expected ending time', WARNING)

!-----------------------------------------------------------------------
    if (Atm%pe) then
      call mpp_set_current_pelist(Atm%pelist)
! The surface first: its grid points to the grid of Atm, which atmos_model_end deallocates
      call flux_exchange_end(Ocean, Land)

      call atmos_model_end(Atm)
    end if
    call mpp_set_current_pelist()

!----- write restart file ------

    if (mpp_pe() .eq. mpp_root_pe()) then
      open (newunit=unit, file='RESTART/coupler.res', form='formatted', action='write', status='replace')
      write (unit, '(i6,8x,a)') calendar_type, &
        '(Calendar: no_calendar=0, thirty_day_months=1, julian=2, gregorian=3, noleap=4)'

      write (unit, '(6i6,8x,a)') date_init, &
        'Model start time:   year, month, day, hour, minute, second'
      write (unit, '(6i6,8x,a)') date, &
        'Current model time: year, month, day, hour, minute, second'
      close (unit)
    end if

!-----------------------------------------------------------------------

  end subroutine coupler_end

!#######################################################################

end program coupler_main
