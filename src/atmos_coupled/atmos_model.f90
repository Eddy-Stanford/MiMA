module atmos_model_mod
!<CONTACT EMAIL="Bruce.Wyman@noaa.gov"> Bruce Wyman
!</CONTACT>
! <REVIEWER EMAIL="Zhi.Liang@noaa.gov">
!  Zhi Liang
! </REVIEWER>
!-----------------------------------------------------------------------
!<OVERVIEW>
!  Driver for the atmospheric model, contains routines to advance the
!  atmospheric model state by one time step.
!</OVERVIEW>

!<DESCRIPTION>
!     This version of atmos_model_mod has been designed around the implicit
!     version diffusion scheme of the GCM. It requires two routines to advance
!     the atmospheric model one time step into the future. These two routines
!     correspond to the down and up sweeps of the standard tridiagonal solver.
!     Most atmospheric processes (dynamics,radiation,etc.) are performed
!     in the down routine. The up routine finishes the vertical diffusion
!     and computes moisture related terms (convection,large-scale condensation,
!     and precipitation).

!     The boundary variables needed by other component models for coupling
!     are contained in a derived data type. A variable of this derived type
!     is returned when initializing the atmospheric model. It is used by other
!     routines in this module and by coupling routines. The contents of
!     this derived type should only be modified by the atmospheric model.

!</DESCRIPTION>

  use mpp_mod, only: mpp_pe, mpp_root_pe, mpp_clock_id, mpp_clock_begin
  use mpp_mod, only: mpp_clock_end, CLOCK_COMPONENT, mpp_error
  use mpp_domains_mod, only: domain2d
  use fms_mod, only: error_mesg, FATAL, NOTE
  use fms_mod, only: write_version_number, stdlog
  use fms_mod, only: clock_flag_default
  use fms_mod, only: input_nml_file, check_nml_error
  use restart_file_mod, only: restart_file_type, open_restart_read, open_restart_write
  use restart_file_mod, only: close_restart, read_restart_field, write_restart_field
  use time_manager_mod, only: time_type, operator(+), get_time
  use field_manager_mod, only: MODEL_ATMOS
  use tracer_manager_mod, only: register_tracers
  use mima_diag_integral_mod, only: diag_integral_init, diag_integral_end
  use mima_diag_integral_mod, only: diag_integral_output
  use atmosphere_mod, only: atmosphere_up, atmosphere_down, atmosphere_init
  use atmosphere_mod, only: atmosphere_end, get_bottom_mass, get_bottom_wind
  use atmosphere_mod, only: atmosphere_resolution, atmosphere_domain
  use atmosphere_mod, only: atmosphere_boundary, get_atmosphere_axes
  use atmosphere_mod, only: surf_diff_type

!-----------------------------------------------------------------------

  implicit none
  private

  public update_atmos_model_down, update_atmos_model_up
  public atmos_model_init, atmos_model_end, atmos_data_type
  public land_ice_atmos_boundary_type, land_atmos_boundary_type
  public ice_atmos_boundary_type
!-----------------------------------------------------------------------

!<PUBLICTYPE >
  type atmos_data_type
    type(domain2d)               :: domain             ! domain decomposition
    integer                       :: axes(4)            ! axis indices (returned by diag_manager) for the atmospheric grid
    ! (they correspond to the x, y, pfull, phalf axes)
    real, pointer, dimension(:)   :: glon_bnd => null() ! global longitude axis grid box boundaries in radians.
    real, pointer, dimension(:)   :: glat_bnd => null() ! global latitude axis grid box boundaries in radians.
    real, pointer, dimension(:)   :: lon_bnd => null() ! local longitude axis grid box boundaries in radians.
    real, pointer, dimension(:)   :: lat_bnd => null() ! local latitude axis grid box boundaries in radians.
    real, pointer, dimension(:, :) :: t_bot => null() ! temperature at lowest model level
    real, pointer, dimension(:, :) :: q_bot => null() ! specific humidity at lowest model level
    real, pointer, dimension(:, :) :: z_bot => null() ! height above the surface for the lowest model level
    real, pointer, dimension(:, :) :: p_bot => null() ! pressure at lowest model level
    real, pointer, dimension(:, :) :: u_bot => null() ! zonal wind component at lowest model level
    real, pointer, dimension(:, :) :: v_bot => null() ! meridional wind component at lowest model level
    real, pointer, dimension(:, :) :: p_surf => null() ! surface pressure
    real, pointer, dimension(:, :) :: gust => null() ! gustiness factor
    real, pointer, dimension(:, :) :: flux_sw => null() ! net shortwave flux (W/m2) at the surface
    real, pointer, dimension(:, :) :: flux_lw => null() ! net longwave flux (W/m2) at the surface
    real, pointer, dimension(:, :) :: lprec => null() ! liquid precipitation rate over the last time step (kg/m2/s)
    real, pointer, dimension(:, :) :: fprec => null() ! frozen precipitation rate over the last time step (kg/m2/s)
    type(surf_diff_type)         :: Surf_diff          ! store data needed by the multi-step version of the diffusion algorithm
    type(time_type)              :: Time               ! current time
    type(time_type)              :: Time_step          ! atmospheric time step.
    type(time_type)              :: Time_init          ! reference time.
    integer, pointer              :: pelist(:) => null() ! pelist where atmosphere is running.
    logical                       :: pe                 ! current pe.
  end type
!</PUBLICTYPE >

!<PUBLICTYPE >
  type land_ice_atmos_boundary_type
    ! variables of this type are declared by coupler_main, allocated by flux_exchange_init.
!quantities going from land+ice to atmos
    real, dimension(:, :), pointer :: t => null() ! surface temperature for radiation calculations
    real, dimension(:, :), pointer :: albedo => null() ! surface albedo for radiation calculations
    real, dimension(:, :), pointer :: land_frac => null() ! fraction amount of land in a grid box
    real, dimension(:, :), pointer :: dt_t => null() ! temperature tendency at the lowest level
    real, dimension(:, :), pointer :: dt_q => null() ! specific humidity tendency at the lowest level
    real, dimension(:, :), pointer :: u_flux => null() ! zonal wind stress
    real, dimension(:, :), pointer :: v_flux => null() ! meridional wind stress
    real, dimension(:, :), pointer :: dtaudu => null() ! derivative of zonal wind stress w.r.t. the lowest zonal level wind speed
    real, dimension(:, :), pointer :: dtaudv => null() ! derivative of meridional wind stress w.r.t. the lowest meridional level wind speed
    real, dimension(:, :), pointer :: u_star => null() ! friction velocity
    real, dimension(:, :), pointer :: b_star => null() ! bouyancy scale
    real, dimension(:, :), pointer :: q_star => null() ! moisture scale
    real, dimension(:, :), pointer :: rough_mom => null() ! surface roughness (used for momentum)
    real, dimension(:, :, :), pointer :: data => null() !collective field for "named" fields above
    integer                         :: xtype                   !REGRID, REDIST or DIRECT
  end type land_ice_atmos_boundary_type
!</PUBLICTYPE >

!<PUBLICTYPE >
  type :: land_atmos_boundary_type
    real, dimension(:, :), pointer :: data => null() ! quantities going from land alone to atmos (none at present)
  end type land_atmos_boundary_type
!</PUBLICTYPE >

!<PUBLICTYPE >
!quantities going from ice alone to atmos (none at present)
  type :: ice_atmos_boundary_type
    real, dimension(:, :), pointer :: data => null() ! quantities going from ice alone to atmos (none at present)
  end type ice_atmos_boundary_type
!</PUBLICTYPE >

!Balaji
  integer :: atmClock
!-----------------------------------------------------------------------

  character(len=128) :: version = '$Id: atmos_model.f90,v 12.0 2005/04/14 15:35:34 fms Exp $'
  character(len=128) :: tagname = '$Name: lima $'

!-----------------------------------------------------------------------
  logical           :: restart_tbot_qbot = .false.
  namelist /atmos_model_nml/ restart_tbot_qbot

contains

!#######################################################################
! <SUBROUTINE NAME="update_atmos_model_down">
!
! <OVERVIEW>
!   compute the atmospheric tendencies for dynamics, radiation,
!   vertical diffusion of momentum, tracers, and heat/moisture.
! </OVERVIEW>
!
!<DESCRIPTION>
!   Called every time step as the atmospheric driver to compute the
!   atmospheric tendencies for dynamics, radiation, vertical diffusion of
!   momentum, tracers, and heat/moisture.  For heat/moisture only the
!   downward sweep of the tridiagonal elimination is performed, hence
!   the name "_down".
!</DESCRIPTION>

!   <TEMPLATE>
!     call  update_atmos_model_down( Surface_boundary, Atmos )
!   </TEMPLATE>

! <IN NAME = "Surface_boundary" TYPE="type(land_ice_atmos_boundary_type)">
!   Derived-type variable that contains quantities going from land+ice to atmos.
! </IN>

! <INOUT NAME="Atmos" TYPE="type(atmos_data_type)">
!   Derived-type variable that contains fields needed by the flux exchange module.
!   These fields describe the atmospheric grid and are needed to
!   compute/exchange fluxes with other component models.  All fields in this
!   variable type are allocated for the global grid (without halo regions).
! </INOUT>

  subroutine update_atmos_model_down(Surface_boundary, Atmos)
!
!-----------------------------------------------------------------------
    type(land_ice_atmos_boundary_type), intent(inout) :: Surface_boundary
    type(atmos_data_type), intent(inout) :: Atmos

!-----------------------------------------------------------------------
    call mpp_clock_begin(atmClock)

    call atmosphere_down(Atmos%Time, Surface_boundary%land_frac, &
                         Surface_boundary%t, Surface_boundary%albedo, &
                         Surface_boundary%rough_mom, &
                         Surface_boundary%u_star, &
                         Surface_boundary%b_star, &
                         Surface_boundary%q_star, &
                         Surface_boundary%dtaudu, &
                         Surface_boundary%dtaudv, &
                         Surface_boundary%u_flux, &
                         Surface_boundary%v_flux, &
                         Atmos%gust, &
                         Atmos%flux_sw, &
                         Atmos%flux_lw, &
                         Atmos%Surf_diff)

!-----------------------------------------------------------------------

    call mpp_clock_end(atmClock)
  end subroutine update_atmos_model_down
! </SUBROUTINE>

!#######################################################################
! <SUBROUTINE NAME="update_atmos_model_up">
!
!-----------------------------------------------------------------------
! <OVERVIEW>
!   upward vertical diffusion of heat/moisture and moisture processes
! </OVERVIEW>

!<DESCRIPTION>
!   Called every time step as the atmospheric driver to finish the upward
!   sweep of the tridiagonal elimination for heat/moisture and compute the
!   convective and large-scale tendencies.  The atmospheric variables are
!   advanced one time step and tendencies set back to zero.
!</DESCRIPTION>

! <TEMPLATE>
!     call  update_atmos_model_up( Surface_boundary, Atmos )
! </TEMPLATE>

! <IN NAME = "Surface_boundary" TYPE="type(land_ice_atmos_boundary_type)">
!   Derived-type variable that contains quantities going from land+ice to atmos.
! </IN>

! <INOUT NAME="Atmos" TYPE="type(atmos_data_type)">
!   Derived-type variable that contains fields needed by the flux exchange module.
!   These fields describe the atmospheric grid and are needed to
!   compute/exchange fluxes with other component models.  All fields in this
!   variable type are allocated for the global grid (without halo regions).
! </INOUT>

  subroutine update_atmos_model_up(Surface_boundary, Atmos)

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

    type(land_ice_atmos_boundary_type), intent(in) :: Surface_boundary
    type(atmos_data_type), intent(inout) :: Atmos

!-----------------------------------------------------------------------
    call mpp_clock_begin(atmClock)

    Atmos%Surf_diff%delta_t = Surface_boundary%dt_t
    Atmos%Surf_diff%delta_q = Surface_boundary%dt_q

    call atmosphere_up(Atmos%Time, Surface_boundary%land_frac, Atmos%Surf_diff, &
                       Atmos%lprec, Atmos%fprec, Atmos%gust)

!   --- advance time ---

    Atmos%Time = Atmos%Time + Atmos%Time_step

    call get_bottom_mass(Atmos%t_bot, Atmos%q_bot, &
                         Atmos%p_bot, Atmos%z_bot, &
                         Atmos%p_surf)

    call get_bottom_wind(Atmos%u_bot, Atmos%v_bot)

!------ global integrals ------

    call diag_integral_output(Atmos%Time)

!-----------------------------------------------------------------------
    call mpp_clock_end(atmClock)

  end subroutine update_atmos_model_up
! </SUBROUTINE>

!#######################################################################
! <SUBROUTINE NAME="atmos_model_init">
!
! <OVERVIEW>
! Routine to initialize the atmospheric model
! </OVERVIEW>

! <DESCRIPTION>
!     This routine allocates storage and returns a variable of type
!     atmos_boundary_data_type, and also reads a namelist input and restart file.
! </DESCRIPTION>

! <TEMPLATE>
!     call atmos_model_init (Atmos, Time_init, Time, Time_step)
! </TEMPLATE>

! <IN NAME="Time_init" TYPE="type(time_type)" >
!   The base (or initial) time of the experiment.
! </IN>

! <IN NAME="Time" TYPE="type(time_type)" >
!   The current time.
! </IN>

! <IN NAME="Time_step" TYPE="type(time_type)" >
!   The atmospheric model/physics time step.
! </IN>

! <INOUT NAME="Atmos" TYPE="type(atmos_data_type)">
!   Derived-type variable that contains fields needed by the flux exchange module.
! </INOUT>

  subroutine atmos_model_init(Atmos, Time_init, Time, Time_step)

    type(atmos_data_type), intent(inout) :: Atmos
    type(time_type), intent(in) :: Time_init, Time, Time_step

    integer :: unit, ntrace, ntprog, ntdiag, ntfamily, i, j
    integer :: mlon, mlat, nlon, nlat
    integer :: ierr, io
    type(restart_file_type) :: rst
!-----------------------------------------------------------------------

!---- set the atmospheric model time ------

    Atmos%Time_init = Time_init
    Atmos%Time = Time
    Atmos%Time_step = Time_step

    read (input_nml_file, nml=atmos_model_nml, iostat=io)
    ierr = check_nml_error(io, 'atmos_model_nml')

!-----------------------------------------------------------------------
! how many tracers have been registered?
!  (will print number below)
    call register_tracers(MODEL_ATMOS, ntrace, ntprog, ntdiag, ntfamily)
    if (ntfamily > 0) call error_mesg('atmos_model', 'ntfamily > 0', FATAL)

!-----------------------------------------------------------------------
!  ----- initialize atmospheric model -----

    call atmosphere_init(Atmos%Time_init, Atmos%Time, Atmos%Time_step, &
                         Atmos%Surf_diff)

!-----------------------------------------------------------------------
!---- allocate space ----

    call atmosphere_resolution(mlon, mlat, global=.true.)
    call atmosphere_resolution(nlon, nlat, global=.false.)
    call atmosphere_domain(Atmos%domain)

    allocate (Atmos%glon_bnd(mlon + 1), &
              Atmos%glat_bnd(mlat + 1), &
              Atmos%lon_bnd(nlon + 1), &
              Atmos%lat_bnd(nlat + 1), &
              Atmos%t_bot(nlon, nlat), &
              Atmos%q_bot(nlon, nlat), &
              Atmos%z_bot(nlon, nlat), &
              Atmos%p_bot(nlon, nlat), &
              Atmos%u_bot(nlon, nlat), &
              Atmos%v_bot(nlon, nlat), &
              Atmos%p_surf(nlon, nlat), &
              Atmos%gust(nlon, nlat), &
              Atmos%flux_sw(nlon, nlat), &
              Atmos%flux_lw(nlon, nlat), &
              Atmos%lprec(nlon, nlat), &
              Atmos%fprec(nlon, nlat))

    do j = 1, nlat
      do i = 1, nlon
        Atmos%flux_sw(i, j) = 0.0
        Atmos%flux_lw(i, j) = 0.0
      end do
    end do
!-----------------------------------------------------------------------
!------ get initial state for dynamics -------

    call get_atmosphere_axes(Atmos%axes)

    call atmosphere_boundary(Atmos%glon_bnd, Atmos%glat_bnd, &
                             global=.true.)
    call atmosphere_boundary(Atmos%lon_bnd, Atmos%lat_bnd, &
                             global=.false.)

    call get_bottom_mass(Atmos%t_bot, Atmos%q_bot, &
                         Atmos%p_bot, Atmos%z_bot, &
                         Atmos%p_surf)

    call get_bottom_wind(Atmos%u_bot, Atmos%v_bot)

!-----------------------------------------------------------------------
!---- print version number to logfile ----

    call write_version_number(version, tagname)
!  write the namelist to a log file
    if (mpp_pe() == 0) then
      unit = stdlog()
      write (unit, nml=atmos_model_nml)
    end if

!  number of tracers
    if (mpp_pe() == mpp_root_pe()) then
      write (stdlog(), '(a,i3)') 'Number of tracers =', ntrace
      write (stdlog(), '(a,i3)') 'Number of prognostic tracers =', ntprog
      write (stdlog(), '(a,i3)') 'Number of diagnostic tracers =', ntdiag
    end if

!------ read initial state for several atmospheric fields ------

    if (open_restart_read(rst, 'INPUT/atmos_coupled.res.nc', Atmos%domain)) then
      if (mpp_pe() == mpp_root_pe()) call mpp_error('atmos_model_mod', &
                                                    'Reading netCDF formatted restart file: INPUT/atmos_coupled.res.nc', NOTE)
      ! lprec and fprec are rates, so a change of time step needs no conversion
      call read_restart_field(rst, 'lprec', Atmos%lprec)
      call read_restart_field(rst, 'fprec', Atmos%fprec)
      call read_restart_field(rst, 'gust', Atmos%gust)

      if (restart_tbot_qbot) then
        call read_restart_field(rst, 't_bot', Atmos%t_bot)
        call read_restart_field(rst, 'q_bot', Atmos%q_bot)
      end if
      call close_restart(rst)
    else
      Atmos%lprec = 0.0
      Atmos%fprec = 0.0
      Atmos%gust = 1.0
    end if

!------ initialize global integral package ------

    call diag_integral_init(Atmos%Time_init, Atmos%Time, &
                            Atmos%lon_bnd, Atmos%lat_bnd)

!-----------------------------------------------------------------------
    atmClock = mpp_clock_id('Atmosphere', flags=clock_flag_default, grain=CLOCK_COMPONENT)
  end subroutine atmos_model_init
! </SUBROUTINE>

!#######################################################################
! <SUBROUTINE NAME="atmos_model_end">
!
! <OVERVIEW>
!  termination routine for atmospheric model
! </OVERVIEW>

! <DESCRIPTION>
!  Call once to terminate this module and any other modules used.
!  This routine writes a restart file and deallocates storage
!  used by the derived-type variable atmos_boundary_data_type.
! </DESCRIPTION>

! <TEMPLATE>
!   call atmos_model_end (Atmos)
! </TEMPLATE>

! <INOUT NAME="Atmos" TYPE="type(atmos_data_type)">
!   Derived-type variable that contains fields needed by the flux exchange module.
! </INOUT>

  subroutine atmos_model_end(Atmos)

    type(atmos_data_type), intent(inout) :: Atmos
    integer :: sec, day, dt
    type(restart_file_type) :: rst
!-----------------------------------------------------------------------
!---- termination routine for atmospheric model ----

    call atmosphere_end(Atmos%Time)

!------ global integrals ------

    call diag_integral_end(Atmos%Time)

!---- compute integer time step (in seconds) ----
    call get_time(Atmos%Time_step, sec, day)
    dt = sec + 86400*day

!------ write several atmospheric fields ------
!        and the time step

    if (mpp_pe() == mpp_root_pe()) then
      call mpp_error('atmos_model_mod', 'Writing netCDF formatted restart file.', NOTE)
    end if
    call open_restart_write(rst, 'RESTART/atmos_coupled.res.nc', Atmos%domain)
    call write_restart_field(rst, 'dt', real(dt))   ! not read any more; kept for older versions
    call write_restart_field(rst, 'lprec', Atmos%lprec)
    call write_restart_field(rst, 'fprec', Atmos%fprec)
    call write_restart_field(rst, 'gust', Atmos%gust)
    if (restart_tbot_qbot) then
      call write_restart_field(rst, 't_bot', Atmos%t_bot)
      call write_restart_field(rst, 'q_bot', Atmos%q_bot)
    end if
    call close_restart(rst)

!-------- deallocate space --------

    deallocate (Atmos%glon_bnd, &
                Atmos%glat_bnd, &
                Atmos%lon_bnd, &
                Atmos%lat_bnd, &
                Atmos%t_bot, &
                Atmos%q_bot, &
                Atmos%z_bot, &
                Atmos%p_bot, &
                Atmos%u_bot, &
                Atmos%v_bot, &
                Atmos%p_surf, &
                Atmos%gust, &
                Atmos%flux_sw, &
                Atmos%flux_lw, &
                Atmos%lprec, &
                Atmos%fprec)

!-----------------------------------------------------------------------

  end subroutine atmos_model_end
! </SUBROUTINE>

!#######################################################################

end module atmos_model_mod

