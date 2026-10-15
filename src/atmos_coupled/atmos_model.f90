!> Driver for the atmospheric model: advances the atmospheric state by one time step.
!>
!> Designed around the implicit vertical diffusion scheme of the GCM, it needs two calls to
!> advance the model one time step. They correspond to the down and up sweeps of the
!> tridiagonal solver: `update_atmos_model_down` computes the radiation and the vertical
!> diffusion down to the surface, and `update_atmos_model_up` finishes the vertical
!> diffusion, computes the moist processes and steps the dynamics.
!>
!> The fields exchanged with the surface are held in derived types. A variable of type
!> `atmos_data_type` is returned by `atmos_model_init`; its contents should only be
!> modified by the atmospheric model. The precipitation and gustiness (and optionally the
!> lowest-level temperature and humidity) are saved in `RESTART/atmos_coupled.res.nc`.
!>
!> Namelist: `atmos_model_nml`
!> ([namelist reference](https://eddy-stanford.github.io/MiMA/Parameters/#atmos_model_nml)).
!>
!> Original authors: Bruce Wyman.
module atmos_model_mod
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

  !> State of the atmosphere seen by the surface and the coupler: grid, lowest-level
  !> fields, surface fluxes computed by the atmosphere, and time.
  type atmos_data_type
    type(domain2d)               :: domain             !! domain decomposition
    integer                       :: axes(4)            !! diagnostic axis ids of the grid (lon, lat, pfull, phalf)
    real, pointer, dimension(:)   :: glon_bnd => null() !! longitudes of the grid-box boundaries, global grid [rad]
    real, pointer, dimension(:)   :: glat_bnd => null() !! latitudes of the grid-box boundaries, global grid [rad]
    real, pointer, dimension(:)   :: lon_bnd => null() !! longitudes of the grid-box boundaries, this PE's subdomain [rad]
    real, pointer, dimension(:)   :: lat_bnd => null() !! latitudes of the grid-box boundaries, this PE's subdomain [rad]
    real, pointer, dimension(:, :) :: t_bot => null() !! temperature at the lowest model level [K]
    real, pointer, dimension(:, :) :: q_bot => null() !! specific humidity at the lowest model level [kg/kg]
    real, pointer, dimension(:, :) :: z_bot => null() !! height of the lowest model level above the surface [m]
    real, pointer, dimension(:, :) :: p_bot => null() !! pressure at the lowest model level [Pa]
    real, pointer, dimension(:, :) :: u_bot => null() !! zonal wind at the lowest model level [m/s]
    real, pointer, dimension(:, :) :: v_bot => null() !! meridional wind at the lowest model level [m/s]
    real, pointer, dimension(:, :) :: p_surf => null() !! surface pressure [Pa]
    real, pointer, dimension(:, :) :: gust => null() !! gustiness [m/s]
    real, pointer, dimension(:, :) :: flux_sw => null() !! net downward shortwave flux at the surface [W/m2]
    real, pointer, dimension(:, :) :: flux_lw => null() !! downward longwave flux at the surface [W/m2]
    real, pointer, dimension(:, :) :: lprec => null() !! liquid precipitation rate over the last time step [kg/m2/s]
    real, pointer, dimension(:, :) :: fprec => null() !! frozen precipitation rate over the last time step [kg/m2/s]
    type(surf_diff_type)         :: Surf_diff          !! data of the implicit vertical diffusion at the surface
    type(time_type)              :: Time               !! current time
    type(time_type)              :: Time_step          !! atmospheric time step
    type(time_type)              :: Time_init          !! initial time of the experiment
    integer, pointer              :: pelist(:) => null() !! PEs on which the atmosphere runs
    logical                       :: pe                 !! `.true.` if the atmosphere runs on this PE
  end type

  !> Fields passed from the surface to the atmosphere.
  !>
  !> Declared and allocated by `coupler_main`.
  type land_ice_atmos_boundary_type
    real, dimension(:, :), pointer :: t => null() !! surface temperature [K]
    real, dimension(:, :), pointer :: albedo => null() !! surface albedo
    real, dimension(:, :), pointer :: land_frac => null() !! land fraction of the grid box
    real, dimension(:, :), pointer :: dt_t => null()
    !! change of the lowest-level temperature for the up sweep of the vertical diffusion [K]
    real, dimension(:, :), pointer :: dt_q => null()
    !! change of the lowest-level specific humidity for the up sweep of the vertical diffusion [kg/kg]
    real, dimension(:, :), pointer :: u_flux => null() !! zonal surface stress on the atmosphere [Pa]
    real, dimension(:, :), pointer :: v_flux => null() !! meridional surface stress on the atmosphere [Pa]
    real, dimension(:, :), pointer :: dtaudu => null()
    !! derivative of the zonal surface stress with respect to the lowest-level zonal wind [kg/m2/s]
    real, dimension(:, :), pointer :: dtaudv => null()
    !! derivative of the meridional surface stress with respect to the lowest-level meridional wind
    !! [kg/m2/s]
    real, dimension(:, :), pointer :: u_star => null() !! friction velocity [m/s]
    real, dimension(:, :), pointer :: b_star => null() !! buoyancy scale [m/s2]
    real, dimension(:, :), pointer :: q_star => null() !! moisture scale [kg/kg] (not set by `simple_surface_mod`)
    real, dimension(:, :), pointer :: rough_mom => null() !! roughness length for momentum [m]
    real, dimension(:, :, :), pointer :: data => null() !! collective field for the named fields above (not used)
    integer                         :: xtype                   !! `REGRID`, `REDIST` or `DIRECT` (not used)
  end type land_ice_atmos_boundary_type

  !> Fields passed from the land alone to the atmosphere (none at present).
  type :: land_atmos_boundary_type
    real, dimension(:, :), pointer :: data => null() !! not used
  end type land_atmos_boundary_type

  !> Fields passed from the sea ice alone to the atmosphere (none at present).
  type :: ice_atmos_boundary_type
    real, dimension(:, :), pointer :: data => null() !! not used
  end type ice_atmos_boundary_type

!Balaji
  integer :: atmClock
!-----------------------------------------------------------------------

  character(len=128) :: version = '$Id: atmos_model.f90,v 12.0 2005/04/14 15:35:34 fms Exp $'
  character(len=128) :: tagname = '$Name: lima $'

!-----------------------------------------------------------------------
  logical           :: restart_tbot_qbot = .false.
  !! also store the lowest-level temperature and humidity (`t_bot`, `q_bot`) in `atmos_coupled.res.nc`
  namelist /atmos_model_nml/ restart_tbot_qbot

contains

!#######################################################################
  !> Computes the atmospheric tendencies of the radiation and of the vertical diffusion of
  !> momentum, heat, moisture and tracers.
  !>
  !> Called every time step. For heat and moisture only the downward sweep of the
  !> tridiagonal elimination is done, hence the name.
  subroutine update_atmos_model_down(Surface_boundary, Atmos)
!
!-----------------------------------------------------------------------
    type(land_ice_atmos_boundary_type), intent(inout) :: Surface_boundary
    !! fields passed from the surface to the atmosphere
    type(atmos_data_type), intent(inout) :: Atmos  !! atmospheric state

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

!#######################################################################
  !> Finishes the vertical diffusion of heat and moisture (upward sweep), computes the
  !> convective and large-scale tendencies and steps the dynamics.
  !>
  !> Called every time step. The atmospheric time is advanced by one step, the
  !> lowest-level fields of `Atmos` are updated and the global integrals are written.
  subroutine update_atmos_model_up(Surface_boundary, Atmos)

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

    type(land_ice_atmos_boundary_type), intent(in) :: Surface_boundary
    !! fields passed from the surface to the atmosphere
    type(atmos_data_type), intent(inout) :: Atmos  !! atmospheric state

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

!#######################################################################
  !> Initializes the atmospheric model.
  !>
  !> Reads `atmos_model_nml`, registers the tracers, initializes the atmosphere, allocates
  !> the fields of `Atmos` and reads `INPUT/atmos_coupled.res.nc` if it exists.
  subroutine atmos_model_init(Atmos, Time_init, Time, Time_step)

    type(atmos_data_type), intent(inout) :: Atmos  !! atmospheric state, allocated here
    type(time_type), intent(in) :: Time_init, Time, Time_step
    !! initial time of the experiment, current time and atmospheric time step

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

!#######################################################################
  !> Terminates the atmospheric model.
  !>
  !> Calls the termination routines of the atmosphere and of the global integrals, writes
  !> `RESTART/atmos_coupled.res.nc` and deallocates the fields of `Atmos`.
  subroutine atmos_model_end(Atmos)

    type(atmos_data_type), intent(inout) :: Atmos  !! atmospheric state
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

!#######################################################################

end module atmos_model_mod

