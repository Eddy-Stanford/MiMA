!> Driver for the atmospheric tracers: calls the source-sink routines of the individual tracers.
!>
!> The tracer and tracer tendency arrays are supplied with the longitude, latitude, wind,
!> temperature and pressure, from which the tracer routines compute the tendencies due to
!> emissions, chemical losses and deposition. The driver applies dry deposition to all
!> prognostic tracers and calls the radon routines for the radon tracers that are in the
!> `field_table`. Other tracers are passive.
!>
!> To add a tracer: add a `use` statement for its module; in `atmos_tracer_driver_init`, get
!> its index with `get_tracer_index(MODEL_ATMOS, 'name')` (positive if the tracer is in the
!> `field_table`) and, if it is positive, call the tracer's initialization routine; in
!> `atmos_tracer_driver`, call its source-sink routine for that index and add the returned
!> tendency to `rdt` (this is the user's responsibility); in `atmos_tracer_driver_end`, call
!> its termination routine.
!>
!> Original authors: William Cooke.
module atmos_tracer_driver_mod

!-----------------------------------------------------------------------

  use fms_mod, only: &
    write_version_number, &
    error_mesg, &
    FATAL, &
    mpp_pe, &
    mpp_root_pe, &
    stdlog
  use time_manager_mod, only: time_type
  use tracer_manager_mod, only: get_number_tracers
  use field_manager_mod, only: MODEL_ATMOS
  use atmos_tracer_utilities_mod, only: dry_deposition, &
                                        atmos_tracer_utilities_init
  use constants_mod, only: grav
  use atmos_radon_mod, only: atmos_radon_sourcesink, &
                             atmos_radon_init, &
                             atmos_radon_end

  implicit none
  private
!-----------------------------------------------------------------------
!----- interfaces -------

  public atmos_tracer_driver, atmos_tracer_driver_init, atmos_tracer_driver_end

!-----------------------------------------------------------------------
!  For each additional tracer, add a variable for its tracer index here.
!-----------------------------------------------------------------------

  integer, dimension(:), pointer :: nradon

  integer :: nt     ! number of activated tracers
  integer :: ntp    ! number of activated prognostic tracers

  logical :: module_is_initialized = .false.

!---- version number -----
  character(len=128) :: version = '$Id: atmos_tracer_driver.f90,v 11.0 2004/09/28 19:26:51 fms Exp $'
  character(len=128) :: tagname = '$Name: lima $'
!-----------------------------------------------------------------------

contains

!#######################################################################

  !> Computes the tracer tendencies due to dry deposition and to the sources and sinks of the
  !> individual tracers, and adds them to `rdt`.
  !>
  !> This is the interface between the dynamical core and the tracer code: it supplies the
  !> information needed to compute the tendency of each tracer due to emissions or chemical
  !> losses.
  subroutine atmos_tracer_driver(is, ie, js, je, Time, lon, lat, land, phalf, pfull, r, &
                                 u, v, t, q, u_star, rdt, rm, &
                                 dt, z_half, z_full, t_surf_rad, albedo, &
                                 Time_next, &
                                 kbot)

!-----------------------------------------------------------------------
    integer, intent(in)                           :: is, ie, js, je  !! local domain boundaries
    type(time_type), intent(in)                   :: Time  !! model time
    real, intent(in), dimension(:, :)           :: lon, lat, u_star
    !! `lon`, `lat`: longitude and latitude of the centres of the grid cells [rad]; `u_star`:
    !! friction velocity [m/s] (the magnitude of the wind stress is density times `u_star**2`)
    real, intent(in), dimension(:, :)           :: land  !! land fraction (land where > 0.5)
    real, intent(in), dimension(:, :, :)         :: phalf, pfull, u, v, t, q
    !! `phalf`, `pfull`: pressure at half and full levels [Pa]; `u`, `v`: zonal and meridional
    !! wind [m/s]; `t`: temperature [K]; `q`: specific humidity [kg/kg]
    real, intent(in), dimension(:, :, :, :)       :: r  !! tracer array
    real, intent(inout), dimension(:, :, :, :)       :: rdt
    !! tendency of the tracer array, to which the tracer tendencies are added
    real, intent(inout), dimension(:, :, :, :)       :: rm  !! tracer array at the previous time step
    real, intent(in)                              :: dt !! time step [s] (not used)
    real, intent(in), dimension(:, :, :)         :: z_half !! height at half levels [m]
    real, intent(in), dimension(:, :, :)         :: z_full !! height at full levels [m]
    real, intent(in), dimension(:, :)           :: t_surf_rad !! surface temperature [K]
    real, intent(in), dimension(:, :)           :: albedo  !! surface albedo
    type(time_type), intent(in)                    :: Time_next  !! time at the end of the step
    integer, intent(in), dimension(:, :), optional :: kbot  !! index of the lowest model level above the surface
!-----------------------------------------------------------------------
    real, dimension(size(r, 1), size(r, 2), size(r, 3)) :: rtnd, pwt
    real, dimension(size(r, 1), size(r, 2)) :: dsinku
    integer :: k, kd

    integer :: nnn
!-----------------------------------------------------------------------

    if (.not. module_is_initialized) &
      call error_mesg('Tracer_driver', 'tracer_driver_init must be called first.', FATAL)

!-----------------------------------------------------------------------
    kd = size(r, 3)

! Flux into the layer is assumed to be kg(tracer)/m2/s
! Tracers are transported as mixing ratio.
! Therefore need to convert flux to a mixing ratio tendency
! kg(tracer)/kg(air)/s = kg(tracer)/m2/s / [dz(m) * density of air(kg(air)/m3)]
! dz = dp/(rho*gravity)
! so dz*density of air = dp/gravity
! kg(tracer)/kg(air)/s = flux(tracer) / (dp/gravity)
    do k = 1, kd
      pwt(:, :, k) = (phalf(:, :, k + 1) - phalf(:, :, k))/grav
    end do
! WARNING pwt is dp/grav!!
! Go do the dry deposition of the tracers

    do k = 1, ntp
      call dry_deposition(k, is, js, u(:, :, kd), v(:, :, kd), t(:, :, kd), &
                          pwt(:, :, kd), pfull(:, :, kd), u_star, &
                          (land > 0.5), dsinku, r(:, :, kd, k), Time)
      rdt(:, :, kd, k) = rdt(:, :, kd, k) - dsinku
    end do

!
!--------------- compute radon source-sink tendency --------------------
    do nnn = 1, size(nradon(:))
      if (nradon(nnn) > 0) then
        if (nradon(nnn) > nt) call error_mesg('Tracer_driver', &
                                              'Number of tracers .lt. number for radon', FATAL)
        call atmos_radon_sourcesink(lon, lat, land, pwt, r(:, :, :, nradon(nnn)), &
                                    rtnd, Time, kbot)
        rdt(:, :, :, nradon(nnn)) = rdt(:, :, :, nradon(nnn)) + rtnd(:, :, :)
      end if
    end do

  end subroutine atmos_tracer_driver

!#######################################################################

  !> Initializes the tracer driver and the individual tracers that are in the `field_table`.
  !>
  !> The arguments are passed on to the individual tracer code, which may set initial values.
  !> The tracer manager provides a simple fixed or exponential profile if this is given in the
  !> field table; a more complicated profile should be set up in the initialization of the
  !> tracer code.
  subroutine atmos_tracer_driver_init(lonb, latb, r, axes, Time, phalf, mask)

!-----------------------------------------------------------------------
    real, intent(in), dimension(:)               :: lonb, latb  !! longitudes and latitudes of the cell corners [rad]
    real, intent(inout), dimension(:, :, :, :)         :: r  !! tracer fields (nlon, nlat, nlev, ntrace)
    type(time_type), intent(in)                                :: Time  !! model time
    integer, intent(in)                                :: axes(4)  !! diagnostic axes (lon, lat, pfull, phalf)
    real, intent(in), dimension(:, :, :)           :: phalf  !! pressure at half levels [Pa] (not used)
    real, intent(in), dimension(:, :, :), optional :: mask
    !! 1. above the ground, 0. below (nlon, nlat, nlev)

!-----------------------------------------------------------------------
!
!  When initializing additional tracers, the user needs to make changes
!
!-----------------------------------------------------------------------

    if (module_is_initialized) return

    call write_version_number(version, tagname)

!If we wish to automatically register diagnostics for wet and dry
! deposition, do it now.
    call atmos_tracer_utilities_init(lonb, latb, axes, Time)

!----- set initial value of radon ------------

    call atmos_radon_init(r, axes, Time, nradon, mask)

    call get_number_tracers(MODEL_ATMOS, num_tracers=nt, &
                            num_prog=ntp)

    module_is_initialized = .true.

  end subroutine atmos_tracer_driver_init

!#######################################################################

  !> Terminates the tracer driver and calls the termination routines of the individual tracers.
  subroutine atmos_tracer_driver_end

!-----------------------------------------------------------------------

    if (mpp_pe() /= mpp_root_pe()) return

    write (stdlog(), '(/,(a))') 'Exiting tracer_driver, have a nice day ...'

    call atmos_radon_end

    module_is_initialized = .false.

!-----------------------------------------------------------------------

  end subroutine atmos_tracer_driver_end

!######################################################################

end module atmos_tracer_driver_mod

