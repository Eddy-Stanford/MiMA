!> Gray radiation of Frierson, Held and Zurita-Gotor (2006).
!>
!> Longwave radiation is a two-stream gray scheme with a surface optical depth that varies
!> from `ir_tau_eq` at the equator to `ir_tau_pole` at the poles as `sin(lat)**2`; the
!> optical depth varies with pressure as `linear_tau * p/p00 + (1 - linear_tau) * (p/p00)**4`,
!> with `p00` = 1000 hPa. The insolation is annual-mean and zonally symmetric,
!> `solar_constant/4 * (1 + del_sol*P2(lat) + del_sw*sin(lat))`, plus two optional localized
!> perturbations; it is absorbed with the optical depth
!> `atm_abs * (1 - sw_diff*sin(lat)**2) * (p/p00)**4`, and the part reflected by the surface
!> leaves the atmosphere without absorption.
!>
!> Namelist: `gray_radiation_nml`
!> ([namelist reference](https://eddy-stanford.github.io/MiMA/Parameters/#gray_radiation_nml)).
!>
!> References:
!>
!> * Frierson, D. M. W., I. M. Held, and P. Zurita-Gotor, 2006: A gray-radiation aquaplanet
!>   moist GCM. Part I: Static stability and eddy scale. J. Atmos. Sci., 63, 2548-2566,
!>   https://doi.org/10.1175/JAS3753.1.
module gray_radiation_mod

! ==================================================================================
! ==================================================================================

!  use utilities_mod,         only: error_mesg, open_file, file_exist, &
!                                   check_nml_error, FATAL, get_my_pe, &
!                                   close_file, set_domain, get_num_pes, &
!                                   read_data, write_data, &
!                                   check_system_clock, NOTE, &
!                                   get_domain_decomp, check_system_clock

  use fms_mod, only: input_nml_file, check_nml_error, &
                     mpp_pe, mpp_root_pe, &
                     write_version_number, stdlog

  use constants_mod, only: stefan, cp_air, grav

  use diag_manager_mod, only: register_diag_field, send_data

  use time_manager_mod, only: time_type, set_date, set_time, &
                              get_time, operator(+), &
                              operator(-), operator(/=), get_date

!==================================================================================
  implicit none
  private
!==================================================================================

! version information

  character(len=128), parameter :: version = &
                                   '$Id: gray_radiation.f90,v 1.1.2.2 2005/05/21 02:01:56 pjp Exp $'

  character(len=128), parameter :: tagname = &
                                   '$Name:  $'

!==================================================================================

! public interfaces

  public :: gray_radiation_init, gray_radiation, gray_radiation_end
!==================================================================================

! module variables
  logical :: initialized = .false.
  real, parameter :: p00 = 1000.e2

  real    :: solar_constant = 1360.0  !! [W/m2] solar constant
  real    :: del_sol = 0.0  !! amplitude of the second Legendre polynomial in the insolation (equator-to-pole contrast)
! modif omp: winter/summer hemisphere
  real    :: del_sw = 0.0  !! amplitude of a hemispherically asymmetric insolation term (winter/summer hemisphere)
  real    :: ir_tau_eq = 4.0  !! longwave optical depth at the surface at the equator
  real    :: ir_tau_pole = 4.0  !! longwave optical depth at the surface at the poles
  real    :: atm_abs = 0.2  !! shortwave optical depth of the atmosphere
  real    :: sw_diff = 0.0  !! equator-to-pole reduction of the shortwave optical depth
  real    :: long_pert = 180.  !! [deg] longitude of the zonally localized insolation perturbation
  real    :: del_long = 30.  !! [deg] width of the zonally localized insolation perturbation
  real    :: size_pert = 0.  !! [W/m2] amplitude of a zonally localized insolation perturbation
  real    :: linear_tau = 0.1  !! fraction of the longwave optical depth that is linear in pressure (the rest goes as p<sup>4</sup>)

  real    :: lat_pert = 0.0  !! [deg] latitude of the centre of the localized (Walker-type) insolation perturbation
  real    :: lon_pert = 180.0  !! [deg] longitude of the centre of the localized (Walker-type) insolation perturbation
  real    :: del_lat = 30.0  !! [deg] latitudinal half-width of the localized (Walker-type) insolation perturbation
  real    :: del_lon = 90.0  !! [deg] longitudinal half-width of the localized (Walker-type) insolation perturbation
  real    :: fcng_pert = 0.0  !! [W/m2] amplitude of a localized (Walker-type) insolation perturbation

  real, save :: pi, deg_to_rad, rad_to_deg

  namelist /gray_radiation_nml/ solar_constant, del_sol, &
    ir_tau_eq, ir_tau_pole, atm_abs, sw_diff, long_pert, del_long, &
    size_pert, linear_tau, del_sw, &
    lat_pert, lon_pert, del_lat, del_lon, fcng_pert

!==================================================================================
!-------------------- diagnostics fields -------------------------------

  integer :: id_olr, id_swdn_sfc, id_swdn_toa, id_swnet_toa, id_lwdn_sfc, id_lwup_sfc, id_albedo, &
             id_tdt_rad, id_flux_rad, id_flux_lw, id_flux_sw, id_entrop_rad
!mj debug
  integer :: id_tau, id_tau_rad

  character(len=9), parameter :: mod_name = 'radiation'

  real :: missing_value = -999.

contains

! ==================================================================================
! ==================================================================================

  !> Initializes the module: reads `gray_radiation_nml` and registers the diagnostics.
  subroutine gray_radiation_init(axes, Time)

!-------------------------------------------------------------------------------------
    integer, intent(in), dimension(4) :: axes  !! diagnostic axes (lon, lat, pfull, phalf)
    type(time_type), intent(in)       :: Time  !! current time
!-------------------------------------------------------------------------------------
    integer, dimension(3) :: half = (/1, 2, 4/)
    integer :: ierr, io, unit
!-----------------------------------------------------------------------------------------
! read namelist and copy to logfile

    read (input_nml_file, nml=gray_radiation_nml, iostat=io)
    ierr = check_nml_error(io, 'gray_radiation_nml')

    call write_version_number(version, tagname)
    if (mpp_pe() == mpp_root_pe()) then
      write (stdlog(), nml=gray_radiation_nml)
    end if

    pi = 4.0*atan(1.)
    deg_to_rad = 2.*pi/360.
    rad_to_deg = 360.0/2./pi

    initialized = .true.

!-----------------------------------------------------------------------
!------- initialize quantities for integral package -------

!       call diag_integral_field_init ('olr',    'f8.3')
!       call diag_integral_field_init ('abs_sw', 'f8.3')

!-----------------------------------------------------------------------
!------------ initialize diagnostic fields ---------------

    id_olr = &
      register_diag_field(mod_name, 'olr', axes(1:2), Time, &
                          'Outgoing longwave radiation at TOA', &
                          'W/m2', missing_value=missing_value)
    id_swdn_sfc = &
      register_diag_field(mod_name, 'swnet_sfc', axes(1:2), Time, &
                          'Net SW flux at surface (positive down)', &
                          'W/m2', missing_value=missing_value)
    id_swdn_toa = &
      register_diag_field(mod_name, 'swdn_toa', axes(1:2), Time, &
                          'SW flux down at TOA', &
                          'W/m2', missing_value=missing_value)
    id_swnet_toa = &
      register_diag_field(mod_name, 'swnet_toa', axes(1:2), Time, &
                          'Net SW flux at TOA (positive down)', &
                          'W/m2', missing_value=missing_value)
    id_lwup_sfc = &
      register_diag_field(mod_name, 'lwup_sfc', axes(1:2), Time, &
                          'LW flux up at surface', &
                          'W/m2', missing_value=missing_value)
    id_albedo = &
      register_diag_field(mod_name, 'albedo_rad', axes(1:2), Time, &
                          'Surface albedo seen by the radiation', &
                          '1', missing_value=missing_value)

    id_lwdn_sfc = &
      register_diag_field(mod_name, 'lwdn_sfc', axes(1:2), Time, &
                          'LW flux down at surface', &
                          'W/m2', missing_value=missing_value)

    id_tdt_rad = &
      register_diag_field(mod_name, 'tdt_rad', axes(1:3), Time, &
                          'Temperature tendency due to radiation', &
                          'K/s', missing_value=missing_value)

    id_flux_rad = &
      register_diag_field(mod_name, 'netrad_half', axes(half), Time, &
                          'Net radiative flux on half levels (positive up)', &
                          'W/m2', missing_value=missing_value)
    id_flux_lw = &
      register_diag_field(mod_name, 'lwnet_half', axes(half), Time, &
                          'Net LW flux on half levels (positive up)', &
                          'W/m2', missing_value=missing_value)
    id_flux_sw = &
      register_diag_field(mod_name, 'swnet_half', axes(half), Time, &
                          'Net SW flux on half levels (positive up)', &
                          'W/m2', missing_value=missing_value)
    id_entrop_rad = &
      register_diag_field(mod_name, 'entrop_rad', axes(1:3), Time, &
                          'Entropy production by radiation', &
                          '1/s', missing_value=missing_value)
    id_tau = &
      register_diag_field(mod_name, 'tau_lw', axes(half), Time, &
                          'LW optical depth on half levels', &
                          '1', missing_value=missing_value)
    id_tau_rad = &
      register_diag_field(mod_name, 'tau_sw', axes(half), Time, &
                          'SW optical depth on half levels', &
                          '1', missing_value=missing_value)

    return
  end subroutine gray_radiation_init

! ==================================================================================

  !> Computes the gray radiative fluxes and heating, adds the heating to `tdt`, returns the
  !> surface fluxes and sends the diagnostics.
  subroutine gray_radiation(is, js, Time_diag, lat, lon, p_half, albedo, t_surf, t, tdt, net_surf_sw_down, surf_lw_down)

    integer, intent(in)                 :: is, js  !! indices of the first point of the physics window in the processor domain
    type(time_type), intent(in)         :: Time_diag  !! time at which the diagnostics are sent
    real, intent(in), dimension(:, :)   :: lat, lon, albedo
    !! `lat`, `lon`: latitudes and longitudes [rad]; `albedo`: surface albedo
    real, intent(in), dimension(:, :)   :: t_surf  !! surface temperature [K]
    real, intent(in), dimension(:, :, :) :: t, p_half  !! `t`: temperature [K]; `p_half`: pressure at half levels [Pa]
    real, intent(inout), dimension(:, :, :) :: tdt  !! temperature tendency, to which the radiative heating is added [K/s]
    real, intent(out), dimension(:, :)   :: net_surf_sw_down, surf_lw_down
    !! net downward shortwave and downward longwave flux at the surface [W/m2]

    real, dimension(size(t, 2)) :: ss, ss2, ss4, ss6, ss8, solar, tau_0, solar_tau_0, p2
    real, dimension(size(t, 1), size(t, 2))              :: b_surf
    real, dimension(size(t, 1), size(t, 2), size(t, 3))   :: b, tdt_rad, entrop_rad
    real, dimension(size(t, 1), size(t, 2), size(t, 3) + 1) :: up, down, net, solar_down, flux_rad, flux_sw
    real, dimension(size(t, 1), size(t, 2), size(t, 3) + 1)   :: tau, solar_tau, dtrans
    real, dimension(size(t, 1), size(t, 2))   :: long_forcing, olr, swin

    real, dimension(size(t, 1), size(t, 2))              :: walker_forcing

    integer :: i, j, k, n

    real :: surf_albedo
    real :: dist
    logical :: used

    n = size(t, 3)

    do i = 1, size(t, 1)
      do j = 1, size(t, 2)
        long_forcing(i, j) = size_pert*exp(-(rad_to_deg*lon(i, j) - long_pert)**2./ &
                                           del_long**2.)
      end do
    end do

    do i = 1, size(t, 1)
      do j = 1, size(t, 2)
        dist = sqrt(((rad_to_deg*lat(i, j) - lat_pert)/del_lat)**2 &
                    + ((rad_to_deg*lon(i, j) - lon_pert)/del_lon)**2)
        if (dist .lt. 1) then
          walker_forcing(i, j) = fcng_pert*cos(pi*dist*0.5)**2; 
        else
          walker_forcing(i, j) = 0.0
        end if
      end do
    end do

    ss = sin(lat(1, :))
    ss2 = ss*ss
    p2 = (1.-3.*ss*ss)/4.

    solar = 0.25*solar_constant*(1.0 + del_sol*p2 + del_sw*ss)

    tau_0 = ir_tau_eq + (ir_tau_pole - ir_tau_eq)*ss*ss

    solar_tau_0 = (1.0 - sw_diff*ss*ss)*atm_abs

    b = stefan*t*t*t*t
    b_surf = stefan*t_surf*t_surf*t_surf*t_surf

    do k = 1, n + 1

!   tau(:,k)       = tau_0(:)      *(p_half(1,:,k)/p_half(1,:,n+1))**2

! modif omp: optical thickness is obtained by combination of linear and quadratic function.
! linear_tau controls the fraction of the optical depth that varies linearily. The
! surface pressure has been replaced by a reference value (p00 = 1000mb) from constant_mod.

!  tau(:,k)       = tau_0(:) * (linear_tau * p_half(1,:,k)/p00 + (1.0 - linear_tau) &
!       * (p_half(1,:,k)/p00)**2)
! modif df 1-23-04: changing profile of IR absorber
      do i = 1, size(t, 1)
        tau(i, :, k) = tau_0(:)*(linear_tau*p_half(i, :, k)/p00 + (1.0 - linear_tau) &
                                 *(p_half(i, :, k)/p00)**4)
!     solar_tau(i,:,k) = solar_tau_0(:)*(p_half(i,:,k)/p_half(i,:,n+1))**4
        solar_tau(i, :, k) = solar_tau_0(:)*(p_half(i, :, k)/p00)**4
      end do
    end do

    do k = 1, n
      do i = 1, size(t, 1)
        dtrans(i, :, k) = exp(-(tau(i, :, k + 1) - tau(i, :, k)))
      end do
    end do

    up(:, :, n + 1) = b_surf
    do k = n, 1, -1
      do j = 1, size(t, 2)
        do i = 1, size(t, 1)
          up(i, j, k) = up(i, j, k + 1)*dtrans(i, j, k) + b(i, j, k)*(1.0 - dtrans(i, j, k))
        end do
      end do
    end do

    down(:, :, 1) = 0.0
    do k = 1, n
      do j = 1, size(t, 2)
        do i = 1, size(t, 1)
          down(i, j, k + 1) = down(i, j, k)*dtrans(i, j, k) + b(i, j, k)*(1.0 - dtrans(i, j, k))
        end do
      end do
    end do

    do k = 1, n + 1
      do j = 1, size(t, 2)
        do i = 1, size(t, 1)
!    solar_down(i,j,k) = (solar(j)+long_forcing(i,j))*exp(-solar_tau(j,k))
          solar_down(i, j, k) = (walker_forcing(i, j) + long_forcing(i, j) &
                                 + solar(j))*exp(-solar_tau(i, j, k))
        end do
      end do
    end do

    do k = 1, n + 1
      do i = 1, size(t, 1)
        net(i, :, k) = up(i, :, k) - down(i, :, k)
        flux_sw(i, :, k) = albedo(i, :)*solar_down(i, :, n + 1) - solar_down(i, :, k)
        flux_rad(i, :, k) = net(i, :, k) + flux_sw(i, :, k)
      end do
    end do

    do k = 1, n
      tdt_rad(:, :, k) = (net(:, :, k + 1) - net(:, :, k) - solar_down(:, :, k + 1) + solar_down(:, :, k)) &
                         *grav/(cp_air*(p_half(:, :, k + 1) - p_half(:, :, k)))
      tdt(:, :, k) = tdt(:, :, k) + tdt_rad(:, :, k)
    end do

    surf_lw_down = down(:, :, n + 1)
    net_surf_sw_down = solar_down(:, :, n + 1)*(1.-albedo(:, :))
    olr = up(:, :, 1)
    swin = solar_down(:, :, 1)

!------- outgoing lw flux toa (olr) -------
    if (id_olr > 0) then
      used = send_data(id_olr, olr, Time_diag, is, js)
    end if
!------- downward sw flux surface -------
    if (id_swdn_sfc > 0) then
      used = send_data(id_swdn_sfc, net_surf_sw_down, Time_diag, is, js)
    end if
!------- incoming sw flux toa -------
    if (id_swdn_toa > 0) then
      used = send_data(id_swdn_toa, swin, Time_diag, is, js)
    end if
!------- net sw flux toa (positive down): no SW absorption above the surface on the way up -------
    if (id_swnet_toa > 0) then
      used = send_data(id_swnet_toa, -flux_sw(:, :, 1), Time_diag, is, js)
    end if
!------- upward lw flux surface -------
    if (id_lwup_sfc > 0) then
      used = send_data(id_lwup_sfc, b_surf, Time_diag, is, js)
    end if

!------- surface albedo -------
    if (id_albedo > 0) then
      used = send_data(id_albedo, albedo, Time_diag, is, js)
    end if
!------- downward lw flux surface -------
    if (id_lwdn_sfc > 0) then
      used = send_data(id_lwdn_sfc, surf_lw_down, Time_diag, is, js)
    end if
!------- temperature tendency due to radiation ------------
    if (id_tdt_rad > 0) then
      used = send_data(id_tdt_rad, tdt_rad, Time_diag, is, js, 1)
    end if
!------- total radiative flux (at half levels) -----------
    if (id_flux_rad > 0) then
      used = send_data(id_flux_rad, flux_rad, Time_diag, is, js, 1)
    end if
!------- longwave radiative flux (at half levels) --------
    if (id_flux_lw > 0) then
      used = send_data(id_flux_lw, net, Time_diag, is, js, 1)
    end if
    if (id_flux_sw > 0) then
      used = send_data(id_flux_sw, flux_sw, Time_diag, is, js, 1)
    end if
    if (id_entrop_rad > 0) then
      do k = 1, n
        entrop_rad(:, :, k) = tdt_rad(:, :, k)/t(:, :, k)*p_half(:, :, n + 1)/1.e5
      end do
      used = send_data(id_entrop_rad, entrop_rad, Time_diag, is, js, 1)
    end if
!mj debug
    if (id_tau > 0) then
      used = send_data(id_tau, tau, Time_diag, is, js, 1)
    end if
    if (id_tau_rad > 0) then
      used = send_data(id_tau_rad, solar_tau, Time_diag, is, js, 1)
    end if

    return
  end subroutine gray_radiation

! ==================================================================================

  !> Finalizes the module (nothing to do).
  subroutine gray_radiation_end()
  end subroutine gray_radiation_end

! ==================================================================================

end module gray_radiation_mod
