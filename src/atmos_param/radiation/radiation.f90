!> Radiation driver: selects the radiation scheme and forwards the physics driver calls to it.
!>
!> `radiation_scheme` chooses between RRTMG clear-sky radiation (`rrtm_radiation`, configured
!> with `rrtm_radiation_nml` and `astro_nml`), the gray radiation of Frierson et al. (2006)
!> (`gray_radiation_mod`, configured with `gray_radiation_nml`), and no radiation at all.
!> See [radiation options](https://eddy-stanford.github.io/MiMA/Configurations/#radiation-options).
!>
!> Namelist: `radiation_nml`
!> ([namelist reference](https://eddy-stanford.github.io/MiMA/Parameters/#radiation_nml)).
!>
!> References:
!>
!> * Frierson, D. M. W., I. M. Held, and P. Zurita-Gotor, 2006: A gray-radiation aquaplanet
!>   moist GCM. Part I: Static stability and eddy scale. J. Atmos. Sci., 63, 2548-2566,
!>   https://doi.org/10.1175/JAS3753.1.
module radiation_mod

  use fms_mod, only: input_nml_file, check_nml_error, &
                     mpp_pe, mpp_root_pe, stdlog, &
                     error_mesg, FATAL, write_version_number
  use time_manager_mod, only: time_type
  use mpp_domains_mod, only: domain2d
  use constants_mod, only: cp_air
  use gray_radiation_mod, only: gray_radiation_init, gray_radiation, gray_radiation_end
  use rrtmg_lw_init, only: rrtmg_lw_ini
  use rrtmg_sw_init, only: rrtmg_sw_ini
  use rrtm_radiation, only: rrtm_radiation_init, interp_temp, run_rrtmg, &
                            rrtm_radiation_end, rrtm_precip_accum

  implicit none
  private

  public :: radiation_init, radiation_down, radiation_precip_accum, radiation_end

  character(len=128) :: version = '$Id: radiation.f90 $'
  character(len=128) :: tagname = '$Name: $'

!-------------------- namelist -----------------------------------------

  character(len=16) :: radiation_scheme = 'rrtm'   !! `'rrtm'`: RRTMG clear-sky radiation
  !! (`rrtm_radiation_nml`, `astro_nml`); `'gray'`: gray radiation (`gray_radiation_nml`); `'none'`:
  !! no radiative heating and no radiative surface fluxes. See
  !! [radiation options](Configurations.md#radiation-options).

  namelist /radiation_nml/ radiation_scheme

  logical :: module_is_initialized = .false.

contains

!#######################################################################

  !> Initializes the module: reads `radiation_nml` and initializes the selected scheme.
  subroutine radiation_init(axes, Time, id, jd, kd, lonb, latb, domain)

    integer, intent(in), dimension(4) :: axes  !! diagnostic axes (lon, lat, pfull, phalf)
    type(time_type), intent(in)               :: Time  !! current time
    integer, intent(in)               :: id, jd, kd  !! numbers of longitudes, latitudes and levels on this processor
    real, intent(in), dimension(:) :: lonb, latb  !! longitudes and latitudes of the cell corners [rad]
    type(domain2d), intent(in)               :: domain   !! grid domain, for restart files

    integer :: unit, ierr, io

    if (module_is_initialized) return

    read (input_nml_file, nml=radiation_nml, iostat=io)
    ierr = check_nml_error(io, 'radiation_nml')

    call write_version_number(version, tagname)
    if (mpp_pe() == mpp_root_pe()) write (stdlog(), nml=radiation_nml)

    select case (trim(radiation_scheme))
    case ('gray')
      call gray_radiation_init(axes, Time)
    case ('rrtm')
      call rrtmg_lw_ini(cp_air)
      call rrtmg_sw_ini(cp_air)
      call rrtm_radiation_init(axes, Time, id*jd, kd, lonb, latb, domain)
    case ('none')
    case default
      call error_mesg('radiation_init', 'radiation_scheme must be ''rrtm'', ''gray'' or ''none'', not '// &
                      trim(radiation_scheme), FATAL)
    end select

    module_is_initialized = .true.

  end subroutine radiation_init

!#######################################################################

  !> Adds the radiative heating to `tdt` and returns the net shortwave and downward
  !> longwave fluxes at the surface.
  !>
  !> `flux_sw` and `flux_lw` must be set by the caller beforehand; they are left unchanged
  !> by `radiation_scheme = 'none'`.
  subroutine radiation_down(is, js, Time, Time_next, lat, lon, p_full, p_half, z_full, z_half, &
                            t, q, t_surf_rad, albedo, tdt, flux_sw, flux_lw)

    integer, intent(in)                     :: is, js  !! indices of the first point of the physics window in the processor domain
    type(time_type), intent(in)                     :: Time, Time_next
    !! `Time`: current time; `Time_next`: time at the end of the step, at which the diagnostics are sent
    real, intent(in), dimension(:, :)  :: lat, lon, t_surf_rad, albedo
    !! `lat`, `lon`: latitudes and longitudes [rad]; `t_surf_rad`: surface temperature for the
    !! radiation [K]; `albedo`: surface albedo
    real, intent(in), dimension(:, :, :):: p_full, p_half, z_full, z_half, t, q
    !! `p_full`, `p_half`: pressure at full and half levels [Pa]; `z_full`, `z_half`: height at
    !! full and half levels [m]; `t`: temperature [K]; `q`: specific humidity [kg/kg]
    real, intent(inout), dimension(:, :, :):: tdt  !! temperature tendency, to which the radiative heating is added [K/s]
    real, intent(inout), dimension(:, :)  :: flux_sw, flux_lw
    !! net downward shortwave and downward longwave flux at the surface [W/m2]

    real, dimension(size(albedo, 1), size(albedo, 2)) :: coszen

    select case (trim(radiation_scheme))
    case ('gray')
      call gray_radiation(is, js, Time_next, lat, lon, p_half, albedo, t_surf_rad, t, tdt, flux_sw, flux_lw)
    case ('rrtm')
      ! RRTM needs the temperature at half levels
      call interp_temp(z_full, z_half, t_surf_rad, t)
      ! RRTM computes at Time and, like the other physics, sends its diagnostics at Time_next
      call run_rrtmg(is, js, Time, Time_next, lat, lon, p_full, p_half, albedo, q, t, t_surf_rad, tdt, coszen, &
                     flux_sw, flux_lw)
    end select

  end subroutine radiation_down

!#######################################################################

  !> Records where it precipitated on this step, for the optional precipitation-dependent
  !> albedo of RRTM (`do_precip_albedo` in `rrtm_radiation_nml`).
  subroutine radiation_precip_accum(precip, rain, snow)

    real, intent(in), dimension(:, :) :: precip, rain, snow
    !! `precip`: total precipitation; `rain`, `snow`: large-scale rain and snow [kg/m2/s]

    if (trim(radiation_scheme) == 'rrtm') call rrtm_precip_accum(precip, rain, snow)

  end subroutine radiation_precip_accum

!#######################################################################

  !> Finalizes the selected scheme (RRTM writes its restart file).
  subroutine radiation_end

    select case (trim(radiation_scheme))
    case ('gray')
      call gray_radiation_end
    case ('rrtm')
      call rrtm_radiation_end
    end select

    module_is_initialized = .false.

  end subroutine radiation_end

!#######################################################################

end module radiation_mod
