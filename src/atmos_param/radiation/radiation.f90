module radiation_mod

!-----------------------------------------------------------------------
!
!   Radiation driver: selects the radiation scheme and forwards the
!   physics_driver calls to it.
!
!   radiation_nml:
!     radiation_scheme = 'rrtm'  RRTMG clear-sky radiation (default)
!                        'gray'  gray radiation (Frierson et al. 2006)
!                        'none'  no radiative heating or surface fluxes
!
!-----------------------------------------------------------------------

use fms_mod,            only: input_nml_file, check_nml_error, &
                              mpp_pe, mpp_root_pe, stdlog,       &
                              error_mesg, FATAL, write_version_number
use time_manager_mod,   only: time_type
use constants_mod,      only: cp_air
use gray_radiation_mod, only: gray_radiation_init, gray_radiation, gray_radiation_end
use rrtmg_lw_init,      only: rrtmg_lw_ini
use rrtmg_sw_init,      only: rrtmg_sw_ini
use rrtm_radiation,     only: rrtm_radiation_init, interp_temp, run_rrtmg, &
                              rrtm_radiation_end, rrtm_precip_accum

implicit none
private

public :: radiation_init, radiation_down, radiation_precip_accum, radiation_end

character(len=128) :: version = '$Id: radiation.f90 $'
character(len=128) :: tagname = '$Name: $'

!-------------------- namelist -----------------------------------------

character(len=16) :: radiation_scheme = 'rrtm'   ! 'rrtm', 'gray' or 'none'

namelist /radiation_nml/ radiation_scheme

logical :: module_is_initialized = .false.

contains

!#######################################################################

subroutine radiation_init(axes, Time, id, jd, kd, lonb, latb)

integer,         intent(in), dimension(4) :: axes
type(time_type), intent(in)               :: Time
integer,         intent(in)               :: id, jd, kd
real,            intent(in), dimension(:) :: lonb, latb

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
   call rrtm_radiation_init(axes, Time, id*jd, kd, lonb, latb)
case ('none')
case default
   call error_mesg('radiation_init', 'radiation_scheme must be ''rrtm'', ''gray'' or ''none'', not '// &
                   trim(radiation_scheme), FATAL)
end select

module_is_initialized = .true.

end subroutine radiation_init

!#######################################################################

subroutine radiation_down(is, js, Time, Time_next, lat, lon, p_full, p_half, z_full, z_half, &
                          t, q, t_surf_rad, albedo, tdt, flux_sw, flux_lw)

!-----------------------------------------------------------------------
!   Adds the radiative heating to tdt and returns the net shortwave and
!   downward longwave fluxes at the surface. flux_sw and flux_lw must be
!   set by the caller beforehand; they are left unchanged by
!   radiation_scheme = 'none'.
!-----------------------------------------------------------------------

integer,         intent(in)                     :: is, js
type(time_type), intent(in)                     :: Time, Time_next
real,            intent(in),    dimension(:,:)  :: lat, lon, t_surf_rad, albedo
real,            intent(in),    dimension(:,:,:):: p_full, p_half, z_full, z_half, t, q
real,            intent(inout), dimension(:,:,:):: tdt
real,            intent(inout), dimension(:,:)  :: flux_sw, flux_lw

real, dimension(size(albedo,1), size(albedo,2)) :: coszen

select case (trim(radiation_scheme))
case ('gray')
   call gray_radiation(is, js, Time_next, lat, lon, p_half, albedo, t_surf_rad, t, tdt, flux_sw, flux_lw)
case ('rrtm')
   ! RRTM needs the temperature at half levels
   call interp_temp(z_full, z_half, t_surf_rad, t)
   call run_rrtmg(is, js, Time, lat, lon, p_full, p_half, albedo, q, t, t_surf_rad, tdt, coszen, flux_sw, flux_lw)
end select

end subroutine radiation_down

!#######################################################################

subroutine radiation_precip_accum(precip, rain, snow)

!-----------------------------------------------------------------------
!   Records where it precipitated on this step, for RRTM's optional
!   precipitation-dependent albedo (rrtm_radiation_nml do_precip_albedo).
!   precip is the total precipitation; rain and snow are the large-scale
!   parts.
!-----------------------------------------------------------------------

real, intent(in), dimension(:,:) :: precip, rain, snow

if (trim(radiation_scheme) == 'rrtm') call rrtm_precip_accum(precip, rain, snow)

end subroutine radiation_precip_accum

!#######################################################################

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
