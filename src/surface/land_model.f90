!> The land surface: where the land is, its roughness, heat capacity and surface
!> temperature.
!>
!> The land points are set by `land_mask`: none, the Navy land-sea mask (the `water_file` of
!> `topography_nml`), a land-sea mask file, longitude-latitude boxes, or the points higher
!> than `land_height_limit`. Over land the momentum and moisture roughness lengths are
!> multiplied by `rough_mom_factor` and `rough_moist_factor` (the moisture roughness can also
!> vary with latitude), which throttles the evaporation of the otherwise saturated surface.
!> The land surface temperature `t_surf` is
!>
!> * `'simple'`: a slab with heat capacity `heat_capacity`, stepped with the implicit surface
!>   energy budget (`slab_tendency_mod`);
!> * `'prescribed'`: read from `t_surf_file` at every step;
!> * `'fixed'`: held at its initial value.
!>
!> In every mode the local surface heating of `local_heating_nml` is added. On a cold start
!> the land takes the ocean's cold-start SST (`t_init = 'ocean'`) or `t_surf_file`. The state
!> is saved in `RESTART/land_model.res.nc`, which takes precedence when it exists in `INPUT/`.
!>
!> Namelist: `land_model_nml`
!> ([namelist reference](https://eddy-stanford.github.io/MiMA/Parameters/#land_model_nml)).
module land_model_mod

  use fms_mod, only: input_nml_file, check_nml_error, error_mesg, FATAL, NOTE, &
                     mpp_pe, mpp_root_pe, stdlog
  use fms2_io_mod, only: FmsNetcdfFile_t, open_file, close_file, read_data
  use restart_file_mod, only: restart_file_type, open_restart_read, open_restart_write, &
                              close_restart, read_restart_field, write_restart_field, &
                              restart_field_exists, check_field_size
  use diag_manager_mod, only: register_static_field, send_data
  use time_manager_mod, only: time_type
  use constants_mod, only: hlf, grav
  use mpp_domains_mod, only: mpp_get_compute_domain, mpp_get_global_domain
! mj know about surface topography
  use spectral_dynamics_mod, only: get_surf_geopotential
  use topography_mod, only: get_ocean_mask
  use mima_interpolator_mod, only: interpolate_type, interpolator_init, CONSTANT, interpolator
  use surface_types_mod, only: surface_grid_type, surface_exchange_type
  use slab_tendency_mod, only: slab_dt

  implicit none
  private

  public :: land_data_type, land_model_read_nml_and_mask, land_model_init_state, &
            land_model_roughness, update_land_model_fast, land_model_end

  !> State of the land surface. All arrays are on this PE's grid; off the land mask they
  !> are 0.
  type land_data_type
    logical, allocatable :: mask(:, :)  !! land points
    real, allocatable :: t_surf(:, :)  !! surface temperature [K]
    real, allocatable :: cap(:, :)  !! heat capacity [J/m2/K]
    real, allocatable :: dt_t_surf(:, :)  !! change of `t_surf` over the last step [K]
  end type land_data_type

  character(len=*), parameter :: mod_name = 'land'

  character(len=32) :: land_model = 'simple'
  !! `'simple'`: slab land surface with an interactive temperature; `'prescribed'`: temperature
  !! read from `t_surf_file` at every step; `'fixed'`: temperature held at its initial value
  character(len=32) :: land_mask = 'navy'
  !! where the land is: `'none'`; `'navy'`: Navy land-sea mask (the `water_file` of
  !! `topography_nml`); `'file'`: land where the variable `land_mask_file` of
  !! `INPUT/<land_mask_file>.nc` (on the model grid) is positive; `'boxes'`: the boxes
  !! `slandlon`..`elandlon`, `slandlat`..`elandlat`; `'height'`: surface height above
  !! `land_height_limit`
  character(len=256) :: land_mask_file = 'lmask'
  !! with `land_mask = 'file'`, name of the land-sea mask file (without `.nc`, in `INPUT/`) and
  !! of its variable
  real, dimension(10) :: slandlon = 0, slandlat = 0, elandlon = -1, elandlat = -1
  !! with `land_mask = 'boxes'`, start and end longitude and latitude [deg] of up to 10
  !! land boxes
  real :: land_height_limit = 10.  !! [m] with `land_mask = 'height'`, points higher than this are land
  real :: heat_capacity = 1.e07  !! [J/m2/K] heat capacity of the land; <= 0: the ocean `heat_capacity`
  character(len=32) :: t_init = 'ocean'
  !! initial land temperature on a cold start: `'ocean'`: the cold-start SST of
  !! `ocean_model_nml` (`sst_init`); `'file'`: from `t_surf_file`. `land_model = 'prescribed'`
  !! always starts from `t_surf_file`
  character(len=256) :: t_surf_file = ''  !! land temperature file name, without `.nc`, in `INPUT/`
  character(len=256) :: t_surf_file_var = ''
  !! name of the temperature variable in `t_surf_file` (empty: `t_surf_file`)
  real :: rough_mom_factor = 5.e3  !! factor on the ocean roughness `rough_ocean` for the land momentum roughness
  real :: rough_moist_factor = 1.e-12  !! factor on `rough_ocean` for the land moisture roughness
  character(len=32) :: rough_moist_profile = 'legacy_tropical'
  !! land moisture roughness: `'uniform'`: `rough_ocean*rough_moist_factor`; `'legacy_tropical'`:
  !! as `'uniform'`, plus a larger value over tropical land (legacy)
  character(len=32) :: hydrology = 'saturated'
  !! land hydrology: `'saturated'`: the land surface is saturated, as the ocean (`'bucket'` is
  !! not yet available)
  logical :: do_snow = .false.  !! snow on land (not yet available; needs `hydrology = 'bucket'`)

  namelist /land_model_nml/ land_model, land_mask, land_mask_file, &
    slandlon, slandlat, elandlon, elandlat, land_height_limit, &
    heat_capacity, t_init, t_surf_file, t_surf_file_var, &
    rough_mom_factor, rough_moist_factor, rough_moist_profile, &
    hydrology, do_snow

  type(interpolate_type), save :: t_interp
  character(len=256) :: t_var
  integer :: id_land_mask

contains

!#######################################################################

  !> Reads `land_model_nml` and sets the land mask.
  subroutine land_model_read_nml_and_mask(Land, grid)
    type(land_data_type), intent(inout) :: Land  !! the land state
    type(surface_grid_type), intent(in) :: grid  !! the surface grid

    integer :: ierr, io, i, j, k, is, ie, js, je, nlon, nlat
    real :: lat, lon, pi
    logical :: ocean_mask_worked
    logical, dimension(grid%nx, grid%ny) :: lmask_navy
    real, dimension(grid%nx, grid%ny) :: zsurf
    real, allocatable, dimension(:, :) :: global_mask
    type(FmsNetcdfFile_t) :: mask_file

    pi = 4.0*atan(1.)

    read (input_nml_file, nml=land_model_nml, iostat=io)
    ierr = check_nml_error(io, 'land_model_nml')

    select case (trim(land_model))
    case ('simple', 'fixed')
    case ('prescribed')
      if (len_trim(t_surf_file) == 0) &
        call error_mesg('land_model_init', "land_model = 'prescribed' needs t_surf_file", FATAL)
    case default
      call error_mesg('land_model_init', 'unknown land_model = "'//trim(land_model)// &
                      '"; use simple, prescribed or fixed', FATAL)
    end select
    select case (trim(t_init))
    case ('ocean')
    case ('file')
      if (len_trim(t_surf_file) == 0) &
        call error_mesg('land_model_init', "t_init = 'file' needs t_surf_file", FATAL)
    case default
      call error_mesg('land_model_init', 'unknown t_init = "'//trim(t_init)//'"; use ocean or file', FATAL)
    end select
    select case (trim(rough_moist_profile))
    case ('uniform', 'legacy_tropical')
    case default
      call error_mesg('land_model_init', 'unknown rough_moist_profile = "'//trim(rough_moist_profile)// &
                      '"; use uniform or legacy_tropical', FATAL)
    end select
    select case (trim(hydrology))
    case ('saturated')
    case ('bucket')
      if (rough_moist_factor /= 1. .or. trim(rough_moist_profile) /= 'uniform') &
        call error_mesg('land_model_init', "hydrology = 'bucket' must not be combined with the legacy "// &
                        "evaporation throttle: set rough_moist_factor = 1. and rough_moist_profile = 'uniform'", FATAL)
      call error_mesg('land_model_init', "hydrology = 'bucket' is not implemented yet "// &
                      '(available from stage S5b)', FATAL)
    case default
      call error_mesg('land_model_init', 'unknown hydrology = "'//trim(hydrology)// &
                      '"; use saturated or bucket', FATAL)
    end select
    if (do_snow) then
      if (trim(hydrology) /= 'bucket') &
        call error_mesg('land_model_init', "do_snow needs hydrology = 'bucket'", FATAL)
      call error_mesg('land_model_init', 'do_snow is not implemented yet (available from stage S5b)', FATAL)
    end if

    allocate (Land%mask(grid%nx, grid%ny))
    Land%mask = .false.

! Each branch records the land points where it sets them
    select case (trim(land_mask))
    case ('none')
      if (trim(land_model) /= 'simple') call error_mesg('land_model_init', &
                                                        "land_mask = 'none': there are no land points, "// &
                                                        'so land_model = "'//trim(land_model)//'" has no effect', NOTE)
    case ('file')
      if (mpp_pe() .eq. mpp_root_pe()) write (*, '(a)') 'Reading land-sea mask from file INPUT/'//trim(land_mask_file)//'.nc'
      if (.not. open_file(mask_file, 'INPUT/'//trim(land_mask_file)//'.nc', 'read')) &
        call error_mesg('land_model_init', &
                        'cannot open INPUT/'//trim(land_mask_file)//'.nc', FATAL)
      call mpp_get_compute_domain(grid%domain, is, ie, js, je)
      call mpp_get_global_domain(grid%domain, xsize=nlon, ysize=nlat)
      call check_field_size(mask_file, trim(land_mask_file), (/nlon, nlat/))
      allocate (global_mask(nlon, nlat))
      call read_data(mask_file, trim(land_mask_file), global_mask)
      call close_file(mask_file)
      where (global_mask(is:ie, js:je) .gt. 0) Land%mask = .true.
      deallocate (global_mask)
! mj use navy land-sea mask
    case ('navy')
      ocean_mask_worked = get_ocean_mask(grid%lon_bnd, grid%lat_bnd, lmask_navy)
      if (.not. ocean_mask_worked) then
        call error_mesg('get_ocean_mask', 'land_mask="'//trim(land_mask)//'"'// &
                        ' and ocean_mask is not present but water data file does not exist', FATAL)
      end if
      where (.not. lmask_navy) Land%mask = .true.
! mj land function of surface topography
    case ('height')
      call get_surf_geopotential(zsurf)
      where (zsurf/grav > land_height_limit) Land%mask = .true.
! mj land given through ?landlon, ?landlat
    case ('boxes')
      do j = 1, grid%ny
        lat = 0.5*180/pi*(grid%lat_bnd(j + 1) + grid%lat_bnd(j))
        do i = 1, grid%nx
          lon = 0.5*180/pi*(grid%lon_bnd(i + 1) + grid%lon_bnd(i))
          do k = 1, size(slandlat)
            if (lon >= slandlon(k) .and. lon <= elandlon(k) &
                 &.and. lat >= slandlat(k) .and. lat <= elandlat(k)) then
              Land%mask(i, j) = .true.
            end if
          end do
        end do
      end do
    case default
      call error_mesg('land_model_init', 'unknown land_mask = "'//trim(land_mask)// &
                      '"; use none, navy, file, boxes or height', FATAL)
    end select

    if (mpp_pe() == mpp_root_pe()) write (stdlog(), nml=land_model_nml)

  end subroutine land_model_read_nml_and_mask

!#######################################################################

  !> Sets the land heat capacity and the initial land temperature (from
  !> `INPUT/land_model.res.nc` if it exists), and registers the diagnostics.
  subroutine land_model_init_state(Land, Time, grid, t_all, heat_capacity_ocean_default)
    type(land_data_type), intent(inout) :: Land  !! the land state
    type(time_type), intent(in) :: Time  !! current time
    type(surface_grid_type), intent(in) :: grid  !! the surface grid
    real, intent(in) :: t_all(:, :)  !! the ocean's cold-start SST on all points [K]
    real, intent(in) :: heat_capacity_ocean_default  !! the ocean `heat_capacity` [J/m2/K]

    real, dimension(grid%nx, grid%ny) :: t_file, t_rst
    type(restart_file_type) :: rst
    logical :: used

    allocate (Land%t_surf(grid%nx, grid%ny), Land%cap(grid%nx, grid%ny), &
              Land%dt_t_surf(grid%nx, grid%ny))
    Land%dt_t_surf = 0.0

    if (heat_capacity .le. 0.) heat_capacity = heat_capacity_ocean_default
    Land%cap = 0.0
    where (Land%mask) Land%cap = heat_capacity

    t_var = t_surf_file_var
    if (len_trim(t_var) == 0) t_var = t_surf_file
    if (trim(land_model) == 'prescribed' .or. trim(t_init) == 'file') &
      call interpolator_init(t_interp, trim(t_surf_file)//'.nc', grid%lon_bnd, grid%lat_bnd, &
                             data_names=(/t_var/), data_out_of_bounds=(/CONSTANT/))

    Land%t_surf = 0.0
    if (trim(land_model) == 'prescribed' .or. trim(t_init) == 'file') then
      call interpolator(t_interp, Time, t_file, trim(t_var))
      where (Land%mask) Land%t_surf = t_file
    else
      where (Land%mask) Land%t_surf = t_all
    end if

! The restart file takes precedence over t_init
    if (open_restart_read(rst, 'INPUT/land_model.res.nc', grid%domain)) then
      if (.not. restart_field_exists(rst, 't_surf')) &
        call error_mesg('land_model_init', 'INPUT/land_model.res.nc has no field t_surf', FATAL)
      call read_restart_field(rst, 't_surf', t_rst)
      call close_restart(rst)
      Land%t_surf = 0.0
      where (Land%mask) Land%t_surf = t_rst
    end if

    id_land_mask = register_static_field(mod_name, 'land_mask', grid%axes, &
                                         'land mask (1: land, 0: ocean)', '1')
    if (id_land_mask > 0) used = send_data(id_land_mask, merge(1., 0., Land%mask))

  end subroutine land_model_init_state

!#######################################################################

  !> Sets the roughness lengths on the land points.
  !>
  !> Momentum: `rough_ocean*rough_mom_factor`; heat: `rough_ocean`; moisture: according to
  !> `rough_moist_profile`. The arrays are left unchanged elsewhere.
  subroutine land_model_roughness(Land, grid, rough_ocean, rough_mom, rough_heat, rough_moist)
    type(land_data_type), intent(in) :: Land  !! the land state
    type(surface_grid_type), intent(in) :: grid  !! the surface grid
    real, intent(in) :: rough_ocean  !! roughness length of the ocean [m]
    real, intent(inout) :: rough_mom(:, :), rough_heat(:, :), rough_moist(:, :)
    !! roughness lengths for momentum, heat and moisture [m]

    integer :: j
    real :: lat, pi

    pi = 4.0*atan(1.)

    !cig: set higher roughness values over land as compared to ocean
    where (Land%mask) rough_mom = rough_ocean*rough_mom_factor
    where (Land%mask) rough_heat = rough_ocean
    select case (trim(rough_moist_profile))
    case ('uniform')
      where (Land%mask) rough_moist = rough_ocean*rough_moist_factor
    case ('legacy_tropical')
      !cig: more evaporation over tropics and midlatitudes as compared to subtropics
      do j = 1, size(rough_moist, 2)
        lat = 0.5*(grid%lat_bnd(j + 1) + grid%lat_bnd(j))*180./pi
        where (Land%mask(:, j)) rough_moist(:, j) = rough_ocean*rough_moist_factor + &
                    &         +(1.e-7)*exp(-abs(lat - 0.)**3./(2*15.)) &
                    &         + (1.e-25)*exp(-abs(lat - 45.)**3./(2*30.)) &
                    &         + (1.e-25)*exp(-abs(lat + 45.)**3./(2*30.))
      end do
    end select

  end subroutine land_model_roughness

!#######################################################################

  !> Steps the land surface temperature over one atmospheric time step.
  !>
  !> `'simple'`: implicit surface energy budget; `'prescribed'`: change to the temperature of
  !> `t_surf_file`; `'fixed'`: no change. The local surface heating is then added in every
  !> mode.
  subroutine update_land_model_fast(Exch, Land)
    type(surface_exchange_type), intent(in) :: Exch  !! fluxes and derivatives from the flux exchange
    type(land_data_type), intent(inout) :: Land  !! the land state

    real, dimension(size(Land%t_surf, 1), size(Land%t_surf, 2)) :: t_new

    select case (trim(land_model))
    case ('simple')
      where (Land%mask)
        Land%dt_t_surf = slab_dt(Exch%flux_lw, Exch%flux_sw, hlf*Exch%fprec, Exch%flux_t, &
                                 Exch%flux_q, Exch%qflux_legacy, Exch%dhdt_surf, Exch%dedt_surf, &
                                 Exch%drdt_surf, Exch%dt, Land%cap)
      elsewhere
        Land%dt_t_surf = 0.0
      end where
    case ('prescribed')
      call interpolator(t_interp, Exch%Time, t_new, trim(t_var))
      where (Land%mask)
        Land%dt_t_surf = t_new - Land%t_surf
      elsewhere
        Land%dt_t_surf = 0.0
      end where
    case ('fixed')
      Land%dt_t_surf = 0.0
    end select

    where (Land%mask)
      Land%dt_t_surf = Land%dt_t_surf + Exch%extra_heating*Exch%dt
      Land%t_surf = Land%t_surf + Land%dt_t_surf
    end where

  end subroutine update_land_model_fast

!#######################################################################

  !> Writes the restart file `RESTART/land_model.res.nc` (the surface temperature, 0 over
  !> the ocean).
  subroutine land_model_end(Land, grid)
    type(land_data_type), intent(in) :: Land  !! the land state
    type(surface_grid_type), intent(in) :: grid  !! the surface grid
    type(restart_file_type) :: rst

    call open_restart_write(rst, 'RESTART/land_model.res.nc', grid%domain)
    call write_restart_field(rst, 't_surf', Land%t_surf)
    call close_restart(rst)

  end subroutine land_model_end

end module land_model_mod
