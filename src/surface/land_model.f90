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
!> In every mode the local surface heating of `local_heating_nml` is added: for `'simple'` as
!> the heat flux `cap*extra_heating` in the implicit step, so that the fluxes returned to the
!> atmosphere are those of the final temperature (energy-consistent); otherwise as an
!> increment `extra_heating*dt` of the temperature. On a cold start
!> the land takes the ocean's cold-start SST (`t_init = 'ocean'`) or `t_surf_file`. The state
!> is saved in `RESTART/land_model.res.nc`, which takes precedence when it exists in `INPUT/`.
!>
!> The hydrology is `'saturated'` (the default: the land evaporates as a wet surface, throttled
!> only by the roughness) or a Manabe `'bucket'`. The bucket holds the soil water `W`
!> [kg/m2], filled by rain (and, without snow, snowfall) and emptied by evaporation; above
!> `water_capacity` the excess runs off. The evaporation is the potential (saturated)
!> evaporation times `beta = min(W/(beta_frac*water_capacity), 1)`. With `do_snow`, snowfall
!> is stored as snow `S` [kg/m2] instead of melting on arrival: snow sublimates first (at the
!> potential rate, `beta = 1` while `S > 0`), takes up frost, melts into the bucket when the
!> surface is above freezing, brightens the albedo, and above `snow_max` runs off (a crude
!> glacier). Sublimation and melt take their latent heat of fusion from the surface
!> temperature after the implicit step, and frost releases it. The scheme conserves water
!> exactly, apart from the diagnosed `water_clip`; it is verified for conservation over short
!> runs, not validated for climate (in particular the soil water takes years to spin up).
!> The buoyancy in `surface_flux_mod` still uses the saturated surface humidity, an
!> approximation. Note that MiMA's moist physics returns all precipitation as rain (the frozen
!> precipitation `fprec` is 0), so at present snow comes only from a restart or from frost.
!> `W` and `S` are saved in `RESTART/land_model.res.nc`; points that were not land when the
!> restart was written start dry.
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
  use diag_manager_mod, only: register_diag_field, register_static_field, send_data
  use time_manager_mod, only: time_type
  use constants_mod, only: hlf, hlv, grav, tfreeze
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
            land_model_read_restart, land_model_roughness, land_model_boundary, &
            update_land_model_fast, land_model_end

  !> State of the land surface. All arrays are on this PE's grid; off the land mask they
  !> are 0.
  type land_data_type
    logical, allocatable :: mask(:, :)  !! land points
    real, allocatable :: t_surf(:, :)  !! surface temperature [K]
    real, allocatable :: cap(:, :)  !! heat capacity [J/m2/K]
    real, allocatable :: dt_t_surf(:, :)
    !! change of `t_surf` over the last step from the implicit step and the local surface
    !! heating, before the snow adjustments (the change seen by the atmosphere) [K]
    logical :: do_bucket = .false.  !! the hydrology is `'bucket'`; the arrays below exist only then
    logical :: do_snow = .false.  !! snow is on (`do_snow`); the snow arrays exist only then
    real, allocatable :: water(:, :)  !! soil water `W` [kg/m2]
    real, allocatable :: beta(:, :)  !! evaporation efficiency (fraction of the potential evaporation)
    real, allocatable :: runoff(:, :)  !! runoff rate over the last step [kg/m2/s]
    real, allocatable :: water_clip(:, :)
    !! rate at which water was added over the last step to stop `W` going negative [kg/m2/s]
    real, allocatable :: snow(:, :)  !! snow `S` [kg/m2]
    real, allocatable :: melt(:, :)  !! snow melt rate over the last step [kg/m2/s]
    real, allocatable :: sublimation(:, :)  !! sublimation rate over the last step [kg/m2/s]
    real, allocatable :: frost(:, :)  !! frost (deposition) rate over the last step [kg/m2/s]
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
  !! land hydrology: `'saturated'`: the land surface is saturated, as the ocean; `'bucket'`:
  !! Manabe bucket, evaporation limited by the soil water (needs `rough_moist_factor = 1.` and
  !! `rough_moist_profile = 'uniform'`)
  real :: water_capacity = 150.
  !! [kg/m2] with `hydrology = 'bucket'`, field capacity of the bucket; more water runs off
  real :: beta_frac = 0.75
  !! with `hydrology = 'bucket'`, fraction of `water_capacity` above which the land evaporates
  !! at the potential rate (below it, in proportion to the water)
  real :: water_init_frac = 1.0
  !! with `hydrology = 'bucket'`, initial soil water on a cold start, as a fraction of
  !! `water_capacity`
  logical :: do_snow = .false.
  !! snow on land: stored snowfall, sublimation, frost, melt and snow albedo (needs
  !! `hydrology = 'bucket'` and `land_model = 'simple'`)
  real :: albedo_snow = 0.60
  !! with `do_snow`, albedo of deep snow (snow never lowers the albedo below the background)
  real :: snow_mask_depth = 10.
  !! [kg/m2] with `do_snow`, snow at which the albedo is halfway between the background and
  !! `albedo_snow`
  real :: snow_max = 1000.  !! [kg/m2] with `do_snow`, maximum snow; more runs off (a crude glacier)

  namelist /land_model_nml/ land_model, land_mask, land_mask_file, &
    slandlon, slandlat, elandlon, elandlat, land_height_limit, &
    heat_capacity, t_init, t_surf_file, t_surf_file_var, &
    rough_mom_factor, rough_moist_factor, rough_moist_profile, &
    hydrology, water_capacity, beta_frac, water_init_frac, &
    do_snow, albedo_snow, snow_mask_depth, snow_max

  type(interpolate_type), save :: t_interp
  character(len=256) :: t_var
  integer :: id_land_mask
  integer :: id_water = 0, id_beta = 0, id_runoff = 0, id_water_clip = 0, &
             id_snow = 0, id_melt = 0, id_sublimation = 0, id_frost = 0, &
             id_f_net_final = 0, id_extra_heating = 0

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
      if (water_capacity <= 0.) &
        call error_mesg('land_model_init', 'water_capacity must be positive', FATAL)
      if (beta_frac <= 0. .or. beta_frac > 1.) &
        call error_mesg('land_model_init', 'beta_frac must be in (0, 1]', FATAL)
      if (water_init_frac < 0. .or. water_init_frac > 1.) &
        call error_mesg('land_model_init', 'water_init_frac must be in [0, 1]', FATAL)
      Land%do_bucket = .true.
    case default
      call error_mesg('land_model_init', 'unknown hydrology = "'//trim(hydrology)// &
                      '"; use saturated or bucket', FATAL)
    end select
    if (do_snow) then
      if (trim(hydrology) /= 'bucket') &
        call error_mesg('land_model_init', "do_snow needs hydrology = 'bucket'", FATAL)
      if (trim(land_model) /= 'simple') &
        call error_mesg('land_model_init', "do_snow needs land_model = 'simple' (the melt, sublimation "// &
                        'and frost change the surface temperature)', FATAL)
      if (albedo_snow < 0. .or. albedo_snow > 1.) &
        call error_mesg('land_model_init', 'albedo_snow must be in [0, 1]', FATAL)
      if (snow_mask_depth <= 0.) &
        call error_mesg('land_model_init', 'snow_mask_depth must be positive', FATAL)
      if (snow_max <= 0.) &
        call error_mesg('land_model_init', 'snow_max must be positive', FATAL)
      Land%do_snow = .true.
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

  end subroutine land_model_read_nml_and_mask

!#######################################################################

  !> Sets the land heat capacity and the cold-start land temperature (and, with the bucket,
  !> soil water and snow), and registers the diagnostics. A restart
  !> (`INPUT/land_model.res.nc`) is applied afterwards by `flux_exchange_init`.
  subroutine land_model_init_state(Land, Time, grid, t_all, heat_capacity_ocean_default)
    type(land_data_type), intent(inout) :: Land  !! the land state
    type(time_type), intent(in) :: Time  !! current time
    type(surface_grid_type), intent(in) :: grid  !! the surface grid
    real, intent(in) :: t_all(:, :)  !! the ocean's cold-start SST on all points [K]
    real, intent(in) :: heat_capacity_ocean_default  !! the ocean `heat_capacity` [J/m2/K]

    real, dimension(grid%nx, grid%ny) :: t_file
    logical :: used

    allocate (Land%t_surf(grid%nx, grid%ny), Land%cap(grid%nx, grid%ny), &
              Land%dt_t_surf(grid%nx, grid%ny))
    Land%dt_t_surf = 0.0

    if (heat_capacity .le. 0.) heat_capacity = heat_capacity_ocean_default
    if (mpp_pe() == mpp_root_pe()) write (stdlog(), nml=land_model_nml)
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

    id_land_mask = register_static_field(mod_name, 'land_mask', grid%axes, &
                                         'land mask (1: land, 0: ocean)', '1')
    if (id_land_mask > 0) used = send_data(id_land_mask, merge(1., 0., Land%mask))

    id_extra_heating = register_diag_field(mod_name, 'extra_heating', grid%axes, Time, &
                                           'local surface heating over land (debug)', 'K/s')
    if (trim(land_model) == 'simple') &
      id_f_net_final = register_diag_field(mod_name, 'f_net_final', grid%axes, Time, &
                                           'net downward surface heat flux over land balanced by the '// &
                                           'implicit step, at the end of the step (debug)', 'W/m2')

! Bucket and snow: the cold-start state (a restart replaces it in land_model_read_restart)
    if (trim(hydrology) == 'bucket') then
      allocate (Land%water(grid%nx, grid%ny), Land%beta(grid%nx, grid%ny), &
                Land%runoff(grid%nx, grid%ny), Land%water_clip(grid%nx, grid%ny))
      Land%water = 0.0
      where (Land%mask) Land%water = water_init_frac*water_capacity
      Land%beta = 0.0
      Land%runoff = 0.0
      Land%water_clip = 0.0
      id_water = register_diag_field(mod_name, 'water', grid%axes, Time, 'soil water', 'kg/m2')
      id_beta = register_diag_field(mod_name, 'beta', grid%axes, Time, &
                                    'evaporation efficiency (fraction of the potential evaporation)', '1')
      id_runoff = register_diag_field(mod_name, 'runoff', grid%axes, Time, 'runoff', 'kg/m2/s')
      id_water_clip = register_diag_field(mod_name, 'water_clip', grid%axes, Time, &
                                          'water added to keep the soil water non-negative', 'kg/m2/s')
    end if
    if (do_snow) then
      allocate (Land%snow(grid%nx, grid%ny), Land%melt(grid%nx, grid%ny), &
                Land%sublimation(grid%nx, grid%ny), Land%frost(grid%nx, grid%ny))
      Land%snow = 0.0
      Land%melt = 0.0
      Land%sublimation = 0.0
      Land%frost = 0.0
      id_snow = register_diag_field(mod_name, 'snow', grid%axes, Time, 'snow', 'kg/m2')
      id_melt = register_diag_field(mod_name, 'melt', grid%axes, Time, 'snow melt', 'kg/m2/s')
      id_sublimation = register_diag_field(mod_name, 'sublimation', grid%axes, Time, &
                                           'sublimation of snow', 'kg/m2/s')
      id_frost = register_diag_field(mod_name, 'frost', grid%axes, Time, &
                                     'frost (deposition on snow)', 'kg/m2/s')
    end if

  end subroutine land_model_init_state

!#######################################################################

  !> Reads the soil water and snow from `INPUT/land_model.res.nc` (called by
  !> `flux_exchange_init` when the surface restarts exist).
  !>
  !> Only with the bucket (`water`) and snow (`snow`). A field missing from the file (for
  !> example after switching on the bucket or the snow in the middle of an experiment) keeps
  !> its cold-start value, with a NOTE. The surface temperature is read by
  !> `flux_exchange_init`.
  subroutine land_model_read_restart(Land, grid)
    type(land_data_type), intent(inout) :: Land  !! the land state
    type(surface_grid_type), intent(in) :: grid  !! the surface grid

    type(restart_file_type) :: rst
    real, dimension(grid%nx, grid%ny) :: field

    if (.not. Land%do_bucket) return
    if (.not. open_restart_read(rst, 'INPUT/land_model.res.nc', grid%domain)) &
      call error_mesg('land_model_read_restart', 'cannot open INPUT/land_model.res.nc', FATAL)
    if (restart_field_exists(rst, 'water')) then
      call read_restart_field(rst, 'water', field)
      where (Land%mask) Land%water = field
    else
      call error_mesg('land_model_read_restart', 'INPUT/land_model.res.nc has no soil water (water): '// &
                      'cold start of the bucket (water_init_frac)', NOTE)
    end if
    if (Land%do_snow) then
      if (restart_field_exists(rst, 'snow')) then
        call read_restart_field(rst, 'snow', field)
        where (Land%mask) Land%snow = field
      else
        call error_mesg('land_model_read_restart', 'INPUT/land_model.res.nc has no snow: '// &
                        'cold start of the snow (no snow)', NOTE)
      end if
    end if
    call close_restart(rst)

  end subroutine land_model_read_restart

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

  !> Sets the evaporation efficiency `beta` of the bucket and the snow albedo, before the
  !> surface fluxes are computed (called by `sfc_boundary_layer` only with the bucket).
  !>
  !> On land points, `beta = 1` where there is snow, else `min(W/(beta_frac*water_capacity), 1)`.
  !> With `do_snow`, the albedo becomes
  !> `max(alb, alb + (albedo_snow - alb)*S/(S + snow_mask_depth))`, where `alb` is the
  !> background albedo: snow brightens the surface but never darkens a background that is
  !> brighter than `albedo_snow`.
  subroutine land_model_boundary(Land, albedo)
    type(land_data_type), intent(inout) :: Land  !! the land state; on return `Land%beta`
    real, intent(inout) :: albedo(:, :)  !! surface albedo: the background on entry, with snow on return

    if (.not. Land%do_bucket) return
    where (Land%mask) Land%beta = min(Land%water/(beta_frac*water_capacity), 1.0)
    if (Land%do_snow) then
      where (Land%mask .and. Land%snow > 0.0) Land%beta = 1.0
      where (Land%mask) &
        albedo = max(albedo, albedo + (albedo_snow - albedo)*Land%snow/(Land%snow + snow_mask_depth))
    end if

  end subroutine land_model_boundary

!#######################################################################

  !> Steps the land surface temperature (and the bucket and snow) over one atmospheric time
  !> step.
  !>
  !> `'simple'`: implicit surface energy budget; `'prescribed'`: change to the temperature of
  !> `t_surf_file`; `'fixed'`: no change. The local surface heating is a heat flux
  !> `cap*extra_heating` in the implicit step of `'simple'`, and is added as
  !> `extra_heating*dt` to the change in the other modes. Without snow the snowfall melts on arrival, taking its latent heat of fusion in
  !> the implicit step; with snow it is stored and the implicit step has no melt term. With
  !> the bucket, the water and snow are then updated (`bucket_step`). `Land%dt_t_surf` is the
  !> change before the snow adjustments of the surface temperature (the change from which the
  !> fluxes returned to the atmosphere are completed).
  subroutine update_land_model_fast(Exch, Land)
    type(surface_exchange_type), intent(in) :: Exch  !! fluxes and derivatives from the flux exchange
    type(land_data_type), intent(inout) :: Land  !! the land state

    real, dimension(size(Land%t_surf, 1), size(Land%t_surf, 2)) :: t_new, melt_heat, field
    logical :: used

    select case (trim(land_model))
    case ('simple')
! the heat taken by the snowfall that melts on arrival (none with snow: it is stored)
      if (Land%do_snow) then
        melt_heat = 0.0
      else
        melt_heat = hlf*Exch%fprec
      end if
      where (Land%mask)
        Land%dt_t_surf = slab_dt(Exch%flux_lw, Exch%flux_sw, melt_heat, Exch%flux_t, &
                                 Exch%flux_q, Exch%qflux_land + Land%cap*Exch%extra_heating, &
                                 Exch%dhdt_surf, Exch%dedt_surf, Exch%drdt_surf, Exch%dt, Land%cap)
      elsewhere
        Land%dt_t_surf = 0.0
      end where
! Debug: the net downward heat flux at the end of the step, from the fluxes; the implicit
! step balances cap*dt_t/dt = f_net_final + cap*extra_heating
      if (id_f_net_final > 0) then
        field = 0.0
        where (Land%mask) field = Exch%flux_lw + Exch%flux_sw - melt_heat &
                                  - (Exch%flux_t + hlv*Exch%flux_q) + Exch%qflux_land &
                                  - (Exch%dhdt_surf + hlv*Exch%dedt_surf + Exch%drdt_surf)*Land%dt_t_surf
        used = send_data(id_f_net_final, field, Exch%Time)
      end if
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

! 'simple' has the local surface heating in its implicit step; the other modes add it here
    if (trim(land_model) /= 'simple') then
      where (Land%mask) Land%dt_t_surf = Land%dt_t_surf + Exch%extra_heating*Exch%dt
    end if
    where (Land%mask) Land%t_surf = Land%t_surf + Land%dt_t_surf

    if (Land%do_bucket) call bucket_step(Exch, Land)

    if (id_extra_heating > 0) then
      field = 0.0
      where (Land%mask) field = Exch%extra_heating
      used = send_data(id_extra_heating, field, Exch%Time)
    end if
    if (Land%do_bucket) then
      if (id_water > 0) used = send_data(id_water, Land%water, Exch%Time)
      if (id_beta > 0) used = send_data(id_beta, Land%beta, Exch%Time)
      if (id_runoff > 0) used = send_data(id_runoff, Land%runoff, Exch%Time)
      if (id_water_clip > 0) used = send_data(id_water_clip, Land%water_clip, Exch%Time)
    end if
    if (Land%do_snow) then
      if (id_snow > 0) used = send_data(id_snow, Land%snow, Exch%Time)
      if (id_melt > 0) used = send_data(id_melt, Land%melt, Exch%Time)
      if (id_sublimation > 0) used = send_data(id_sublimation, Land%sublimation, Exch%Time)
      if (id_frost > 0) used = send_data(id_frost, Land%frost, Exch%Time)
    end if

  end subroutine update_land_model_fast

!#######################################################################

  ! Updates the soil water W and the snow S on the land points after the surface temperature
  ! step, conserving water exactly (apart from water_clip):
  !
  !   E = flux_q + dt_t*dedt_surf  (the final evaporation, as in flux_up_to_atmos)
  !   with snow:  S += fprec*dt; E >= 0: sublimation Es = min(E*dt, S)/dt from the snow, the
  !               rest Ew = E - Es from the bucket; E < 0 on snow: frost (S -= E*dt); E < 0
  !               without snow: dew into the bucket (Ew = E). Sublimation takes hlf*Es*dt/cap
  !               from T and frost gives hlf*(-E)*dt/cap. Then, above tfreeze, the snow melts
  !               into the bucket, limited by the heat cap*(T - tfreeze) and by the snow, and
  !               snow above snow_max runs off.
  !   without snow: W += fprec*dt (the snowfall melted in the implicit step), Ew = E
  !   W += (lprec - Ew)*dt; W above water_capacity runs off; W < 0 is set to 0 (water_clip).
  !
  ! These snow adjustments of T are explicit and surface-only: dt_t_surf, from which the
  ! atmosphere's fluxes are completed, is not changed.
  subroutine bucket_step(Exch, Land)
    type(surface_exchange_type), intent(in) :: Exch
    type(land_data_type), intent(inout) :: Land

    integer :: i, j
    real :: dt, e, ew, es_mass, m, w, s, t, cap

    dt = Exch%dt
    do j = 1, size(Land%t_surf, 2)
      do i = 1, size(Land%t_surf, 1)
        if (.not. Land%mask(i, j)) cycle
        e = Exch%flux_q(i, j) + Land%dt_t_surf(i, j)*Exch%dedt_surf(i, j)
        w = Land%water(i, j)
        Land%runoff(i, j) = 0.0
        Land%water_clip(i, j) = 0.0
        if (Land%do_snow) then
          cap = Land%cap(i, j)
          t = Land%t_surf(i, j)
          s = Land%snow(i, j) + Exch%fprec(i, j)*dt
          Land%melt(i, j) = 0.0
          Land%sublimation(i, j) = 0.0
          Land%frost(i, j) = 0.0
          if (e >= 0.0) then
            ! sublimation first, the rest from the bucket
            es_mass = min(e*dt, s)
            s = s - es_mass
            Land%sublimation(i, j) = es_mass/dt
            ew = e - Land%sublimation(i, j)
            t = t - hlf*es_mass/cap
          else if (s > 0.0) then
            ! frost on the snow
            s = s - e*dt
            Land%frost(i, j) = -e
            ew = 0.0
            t = t + hlf*(-e)*dt/cap
          else
            ! dew into the bucket
            ew = e
          end if
          if (t > tfreeze .and. s > 0.0) then
            m = min(s, cap*(t - tfreeze)/hlf)
            t = t - hlf*m/cap
            s = s - m
            w = w + m
            Land%melt(i, j) = m/dt
          end if
          if (s > snow_max) then
            Land%runoff(i, j) = (s - snow_max)/dt
            s = snow_max
          end if
          Land%snow(i, j) = s
          Land%t_surf(i, j) = t
        else
          w = w + Exch%fprec(i, j)*dt
          ew = e
        end if
        w = w + (Exch%lprec(i, j) - ew)*dt
        if (w > water_capacity) then
          Land%runoff(i, j) = Land%runoff(i, j) + (w - water_capacity)/dt
          w = water_capacity
        end if
        if (w < 0.0) then
          Land%water_clip(i, j) = -w/dt
          w = 0.0
        end if
        Land%water(i, j) = w
      end do
    end do

  end subroutine bucket_step

!#######################################################################

  !> Writes the restart file `RESTART/land_model.res.nc`: the surface temperature, and with
  !> the bucket the soil water (`water`) and with snow the snow (`snow`); all 0 over the ocean.
  subroutine land_model_end(Land, grid)
    type(land_data_type), intent(in) :: Land  !! the land state
    type(surface_grid_type), intent(in) :: grid  !! the surface grid
    type(restart_file_type) :: rst

    call open_restart_write(rst, 'RESTART/land_model.res.nc', grid%domain)
    call write_restart_field(rst, 't_surf', Land%t_surf)
    if (Land%do_bucket) call write_restart_field(rst, 'water', Land%water)
    if (Land%do_snow) call write_restart_field(rst, 'snow', Land%snow)
    call close_restart(rst)

  end subroutine land_model_end

end module land_model_mod
