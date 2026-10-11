!> The ocean surface: a slab (mixed-layer) ocean, or a prescribed or fixed sea surface
!> temperature.
!>
!> The ocean covers the points that are not land (`land_model_mod`). Its surface temperature
!> `t_surf` is
!>
!> * `'slab'`: stepped with the implicit surface energy budget (`slab_tendency_mod`), with a
!>   heat capacity that can vary with latitude and the prescribed ocean heat fluxes of
!>   `ocean_qflux_mod` (by default applied over the ocean only, with their mean removed);
!> * `'prescribed'`: read from `sst_file` at every step;
!> * `'fixed'`: held at its initial value.
!>
!> In every mode the local surface heating of `local_heating_nml` is added. On a cold start
!> the SST is set by `sst_init`: the profile `Tm - deltaT*(3 sin^2(lat) - 1)/3`, a uniform
!> `Tm`, one of the Aqua-Planet Experiment profiles, or `sst_file`. The cold-start field and
!> the heat-capacity profile are computed on all points, so that the land model can start
!> from the same field. The state is saved in `RESTART/ocean_model.res.nc`, which takes
!> precedence over `sst_init` when it exists in `INPUT/`.
!>
!> Namelist: `ocean_model_nml`
!> ([namelist reference](https://eddy-stanford.github.io/MiMA/Parameters/#ocean_model_nml)).
module ocean_model_mod

  use fms_mod, only: input_nml_file, check_nml_error, error_mesg, FATAL, NOTE, &
                     mpp_pe, mpp_root_pe, stdlog
  use restart_file_mod, only: restart_file_type, open_restart_write, close_restart, &
                              write_restart_field
  use diag_manager_mod, only: register_diag_field, send_data
  use time_manager_mod, only: time_type
  use constants_mod, only: hlf
  use mima_interpolator_mod, only: interpolate_type, interpolator_init, CONSTANT, interpolator
  use surface_types_mod, only: surface_grid_type, surface_exchange_type
  use slab_tendency_mod, only: slab_dt
  use ocean_qflux_mod, only: ocean_qflux_init

  implicit none
  private

  public :: ocean_data_type, ocean_model_init, update_ocean_model_fast, &
            update_ocean_model_slow, ocean_model_end

  !> State of the ocean surface. All arrays are on this PE's grid; off the ocean mask they
  !> are 0.
  type ocean_data_type
    logical, allocatable :: mask(:, :)  !! ocean points
    real, allocatable :: t_surf(:, :)  !! surface temperature [K]
    real, allocatable :: cap(:, :)  !! mixed-layer heat capacity [J/m2/K]
    real, allocatable :: qflux(:, :)
    !! prescribed ocean heat flux convergence (Q-flux) applied at each point: by the slab over
    !! the ocean and, with `qflux_ocean_only = .false.` (v1), by the land model over land;
    !! otherwise 0 over land [W/m2]
    real, allocatable :: dt_t_surf(:, :)  !! change of `t_surf` over the last step [K]
    real :: rough_ocean  !! roughness length of the ocean for momentum, heat and moisture [m]
  end type ocean_data_type

  character(len=*), parameter :: mod_name = 'ocean'

  character(len=32) :: ocean_model = 'slab'
  !! `'slab'`: mixed layer with an interactive SST; `'prescribed'`: SST read from `sst_file` at
  !! every step (implies `sst_init = 'file'`); `'fixed'`: SST held at its initial value
  character(len=32) :: sst_init = 'profile'
  !! initial SST on a cold start: `'profile'`: `Tm - deltaT*(3 sin^2(lat) - 1)/3`; `'uniform'`:
  !! `Tm` everywhere; `'ape_control_plus2'`, `'ape_control'`, `'ape_qobs'`, `'ape_qobs5'`: the
  !! Aqua-Planet Experiment profiles "control+2", "control", "qobs" and "qobs+5"; `'file'`:
  !! from `sst_file`
  real :: Tm = 285.  !! [K] mean initial SST of `'profile'`; the SST of `'uniform'`
  real :: deltaT = 40.  !! [K] equator-to-pole difference of the initial SST of `'profile'`
  character(len=256) :: sst_file = ''  !! SST file name, without `.nc`, in `INPUT/`
  character(len=256) :: sst_file_var = ''  !! name of the SST variable in `sst_file` (empty: `sst_file`)
  real :: rough_ocean = 3.21e-05  !! [m] roughness length of the ocean (momentum, heat and moisture)
  real :: heat_capacity = 3.e08  !! [J/m2/K] mixed-layer heat capacity (poleward of `heat_cap_limit`)
  real :: trop_capacity = 1.e08
  !! [J/m2/K] heat capacity equatorward of `trop_cap_limit`, varying linearly to
  !! `heat_capacity` at `heat_cap_limit`; <= 0: `heat_capacity`
  real :: trop_cap_limit = 20.  !! [deg] latitude of the tropical heat capacity `trop_capacity`
  real :: heat_cap_limit = 60.  !! [deg] latitude of the extratropical heat capacity `heat_capacity`
  real :: np_cap_factor = 1.  !! factor on `heat_capacity` in the Northern Hemisphere

  namelist /ocean_model_nml/ ocean_model, sst_init, Tm, deltaT, sst_file, sst_file_var, &
    rough_ocean, heat_capacity, trop_capacity, trop_cap_limit, heat_cap_limit, np_cap_factor

  type(interpolate_type), save :: sst_interp
  character(len=256) :: sst_var
  integer :: id_qflux, id_heat

contains

!#######################################################################

  !> Reads `ocean_model_nml`, sets up the heat capacity, the Q-flux and the cold-start SST,
  !> and registers the diagnostics.
  !>
  !> The cold-start SST and the heat-capacity profile are also returned on all points, for
  !> the land model. A restart (`INPUT/ocean_model.res.nc`) is applied afterwards by
  !> `flux_exchange_init`.
  subroutine ocean_model_init(Ocean, Time, grid, mask, t_all, cap_all, cap_default)
    type(ocean_data_type), intent(inout) :: Ocean  !! the ocean state
    type(time_type), intent(in) :: Time  !! current time
    type(surface_grid_type), intent(in) :: grid  !! the surface grid
    logical, intent(in) :: mask(:, :)  !! ocean points
    real, intent(out) :: t_all(:, :)  !! cold-start SST on all points [K]
    real, intent(out) :: cap_all(:, :)  !! mixed-layer heat-capacity profile on all points [J/m2/K]
    real, intent(out) :: cap_default  !! `heat_capacity` [J/m2/K]

    integer :: ierr, io, j
    real :: lat, pi, xx, loc_cap

    pi = 4.0*atan(1.)

    read (input_nml_file, nml=ocean_model_nml, iostat=io)
    ierr = check_nml_error(io, 'ocean_model_nml')

    select case (trim(ocean_model))
    case ('slab', 'fixed')
    case ('prescribed')
      if (trim(sst_init) /= 'file') then
        call error_mesg('ocean_model_init', "ocean_model = 'prescribed': the initial SST is "// &
                        "taken from sst_file (sst_init = 'file')", NOTE)
        sst_init = 'file'
      end if
    case default
      call error_mesg('ocean_model_init', 'unknown ocean_model = "'//trim(ocean_model)// &
                      '"; use slab, prescribed or fixed', FATAL)
    end select
    select case (trim(sst_init))
    case ('profile', 'uniform', 'ape_control_plus2', 'ape_control', 'ape_qobs', 'ape_qobs5')
    case ('file')
      if (len_trim(sst_file) == 0) call error_mesg('ocean_model_init', &
                                                   "ocean_model = 'prescribed' and sst_init = 'file' need sst_file", FATAL)
    case default
      call error_mesg('ocean_model_init', 'unknown sst_init = "'//trim(sst_init)//'"; use profile, '// &
                      'uniform, ape_control_plus2, ape_control, ape_qobs, ape_qobs5 or file', FATAL)
    end select
    if (trim(sst_init) == 'profile' .and. (Tm .le. 0. .or. Tm .ge. 399.)) &
      call error_mesg('ocean_model_init', "sst_init = 'profile' with Tm <= 0 or Tm >= 399: these "// &
                      "v1 magic values are no longer used. For a uniform SST set sst_init = 'uniform' "// &
                      "and Tm to the (positive) SST; for Tm >= 399 use sst_init = 'ape_control_plus2', "// &
                      "'ape_control', 'ape_qobs' or 'ape_qobs5' (see the table in Migration_v2.md)", FATAL)
    if (trop_capacity .le. 0.) trop_capacity = heat_capacity
    if (mpp_pe() == mpp_root_pe()) write (stdlog(), nml=ocean_model_nml)

    sst_var = sst_file_var
    if (len_trim(sst_var) == 0) sst_var = sst_file

    allocate (Ocean%mask(grid%nx, grid%ny), Ocean%t_surf(grid%nx, grid%ny), &
              Ocean%cap(grid%nx, grid%ny), Ocean%qflux(grid%nx, grid%ny), &
              Ocean%dt_t_surf(grid%nx, grid%ny))
    Ocean%mask = mask
    Ocean%rough_ocean = rough_ocean
    Ocean%dt_t_surf = 0.0

!mj read fixed SSTs
    if (trim(sst_init) == 'file') then
      call interpolator_init(sst_interp, trim(sst_file)//'.nc', grid%lon_bnd, grid%lat_bnd, &
                             data_names=(/sst_var/), data_out_of_bounds=(/CONSTANT/))
    end if

! Heat capacity profile (all points)
    cap_default = heat_capacity
    cap_all = heat_capacity
    !mj ocean depth function of latitude
    if (trop_capacity .ne. heat_capacity .or. np_cap_factor .ne. 1.0) then
      do j = 1, size(cap_all, 2)
        lat = 0.5*180/pi*(grid%lat_bnd(j + 1) + grid%lat_bnd(j))
        if (lat > 0.) then
          loc_cap = heat_capacity*np_cap_factor
        else
          loc_cap = heat_capacity
        end if
        if (abs(lat) < trop_cap_limit) then
          cap_all(:, j) = trop_capacity
        elseif (abs(lat) < heat_cap_limit) then
          cap_all(:, j) = trop_capacity*(1.-(abs(lat) - trop_cap_limit)/(heat_cap_limit - trop_cap_limit)) &
                          + (abs(lat) - trop_cap_limit)/(heat_cap_limit - trop_cap_limit)*loc_cap
        elseif (lat > heat_cap_limit) then
          cap_all(:, j) = loc_cap
        end if
      end do
    end if

! Cold-start SST (all points)
    if (trim(sst_init) == 'file') then
      call interpolator(sst_interp, Time, t_all, trim(sst_var))
    else
      do j = 1, size(t_all, 2)
        lat = 0.5*(grid%lat_bnd(j + 1) + grid%lat_bnd(j))
        xx = sin(lat)*sin(lat)
        select case (trim(sst_init))
        case ('profile')
! CL functional form
          t_all(:, j) = Tm - deltaT*(3.*xx - 1.)/3.
        case ('uniform')
! mj equal temperature
          t_all(:, j) = Tm
        case ('ape_control_plus2')
! APE "control+2" experiment
          t_all(:, j) = 273.15 + 29.0*(1.0 - sin(1.40369*lat)*sin(1.40369*lat))
          if (abs(lat) .ge. pi/3.) t_all(:, j) = 273.15
        case ('ape_control')
! APE "control" experiment
          t_all(:, j) = 273.15 + 27.0*(1.0 - sin(1.5*lat)*sin(1.5*lat))
          if (abs(lat) .ge. pi/3.) t_all(:, j) = 273.15
        case ('ape_qobs', 'ape_qobs5')
! APE "qobs" experiment
          t_all(:, j) = 273.15 + .5*27.*(1.-sin(1.5*lat)*sin(1.5*lat)) &
                        + .5*27.*(1.-sin(1.5*lat)*sin(1.5*lat)*sin(1.5*lat)*sin(1.5*lat))
          if (abs(lat) .ge. pi/3.) t_all(:, j) = 273.15
! APE "qobs+5" experiment
          if (trim(sst_init) == 'ape_qobs5') t_all(:, j) = t_all(:, j) + 5.
        end select
      end do
    end if

    Ocean%cap = 0.0
    where (mask) Ocean%cap = cap_all
    Ocean%t_surf = 0.0
    where (mask) Ocean%t_surf = t_all

    call ocean_qflux_init(grid, mask, Ocean%qflux)

    id_qflux = register_diag_field(mod_name, 'qflux', grid%axes, Time, &
                                   'prescribed ocean heat divergence', 'W/m2')
    id_heat = register_diag_field(mod_name, 'heat_capacity_ocean', grid%axes, Time, &
                                  'ocean mixed layer heat capacity (0 over land)', 'J/m2/K')

  end subroutine ocean_model_init

!#######################################################################

  !> Steps the ocean surface temperature over one atmospheric time step.
  !>
  !> `'slab'`: implicit surface energy budget; `'prescribed'`: change to the SST of
  !> `sst_file`; `'fixed'`: no change. The local surface heating is then added in every
  !> mode. Sends the ocean diagnostics.
  subroutine update_ocean_model_fast(Exch, Ocean)
    type(surface_exchange_type), intent(in) :: Exch  !! fluxes and derivatives from the flux exchange
    type(ocean_data_type), intent(inout) :: Ocean  !! the ocean state

    real, dimension(size(Ocean%t_surf, 1), size(Ocean%t_surf, 2)) :: sst_new
    logical :: used

    select case (trim(ocean_model))
    case ('slab')
      where (Ocean%mask)
        Ocean%dt_t_surf = slab_dt(Exch%flux_lw, Exch%flux_sw, hlf*Exch%fprec, Exch%flux_t, &
                                  Exch%flux_q, Ocean%qflux, Exch%dhdt_surf, Exch%dedt_surf, &
                                  Exch%drdt_surf, Exch%dt, Ocean%cap)
      elsewhere
        Ocean%dt_t_surf = 0.0
      end where
    case ('prescribed') !mj sst read from input file
      call interpolator(sst_interp, Exch%Time, sst_new, trim(sst_var))
      where (Ocean%mask)
        Ocean%dt_t_surf = sst_new - Ocean%t_surf
      elsewhere
        Ocean%dt_t_surf = 0.0
      end where
    case ('fixed')
      Ocean%dt_t_surf = 0.0
    end select

    where (Ocean%mask)
      Ocean%dt_t_surf = Ocean%dt_t_surf + Exch%extra_heating*Exch%dt
      Ocean%t_surf = Ocean%t_surf + Ocean%dt_t_surf
    end where

    if (id_qflux > 0) used = send_data(id_qflux, Ocean%qflux, Exch%Time)
    if (id_heat > 0) used = send_data(id_heat, Ocean%cap, Exch%Time)

  end subroutine update_ocean_model_fast

!#######################################################################

  !> Hook for slow ocean processes such as advection (none at present).
  subroutine update_ocean_model_slow(Time, Ocean)
    type(time_type), intent(in) :: Time  !! current time
    type(ocean_data_type), intent(inout) :: Ocean  !! the ocean state
  end subroutine update_ocean_model_slow

!#######################################################################

  !> Writes the restart file `RESTART/ocean_model.res.nc` (the surface temperature, 0 over
  !> land).
  subroutine ocean_model_end(Ocean, grid)
    type(ocean_data_type), intent(in) :: Ocean  !! the ocean state
    type(surface_grid_type), intent(in) :: grid  !! the surface grid
    type(restart_file_type) :: rst

    call open_restart_write(rst, 'RESTART/ocean_model.res.nc', grid%domain)
    call write_restart_field(rst, 't_surf', Ocean%t_surf)
    call close_restart(rst)

  end subroutine ocean_model_end

end module ocean_model_mod
