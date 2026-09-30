!> Initial fields of the spectral dynamical core on a cold start.
!>
!> Starts from an isothermal atmosphere at rest whose surface pressure is in hydrostatic
!> balance with the topography, then, depending on `choice_of_init`, perturbs the
!> temperature at one point, adds a small vorticity perturbation, or reads the winds,
!> temperature and surface pressure from `INPUT/initial_conditions.nc`. The fields are
!> transformed to spectral space and back, so the grid fields are spectrally truncated.
module spectral_initialize_fields_mod

  use fms_mod, only: mpp_pe, mpp_root_pe, write_version_number, FATAL, error_mesg

  use netcdf, only: nf90_open, nf90_close, nf90_inq_varid, nf90_inquire_variable, &
                    nf90_inquire_dimension, nf90_get_var, nf90_strerror, nf90_noerr, &
                    nf90_nowrite, nf90_max_var_dims

  use constants_mod, only: rdgas

  use transforms_mod, only: trans_grid_to_spherical, trans_spherical_to_grid, vor_div_from_uv_grid, &
                            uv_grid_from_vor_div, get_grid_domain, get_spec_domain, area_weighted_global_mean, &
                            get_lon_max, get_lat_max

  implicit none
  private

  public :: spectral_initialize_fields, read_initial_condition

  !> Reads one field of `INPUT/initial_conditions.nc` on this PE's grid subdomain.
  !>
  !> The variable must have the dimensions (lon, lat) or (lon, lat, level) of the model
  !> grid; any further dimensions (e.g. time) must have length 1.
  interface read_initial_condition
    module procedure read_initial_condition_2d, read_initial_condition_3d
  end interface

  character(len=*), parameter :: ic_file = 'INPUT/initial_conditions.nc'

  character(len=128), parameter :: version = &
                                   '$Id: spectral_initialize_fields.f90,v 10.0 2003/10/24 22:00:59 fms Exp $'

  character(len=128), parameter :: tagname = &
                                   '$Name: lima $'

  logical :: entry_to_logfile_done = .false.

contains

!-------------------------------------------------------------------------------------------------
  !> Computes the initial grid and spectral fields.
  subroutine spectral_initialize_fields(reference_sea_level_press, triang_trunc, choice_of_init, initial_temperature, &
                                        surf_geopotential, ln_ps, vors, divs, ts, psg, ug, vg, tg, vorg, divg)

    real, intent(in) :: reference_sea_level_press  !! initial surface pressure where the surface height is 0 [Pa]
    logical, intent(in) :: triang_trunc  !! triangular (`.true.`) or rhomboidal truncation
    integer, intent(in) :: choice_of_init
    !! 1: add 1 K to the temperature of the first grid column; 2: small vorticity perturbation
    !! in the lowest three levels; 3: read `ucomp`, `vcomp`, `temp` and `ps` from
    !! `INPUT/initial_conditions.nc`
    real, intent(in) :: initial_temperature  !! temperature of the isothermal atmosphere [K]

    real, intent(in), dimension(:, :) :: surf_geopotential  !! surface geopotential [m2/s2]
    complex, intent(out), dimension(:, :) :: ln_ps  !! spectral log of surface pressure
    complex, intent(out), dimension(:, :, :) :: vors, divs, ts
    !! spectral vorticity [1/s], divergence [1/s] and temperature [K]
    real, intent(out), dimension(:, :) :: psg  !! grid surface pressure [Pa]
    real, intent(out), dimension(:, :, :) :: ug, vg, tg
    !! grid zonal and meridional wind [m/s] and temperature [K]
    real, intent(out), dimension(:, :, :) :: vorg, divg  !! grid vorticity and divergence [1/s]

    real, allocatable, dimension(:, :) :: ln_psg

    real :: initial_sea_level_press, global_mean_psg
    real :: initial_perturbation = 1.e-7

    integer :: ms, me, ns, ne, is, ie, js, je, num_levels

    if (.not. entry_to_logfile_done) then
      call write_version_number(version, tagname)
      entry_to_logfile_done = .true.
    end if

    num_levels = size(ug, 3)
    call get_grid_domain(is, ie, js, je)
    call get_spec_domain(ms, me, ns, ne)
    allocate (ln_psg(is:ie, js:je))

    initial_sea_level_press = reference_sea_level_press

    ug = 0.
    vg = 0.
    tg = 0.
    psg = 0.
    vorg = 0.
    divg = 0.

    vors = (0., 0.)
    divs = (0., 0.)
    ts = (0., 0.)
    ln_ps = (0., 0.)

    tg = initial_temperature
    ln_psg = log(initial_sea_level_press) - surf_geopotential/(rdgas*initial_temperature)
    psg = exp(ln_psg)

    if (choice_of_init == 1) then  ! perturb temperature field
      if (is <= 1 .and. ie >= 1 .and. js <= 1 .and. je >= 1) then
        tg(1, 1, :) = tg(1, 1, :) + 1.0
      end if
    end if

    if (choice_of_init == 2) then   ! initial vorticity perturbation used in benchmark code
      if (ms <= 1 .and. me >= 1 .and. ns <= 3 .and. ne >= 3) then
        vors(2 - ms, 4 - ns, num_levels) = initial_perturbation
        vors(2 - ms, 4 - ns, num_levels - 1) = initial_perturbation
        vors(2 - ms, 4 - ns, num_levels - 2) = initial_perturbation
      end if
      if (ms <= 5 .and. me >= 5 .and. ns <= 3 .and. ne >= 3) then
        vors(6 - ms, 4 - ns, num_levels) = initial_perturbation
        vors(6 - ms, 4 - ns, num_levels - 1) = initial_perturbation
        vors(6 - ms, 4 - ns, num_levels - 2) = initial_perturbation
      end if
      if (ms <= 1 .and. me >= 1 .and. ns <= 2 .and. ne >= 2) then
        vors(2 - ms, 3 - ns, num_levels) = initial_perturbation
        vors(2 - ms, 3 - ns, num_levels - 1) = initial_perturbation
        vors(2 - ms, 3 - ns, num_levels - 2) = initial_perturbation
      end if
      if (ms <= 5 .and. me >= 5 .and. ns <= 2 .and. ne >= 2) then
        vors(6 - ms, 3 - ns, num_levels) = initial_perturbation
        vors(6 - ms, 3 - ns, num_levels - 1) = initial_perturbation
        vors(6 - ms, 3 - ns, num_levels - 2) = initial_perturbation
      end if
      call uv_grid_from_vor_div(vors, divs, ug, vg)
    end if

! Initial state read from INPUT/initial_conditions.nc (after Lorenzo Polvani)
    if (choice_of_init == 3) then
      call read_initial_condition('ucomp', ug)
      call read_initial_condition('vcomp', vg)
      call read_initial_condition('temp', tg)
      call read_initial_condition('ps', psg)
      ln_psg = log(psg(:, :))
      if (mpp_pe() == mpp_root_pe()) then
        print *, 'initial dynamical fields read in from initial_conditions.nc'
      end if
    end if

!  initial spectral fields (and spectrally-filtered) grid fields

    call trans_grid_to_spherical(tg, ts)
    call trans_spherical_to_grid(ts, tg)

    call trans_grid_to_spherical(ln_psg, ln_ps)
    call trans_spherical_to_grid(ln_ps, ln_psg)
    psg = exp(ln_psg)

    call vor_div_from_uv_grid(ug, vg, vors, divs, triang=triang_trunc)
    call uv_grid_from_vor_div(vors, divs, ug, vg)
    call trans_spherical_to_grid(vors, vorg)
    call trans_spherical_to_grid(divs, divg)

!  compute and print mean surface pressure
    global_mean_psg = area_weighted_global_mean(psg)
    if (mpp_pe() == mpp_root_pe()) then
      print '("mean surface pressure=",f9.4," mb")', .01*global_mean_psg
    end if

    return
  end subroutine spectral_initialize_fields
!================================================================================

  subroutine read_initial_condition_3d(name, field)
    character(len=*), intent(in) :: name  !! variable name in the file
    real, intent(out), dimension(:, :, :) :: field  !! the field on this PE's subdomain
    integer :: ncid, varid, is, ie, js, je

    call open_initial_condition(name, 3, size(field, 3), ncid, varid)
    call get_grid_domain(is, ie, js, je)
    call check_nc(nf90_get_var(ncid, varid, field, start=(/is, js, 1/)), 'Could not read '//name//' from')
    call check_nc(nf90_close(ncid), 'Could not close')

  end subroutine read_initial_condition_3d
!================================================================================

  subroutine read_initial_condition_2d(name, field)
    character(len=*), intent(in) :: name  !! variable name in the file
    real, intent(out), dimension(:, :) :: field  !! the field on this PE's subdomain
    integer :: ncid, varid, is, ie, js, je

    call open_initial_condition(name, 2, 1, ncid, varid)
    call get_grid_domain(is, ie, js, je)
    call check_nc(nf90_get_var(ncid, varid, field, start=(/is, js/)), 'Could not read '//name//' from')
    call check_nc(nf90_close(ncid), 'Could not close')

  end subroutine read_initial_condition_2d
!================================================================================

  !> Opens the file and finds variable `name`, checking that its dimensions are
  !> (lon, lat[, level]) of the model grid. Any further (e.g. time) dimensions
  !> must have length 1.
  subroutine open_initial_condition(name, rank, num_levels, ncid, varid)
    character(len=*), intent(in) :: name
    integer, intent(in) :: rank, num_levels
    integer, intent(out) :: ncid, varid
    integer :: ndims, n, lon_max, lat_max
    integer, dimension(nf90_max_var_dims) :: dimids, file_len, model_len
    character(len=64) :: file_shape, model_shape, axes

    call check_nc(nf90_open(ic_file, nf90_nowrite, ncid), 'Could not open')
    call check_nc(nf90_inq_varid(ncid, name, varid), 'Could not find variable '//name//' in')
    call check_nc(nf90_inquire_variable(ncid, varid, ndims=ndims, dimids=dimids), 'Could not inquire '//name//' in')

    call get_lon_max(lon_max)
    call get_lat_max(lat_max)
    model_len = 1
    model_len(1:3) = (/lon_max, lat_max, num_levels/)
    do n = 1, ndims
      call check_nc(nf90_inquire_dimension(ncid, dimids(n), len=file_len(n)), 'Could not inquire dimensions of '//name//' in')
    end do
    if (ndims < rank .or. any(file_len(1:ndims) /= model_len(1:ndims))) then
      write (file_shape, '(*(i0,:," x "))') file_len(1:ndims)
      write (model_shape, '(*(i0,:," x "))') model_len(1:rank)
      axes = 'lon x lat'
      if (rank == 3) axes = 'lon x lat x level'
      call error_mesg('spectral_initialize_fields', 'Variable '//name//' in '//ic_file//' has shape '// &
                      trim(file_shape)//', but the model grid needs '//trim(model_shape)//' ('//trim(axes)//')', FATAL)
    end if

  end subroutine open_initial_condition
!================================================================================

  subroutine check_nc(status, action)
    integer, intent(in) :: status
    character(len=*), intent(in) :: action

    if (status /= nf90_noerr) then
      call error_mesg('spectral_initialize_fields', action//' '//ic_file//': '//trim(nf90_strerror(status)), FATAL)
    end if

  end subroutine check_nc
!================================================================================

end module spectral_initialize_fields_mod
