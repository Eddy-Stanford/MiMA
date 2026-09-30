module spectral_init_cond_mod

  use fms_mod, only: mpp_pe, mpp_root_pe, error_mesg, FATAL, stdlog, &
                     write_version_number, check_nml_error, input_nml_file
  use fms2_io_mod, only: file_exists, FmsNetcdfFile_t, open_file, close_file, read_data, &
                         get_variable_num_dimensions, get_variable_size

  use mpp_domains_mod, only: mpp_get_global_domain

  use constants_mod, only: grav, pi

  use vert_coordinate_mod, only: compute_vert_coord

  use transforms_mod, only: get_grid_boundaries, get_deg_lon, get_deg_lat, trans_grid_to_spherical, &
                            trans_spherical_to_grid, get_grid_domain, get_spec_domain
  use spec_mpp_mod, only: grid_domain, spectral_domain

  use press_and_geopot_mod, only: press_and_geopot_init, pressure_variables

  use spectral_initialize_fields_mod, only: spectral_initialize_fields

  use topog_regularization_mod, only: compute_lambda, regularize

  use topography_mod, only: gaussian_topog_init, get_topog_mean, get_ocean_mask

  implicit none
  private

  character(len=128), parameter :: version = &
                                   '$Id: spectral_init_cond.f90,v 10.0 2003/10/24 22:00:59 fms Exp $'

  character(len=128), parameter :: tagname = &
                                   '$Name: lima $'

  public :: spectral_init_cond

  real :: initial_temperature = 264.

  namelist /spectral_init_cond_nml/ initial_temperature

contains

!=========================================================================================================================

  subroutine spectral_init_cond(reference_sea_level_press, triang_trunc, use_virtual_temperature, topography_option, &
                                vert_coord_option, vert_difference_option, scale_heights, surf_res, &
                                p_press, p_sigma, exponent, ocean_topog_smoothing, pk, bk, vors, divs, &
                                ts, ln_ps, ug, vg, tg, psg, vorg, divg, surf_geopotential, ocean_mask, specify_initial_conditions)

    real, intent(in) :: reference_sea_level_press
    logical, intent(in) :: triang_trunc, use_virtual_temperature
    character(len=*), intent(in) :: topography_option, vert_coord_option, vert_difference_option
    real, intent(in) :: scale_heights, surf_res, p_press, p_sigma, exponent, ocean_topog_smoothing
    real, intent(out), dimension(:) :: pk, bk
    complex, intent(out), dimension(:, :, :) :: vors, divs, ts
    complex, intent(out), dimension(:, :) :: ln_ps
    real, intent(out), dimension(:, :, :) :: ug, vg, tg
    real, intent(out), dimension(:, :) :: psg
    real, intent(out), dimension(:, :, :) :: vorg, divg
    real, intent(out), dimension(:, :) :: surf_geopotential
    logical, optional, intent(in), dimension(:, :) :: ocean_mask
    logical, intent(in) :: specify_initial_conditions   !epg+ray

! epg+ray: choice_of_init is used by spectral_initialize_fields to actually set up initial conditions
    integer :: choice_of_init = 2
    integer :: unit, ierr, io

!------------------------------------------------------------------------------------------------

! epg+ray: if we want to specify the initial conditions, set choice_of_init to 3
    if (specify_initial_conditions) then
      choice_of_init = 3
    end if

    read (input_nml_file, nml=spectral_init_cond_nml, iostat=io)
    ierr = check_nml_error(io, 'spectral_init_cond_nml')
    call write_version_number(version, tagname)
    if (mpp_pe() == mpp_root_pe()) write (stdlog(), nml=spectral_init_cond_nml)

    call compute_vert_coord(vert_coord_option, scale_heights, surf_res, exponent, p_press, p_sigma, reference_sea_level_press, &
                            pk, bk)

    call get_topography(topography_option, ocean_topog_smoothing, surf_geopotential, ocean_mask)
    call press_and_geopot_init(pk, bk, use_virtual_temperature, vert_difference_option, surf_geopotential)

    call spectral_initialize_fields(reference_sea_level_press, triang_trunc, choice_of_init, initial_temperature, &
                                    surf_geopotential, ln_ps, vors, divs, ts, psg, ug, vg, tg, vorg, divg)

    call check_vert_coord(size(ug, 3), psg)

    return
  end subroutine spectral_init_cond

!================================================================================

  subroutine check_vert_coord(num_levels, psg)
    integer, intent(in) :: num_levels
    real, intent(in), dimension(:, :) :: psg
    real, dimension(size(psg, 1), size(psg, 2), num_levels) :: p_full, ln_p_full
    real, dimension(size(psg, 1), size(psg, 2), num_levels + 1) :: p_half, ln_p_half
    integer :: i, j, k

    call pressure_variables(p_half, ln_p_half, p_full, ln_p_full, psg)
    do k = 1, size(p_full, 3)
      do j = 1, size(p_full, 2)
        do i = 1, size(p_full, 1)
          if (p_half(i, j, k + 1) < p_half(i, j, k)) then
            call error_mesg('check_vert_coord', 'Pressure levels intersect.', FATAL)
          end if
        end do
      end do
    end do

    return
  end subroutine check_vert_coord
!================================================================================

  subroutine get_topography(topography_option, ocean_topog_smoothing, surf_geopotential, ocean_mask_in)

    character(len=*), intent(in) :: topography_option
    real, intent(in) :: ocean_topog_smoothing
    real, intent(out), dimension(:, :) :: surf_geopotential
    logical, intent(in), optional, dimension(:, :) :: ocean_mask_in
    real, dimension(size(surf_geopotential, 1)) :: deg_lon
    real, dimension(size(surf_geopotential, 2)) :: deg_lat
    real, dimension(size(surf_geopotential, 1), size(surf_geopotential, 2)) :: surf_height
    logical, dimension(size(surf_geopotential, 1), size(surf_geopotential, 2)) :: ocean_mask
    complex, allocatable, dimension(:, :) :: spec_tmp
    real :: fraction_smoothed, lambda
    integer :: is, ie, js, je, ms, me, ns, ne, global_num_lon, global_num_lat
    real, allocatable, dimension(:) :: blon, blat
    logical :: topo_file_exists, water_file_exists
    integer, dimension(2) :: siz
    type(FmsNetcdfFile_t) :: topog_file
    real, allocatable, dimension(:, :) :: global_height
    character(len=12) :: ctmp1 = '     by     ', ctmp2 = '     by     '

    if (trim(topography_option) == 'flat') then
      surf_geopotential = 0.

    else if (trim(topography_option) == 'input') then
      if (file_exists('INPUT/topography.data.nc')) then
        call mpp_get_global_domain(grid_domain, xsize=global_num_lon, ysize=global_num_lat)
        if (.not. open_file(topog_file, 'INPUT/topography.data.nc', 'read')) &
          call error_mesg('get_topography', 'cannot open INPUT/topography.data.nc', FATAL)
        if (get_variable_num_dimensions(topog_file, 'zsurf') /= 2) &
          call error_mesg('get_topography', 'zsurf in INPUT/topography.data.nc must be 2-D (lon, lat)', FATAL)
        call get_variable_size(topog_file, 'zsurf', siz)
        if (siz(1) == global_num_lon .and. siz(2) == global_num_lat) then
!      Read the global field and keep this PE's part of it
          allocate (global_height(global_num_lon, global_num_lat))
          call read_data(topog_file, 'zsurf', global_height)
          call get_grid_domain(is, ie, js, je)
          surf_height = global_height(is:ie, js:je)
          deallocate (global_height)
          call close_file(topog_file)
        else
          write (ctmp1(1:4), '(i4)') siz(1)
          write (ctmp1(9:12), '(i4)') siz(2)
          write (ctmp2(1:4), '(i4)') global_num_lon
          write (ctmp2(9:12), '(i4)') global_num_lat
          call error_mesg('get_topography', 'Topography file contains data on a '// &
                          ctmp1//' grid, but atmos model grid is '//ctmp2, FATAL)
        end if

!    Spectrally truncate the topography
        call get_spec_domain(ms, me, ns, ne)
        allocate (spec_tmp(ms:me, ns:ne))
        call trans_grid_to_spherical(surf_height, spec_tmp)
        call trans_spherical_to_grid(spec_tmp, surf_height)
        deallocate (spec_tmp)
        surf_geopotential = grav*surf_height
      else
        call error_mesg('get_topography', 'topography_option="'//trim(topography_option)//'"'// &
                        ' but INPUT/topography.data.nc does not exist', FATAL)
      end if

    else if (trim(topography_option) == 'interpolated') then

!  Get realistic topography
      call get_grid_domain(is, ie, js, je)
      allocate (blon(is:ie + 1), blat(js:je + 1))
      call get_grid_boundaries(blon, blat)
      topo_file_exists = get_topog_mean(blon, blat, surf_height)
      if (.not. topo_file_exists) then
        call error_mesg('get_topography', 'topography_option="'//trim(topography_option)//'"'// &
                        ' but topography data file does not exist', FATAL)
      end if
      surf_geopotential = grav*surf_height

      if (ocean_topog_smoothing == 0.) then
!    Spectrally truncate the realistic topography
        call get_spec_domain(ms, me, ns, ne)
        allocate (spec_tmp(ms:me, ns:ne))
        call trans_grid_to_spherical(surf_geopotential, spec_tmp)
        call trans_spherical_to_grid(spec_tmp, surf_geopotential)
        deallocate (spec_tmp)
      else
!    Do topography regularization
        if (present(ocean_mask_in)) then
          ocean_mask = ocean_mask_in
        else
          water_file_exists = get_ocean_mask(blon, blat, ocean_mask)
          if (.not. water_file_exists) then
            call error_mesg('get_topography', 'topography_option="'//trim(topography_option)//'"'// &
                            ' and ocean_mask is not present but water data file does not exist', FATAL)
          end if
        end if
        call compute_lambda(ocean_topog_smoothing, ocean_mask, surf_geopotential, lambda, fraction_smoothed)

!  Note that the array surf_height is used here for the smoothed surf_geopotential,
!  then immediately loaded back into surf_geopotential
        call regularize(lambda, ocean_mask, surf_geopotential, surf_height, fraction_smoothed)
        surf_geopotential = surf_height

        if (mpp_pe() == mpp_root_pe()) then
          print '(/,"Message from subroutine get_topography:")'
          print '("lambda=",1pe16.8,"  fraction_smoothed=",1pe16.8,/)', lambda, fraction_smoothed
        end if
      end if
      deallocate (blon, blat)

    else if (trim(topography_option) == 'gaussian') then
      call get_deg_lon(deg_lon)
      call get_deg_lat(deg_lat)
      call gaussian_topog_init(deg_lon*pi/180, deg_lat*pi/180, surf_height)
      surf_geopotential = grav*surf_height
    else
      call error_mesg('get_topography', '"'//trim(topography_option)//'" is an invalid value for topography_option.', FATAL)
    end if

    return
  end subroutine get_topography
!================================================================================
end module spectral_init_cond_mod
