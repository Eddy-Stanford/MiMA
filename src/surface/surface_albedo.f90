!> The background surface albedo, a function of latitude (and, for `'tanh_desert'`, of
!> longitude).
!>
!> The albedo is set once at initialization from one of the profiles of
!> `albedo_profile`, which increase from `const_albedo` at low latitudes to
!> `higher_albedo` at high latitudes. It is the same over ocean and land.
!>
!> Namelist: `surface_albedo_nml`
!> ([namelist reference](https://eddy-stanford.github.io/MiMA/Parameters/#surface_albedo_nml)).
module surface_albedo_mod

  use fms_mod, only: input_nml_file, check_nml_error, error_mesg, FATAL, &
                     mpp_pe, mpp_root_pe, stdlog
  use surface_types_mod, only: surface_grid_type

  implicit none
  private

  public :: surface_albedo_init

  character(len=32) :: albedo_profile = 'tanh_desert'
  !! latitude profile of the albedo: `'constant'`: `const_albedo`; `'step'`: `higher_albedo`
  !! where the latitude is above `step_lat_nh` or below `step_lat_sh`; `'power'`: increase as
  !! `(|lat|/90)^albedo_exp`; `'tanh'`: tanh increase centred at `albedo_cntrNH`, `-albedo_cntrSH`
  !! with width `albedo_wdth`; `'sin2'`: `sin^2(lat)` increase from equator to pole;
  !! `'tanh_desert'`: as `'tanh'`, plus `albedo_desert` over the Sahara, Gobi and Australian
  !! deserts
  real :: const_albedo = 0.23  !! surface albedo (the low-latitude value of the profiles)
  real :: higher_albedo = 0.80  !! high-latitude albedo of the profiles
  real :: albedo_exp = 2.  !! exponent of the `'power'` profile
  real :: albedo_cntrNH = 68.
  !! [deg] centre latitude of the Northern Hemisphere albedo increase (`'tanh'`, `'tanh_desert'`)
  real :: albedo_cntrSH = 64.
  !! [deg] centre latitude (positive) of the Southern Hemisphere albedo increase (`'tanh'`,
  !! `'tanh_desert'`)
  real :: albedo_wdth = 5.  !! [deg] width of the albedo increase (`'tanh'`, `'tanh_desert'`)
  real :: step_lat_nh = 90.
  !! [deg] with `'step'`, `higher_albedo` where the latitude is above this value (90: none)
  real :: step_lat_sh = -90.
  !! [deg] with `'step'`, `higher_albedo` where the latitude is below this value (-90: none)
  real :: albedo_desert = 0.20
  !! with `'tanh_desert'`, the albedo over the desert boxes is `const_albedo + albedo_desert`

  namelist /surface_albedo_nml/ albedo_profile, const_albedo, higher_albedo, albedo_exp, &
    albedo_cntrNH, albedo_cntrSH, albedo_wdth, step_lat_nh, step_lat_sh, albedo_desert

contains

  !> Reads `surface_albedo_nml` and returns the background albedo on this PE's grid.
  subroutine surface_albedo_init(grid, albedo)
    type(surface_grid_type), intent(in) :: grid  !! the surface grid
    real, intent(out) :: albedo(:, :)  !! background surface albedo

    integer :: ierr, io, i, j
    real :: lat, lon, pi

    pi = 4.0*atan(1.)

    read (input_nml_file, nml=surface_albedo_nml, iostat=io)
    ierr = check_nml_error(io, 'surface_albedo_nml')
    if (mpp_pe() == mpp_root_pe()) write (stdlog(), nml=surface_albedo_nml)

    select case (trim(albedo_profile))
    case ('constant')
      albedo = const_albedo
    case ('step')
      do j = 1, size(albedo, 2)
        lat = 0.5*(grid%lat_bnd(j + 1) + grid%lat_bnd(j))*180/pi
        if (lat > step_lat_nh .or. lat < step_lat_sh) then
          albedo(:, j) = higher_albedo
        else
          albedo(:, j) = const_albedo
        end if
      end do
    case ('power')
      do j = 1, size(albedo, 2)
        lat = 0.5*(grid%lat_bnd(j + 1) + grid%lat_bnd(j))*180/pi
        lat = abs(lat)
        albedo(:, j) = const_albedo + (higher_albedo - const_albedo)*(lat/90.)**albedo_exp
      end do
    case ('tanh')
      do j = 1, size(albedo, 2)
        lat = 0.5*(grid%lat_bnd(j + 1) + grid%lat_bnd(j))*180/pi
        albedo(:, j) = const_albedo + &
                       (higher_albedo - const_albedo)*0.5*(1 + tanh((lat - albedo_cntrNH)/albedo_wdth)) + &
                       (higher_albedo - const_albedo)*0.5*(1 - tanh((lat + albedo_cntrSH)/albedo_wdth))
      end do
    case ('sin2')
      do j = 1, size(albedo, 2)
        lat = 0.5*(grid%lat_bnd(j + 1) + grid%lat_bnd(j))
        albedo(:, j) = const_albedo + (higher_albedo - const_albedo)* &
                       sin(lat)**2
      end do
    case ('tanh_desert')
      !cig: as 'tanh', but the Gobi, Sahara and Australian deserts have a higher albedo
      do j = 1, size(albedo, 2)
        lat = 0.5*(grid%lat_bnd(j + 1) + grid%lat_bnd(j))*180/pi

        albedo(:, j) = const_albedo + &
                       (higher_albedo - const_albedo)*0.5*(1 + tanh((lat - albedo_cntrNH)/albedo_wdth)) + &
                       (higher_albedo - const_albedo)*0.5*(1 - tanh((lat + albedo_cntrSH)/albedo_wdth))

        do i = 1, size(albedo, 1)
          lon = 0.5*(grid%lon_bnd(i + 1) + grid%lon_bnd(i))*180/pi

          if ((lon .gt. 118. .and. lon .lt. 145. .and. lat .gt. -30. .and. lat .lt. -19.) .or. &
              (lon .gt. 80. .and. lon .lt. 105. .and. lat .gt. 32. .and. lat .lt. 40.) .or. &
              (lon .gt. 80. .and. lon .lt. 115. .and. lat .gt. 40. .and. lat .lt. 52.) .or. &
              ((lon .gt. 345. .or. lon .lt. 50.) .and. lat .gt. 13. .and. lat .lt. 30.)) then
            albedo(i, j) = const_albedo + albedo_desert
          end if
        end do
      end do
    case default
      call error_mesg('surface_albedo_init', 'unknown albedo_profile = "'//trim(albedo_profile)// &
                      '"; use constant, step, power, tanh, sin2 or tanh_desert', FATAL)
    end select

  end subroutine surface_albedo_init

end module surface_albedo_mod
