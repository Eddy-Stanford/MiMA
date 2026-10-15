
!> Diffusivities in the planetary boundary layer from a non-local K-profile scheme.
!>
!> The boundary-layer depth is the height where the bulk Richardson number, computed from
!> the dry static energy and the wind relative to the lowest level, reaches
!> `rich_crit_pbl` (or the fixed depth `depth_0`). In the surface layer (the lowest
!> fraction `frac_inner` of the boundary layer) the diffusivities follow Monin-Obukhov
!> similarity (`mo_diff`); above it they decrease to zero at the boundary-layer top with a
!> quadratic profile. `molecular_diff` gives the molecular diffusivities.
!>
!> Namelist: `diffusivity_nml`
!> ([namelist reference](https://eddy-stanford.github.io/MiMA/Parameters/#diffusivity_nml)).
module diffusivity_mod

  use constants_mod, only: grav, vonkarm, cp_air, rdgas, rvgas

  use fms_mod, only: error_mesg, FATAL, &
                     check_nml_error, input_nml_file, &
                     mpp_pe, mpp_root_pe, &
                     write_version_number, stdlog

  use mima_monin_obukhov_mod, only: mo_diff

  implicit none
  private

! public interfaces
!=======================================================================

  public diffusivity, pbl_depth, molecular_diff

!--------------------- version number ----------------------------------

  character(len=128) :: version = '$Id: diffusivity.f90,v 10.0.6.1 2005/05/13 18:16:36 pjp Exp $'
  character(len=128) :: tagname = '$Name:  $'

!=======================================================================

!  DEFAULT VALUES OF NAMELIST PARAMETERS:

  logical :: fixed_depth = .false.  !! use a fixed boundary-layer depth `depth_0`
  real    :: depth_0 = 5000.0  !! [m] boundary-layer depth if `fixed_depth`
  real    :: frac_inner = 0.1
  !! fraction of the boundary layer that is the surface layer (between 0 and 1)
  real    :: rich_crit_pbl = 1.0  !! critical bulk Richardson number defining the boundary-layer top
  real    :: background_m = 0.0  !! [m2/s] minimum diffusivity for momentum
  real    :: background_t = 0.0  !! [m2/s] minimum diffusivity for heat

  namelist /diffusivity_nml/ fixed_depth, depth_0, frac_inner, &
    rich_crit_pbl, &
    background_m, background_t

!=======================================================================

!  OTHER MODULE VARIABLES

  real    :: small = 1.e-04
  real    :: gcp = grav/cp_air
  logical :: module_is_initialized = .false.
  real    :: beta = 1.458e-06
  real    :: rbop1 = 110.4
  real    :: rbop2 = 1.405

  real, parameter :: d608 = (rvgas - rdgas)/rdgas

contains

!=======================================================================

  !> Reads `diffusivity_nml` and checks its values.
  subroutine diffusivity_init

    integer :: unit, ierr, io

!------------------- read namelist input -------------------------------

    read (input_nml_file, nml=diffusivity_nml, iostat=io)
    ierr = check_nml_error(io, 'diffusivity_nml')

!------------------- dummy checks --------------------------------------
    if (frac_inner .le. 0. .or. frac_inner .ge. 1.) &
      call error_mesg('diffusivity_init', &
                      'frac_inner must be between 0 and 1', FATAL)
    if (rich_crit_pbl .lt. 0.) &
      call error_mesg('diffusivity_init', &
                      'rich_crit_pbl must be greater than or equal to zero', FATAL)
    if (background_m .lt. 0.) &
      call error_mesg('diffusivity_init', &
                      'background_m must be greater than or equal to zero', FATAL)
    if (background_t .lt. 0.) &
      call error_mesg('diffusivity_init', &
                      'background_t must be greater than or equal to zero', FATAL)

!---------- output namelist to log-------------------------------------

    if (mpp_pe() == mpp_root_pe()) then
      call write_version_number(version, tagname)
      write (stdlog(), nml=diffusivity_nml)
    end if

    module_is_initialized = .true.

    return
  end subroutine diffusivity_init

!=======================================================================

  subroutine diffusivity_end

    module_is_initialized = .false.

  end subroutine diffusivity_end

!=======================================================================

  !> Computes the boundary-layer depth and adds the boundary-layer diffusivities for
  !> momentum and heat to `k_m` and `k_t`.
  !>
  !> The diffusivities are defined at half levels; only levels 2 to `size(t,3)` are
  !> meaningful. `background_m` and `background_t` are then imposed as minimum values.
  subroutine diffusivity(t, q, u, v, p_full, p_half, z_full, z_half, &
                         u_star, b_star, h, k_m, k_t, kbot)

    real, intent(in), dimension(:, :, :) :: t, q, u, v
    !! `t`: temperature [K]; `q`: specific humidity [kg/kg] (not used); `u`, `v`: zonal and
    !! meridional wind [m/s]; the third index runs from the top of the atmosphere down
    real, intent(in), dimension(:, :, :) :: p_full, p_half
    !! pressure at full and half levels [Pa] (not used)
    real, intent(in), dimension(:, :, :) :: z_full, z_half
    !! height of full and half levels [m]; `z_half(:,:,size(z_half,3))` must be the height of
    !! the surface (if `kbot` is not present)
    real, intent(in), dimension(:, :)   :: u_star, b_star
    !! friction velocity [m/s] and buoyancy scale [m/s2] (from `mo_drag`)
    real, intent(inout), dimension(:, :, :) :: k_m, k_t
    !! diffusivities for momentum and for temperature and scalars at half levels [m2/s]
    real, intent(out), dimension(:, :)   :: h  !! depth of the planetary boundary layer [m]
    integer, intent(in), optional, dimension(:, :)   :: kbot  !! index of the lowest model level

    real, dimension(size(t, 1), size(t, 2), size(t, 3))  :: svcp, z_full_ag, &
                                                            k_m_save, k_t_save
    real, dimension(size(t, 1), size(t, 2), size(t, 3) + 1):: z_half_ag
    real, dimension(size(t, 1), size(t, 2))            :: z_surf
    integer                                         :: i, j, k, nlev, nlat, nlon

    if (.not. module_is_initialized) call diffusivity_init

    nlev = size(t, 3)

    k_m_save = k_m
    k_t_save = k_t

!compute height of surface
    if (present(kbot)) then
      nlat = size(t, 2)
      nlon = size(t, 1)
      do j = 1, nlat
      do i = 1, nlon
        z_surf(i, j) = z_half(i, j, kbot(i, j) + 1)
      end do
      end do
    else
      z_surf(:, :) = z_half(:, :, nlev + 1)
    end if

!compute density profile, and heights relative to surface
    do k = 1, nlev
      z_full_ag(:, :, k) = z_full(:, :, k) - z_surf(:, :)
      z_half_ag(:, :, k) = z_half(:, :, k) - z_surf(:, :)
      svcp(:, :, k) = t(:, :, k) + gcp*(z_full_ag(:, :, k))
    end do
    z_half_ag(:, :, nlev + 1) = z_half(:, :, nlev + 1) - z_surf(:, :)

    if (fixed_depth) then
      h = depth_0
    else
      call pbl_depth(svcp, u, v, z_full_ag, u_star, b_star, h, kbot=kbot)
    end if

    call diffusivity_pbl(svcp, u, v, z_half_ag, h, u_star, b_star, &
                         k_m, k_t, kbot=kbot)

    k_m = k_m + k_m_save
    k_t = k_t + k_t_save

!set background diffusivities
    if (background_m .gt. 0.0) k_m = max(k_m, background_m)
    if (background_t .gt. 0.0) k_t = max(k_t, background_t)

    return
  end subroutine diffusivity

!=======================================================================

  !> Computes the boundary-layer depth as the height where the bulk Richardson number
  !> relative to the lowest level first exceeds `rich_crit_pbl` (interpolated linearly
  !> between levels).
  subroutine pbl_depth(t, u, v, z, u_star, b_star, h, kbot)

    real, intent(in), dimension(:, :, :) :: t, u, v, z
    !! `t`: temperature variable of the Richardson number (dry static energy divided by
    !! `cp_air` when called from `diffusivity`) [K]; `u`, `v`: zonal and meridional wind [m/s];
    !! `z`: height of full levels above the surface [m]
    real, intent(in), dimension(:, :)   :: u_star, b_star
    !! friction velocity [m/s] and buoyancy scale [m/s2] (not used)
    real, intent(out), dimension(:, :)   :: h  !! depth of the planetary boundary layer [m]
    integer, intent(in), optional, dimension(:, :)   :: kbot  !! index of the lowest model level

    real, dimension(size(t, 1), size(t, 2), size(t, 3))  :: rich
    real, dimension(size(t, 1), size(t, 2))            :: tbot
    real                                               :: rich1, rich2, &
                                                          h1, h2
    integer, dimension(size(t, 1), size(t, 2))            :: ibot
    integer                                            :: i, j, k, nlon, &
                                                          nlat, nlev

    nlev = size(t, 3)
    nlat = size(t, 2)
    nlon = size(t, 1)

!assign ibot, compute tbot (virtual temperature at lowest level)
    if (present(kbot)) then
      ibot(:, :) = kbot
      do j = 1, nlat
      do i = 1, nlon
        tbot(i, j) = t(i, j, ibot(i, j))
      end do
      end do
    else
      ibot(:, :) = nlev
      tbot(:, :) = t(:, :, nlev)
    end if

!compute richardson number for use in pbl depth of neutral/stable side
    do k = 1, nlev
      rich(:, :, k) = z(:, :, k)*grav*(t(:, :, k) - tbot(:, :))/tbot(:, :) &
                      /(u(:, :, k)*u(:, :, k) + v(:, :, k)*v(:, :, k) + small)
    end do

    do j = 1, nlat
      do i = 1, nlon

        !neutral/stable Richardson-number method in all columns

        h1 = z(i, j, ibot(i, j))
        h(i, j) = h1
        rich1 = rich(i, j, ibot(i, j))
        do k = ibot(i, j) - 1, 1, -1
          rich2 = rich(i, j, k)
          h2 = z(i, j, k)
          if (rich2 .gt. rich_crit_pbl) then
            h(i, j) = h2 + (h1 - h2)*(rich2 - rich_crit_pbl) &
                      /(rich2 - rich1)
            go to 10
          end if
          rich1 = rich2
          h1 = h2
        end do

10      continue
      end do
    end do

    return
  end subroutine pbl_depth

!=======================================================================

  !> Computes the K-profile diffusivities in the boundary layer.
  subroutine diffusivity_pbl(t, u, v, z_half, h, u_star, b_star, &
                             k_m, k_t, kbot)

    real, intent(in), dimension(:, :, :) :: t, u, v, z_half
    real, intent(in), dimension(:, :)   :: h, u_star, b_star
    real, intent(inout), dimension(:, :, :) :: k_m, k_t
    integer, intent(in), optional, dimension(:, :)   :: kbot

    real, dimension(size(t, 1), size(t, 2))              :: h_inner, k_m_ref, &
                                                            k_t_ref, factor
    real, dimension(size(t, 1), size(t, 2), size(t, 3) + 1)  :: zm
    real                                              :: h_inner_max
    integer                                           :: i, j, k, kk, nlev

    nlev = size(t, 3)

!assign z_half to zm, and set to zero any values of zm < 0.
!the setting to zero is necessary so that when using eta model
!below ground half levels will have zero k_m and k_t
    zm = z_half
    if (present(kbot)) then
      where (zm < 0.)
        zm = 0.
      end where
    end if

    h_inner = frac_inner*h
    h_inner_max = maxval(h_inner)

    kk = nlev
    do k = 2, nlev
      if (minval(zm(:, :, k)) < h_inner_max) then
        kk = k
        exit
      end if
    end do

    k_m = 0.0
    k_t = 0.0

    call mo_diff(h_inner, u_star, b_star, k_m_ref, k_t_ref)
    call mo_diff(zm(:, :, kk:nlev), u_star, b_star, k_m(:, :, kk:nlev), k_t(:, :, kk:nlev))

    do k = 2, nlev
      where (zm(:, :, k) >= h_inner .and. zm(:, :, k) < h)
        factor = (zm(:, :, k)/h_inner)* &
                 (1.0 - (zm(:, :, k) - h_inner)/(h - h_inner))**2
        k_m(:, :, k) = k_m_ref*factor
        k_t(:, :, k) = k_t_ref*factor
      end where

! POG change: avoid possibility of k_m and k_t set to non-zero values above PBL due to use of maxval(h_inner) above
      where (zm(:, :, k) >= h)
        k_m(:, :, k) = 0
        k_t(:, :, k) = 0
      end where
! end POG change

    end do

    return
  end subroutine diffusivity_pbl

!=======================================================================

  !> Computes the molecular diffusivities for momentum and heat at half levels.
  !>
  !> `k_m` is the kinematic viscosity of air from Sutherland's formula; `k_t` is 1.405
  !> times `k_m`. Both are zero at the top half level.
  subroutine molecular_diff(temp, press, k_m, k_t)

    real, intent(in), dimension(:, :, :)  ::  temp, press
    !! temperature at full levels [K] and pressure at half levels [Pa]
    real, intent(inout), dimension(:, :, :)  ::  k_m, k_t
    !! molecular diffusivities for momentum and heat at half levels [m2/s]

    real, dimension(size(temp, 1), size(temp, 2)) :: temp_half, &
                                                     rho_half, rbop2d
    integer      :: k

!---------------------------------------------------------------------

    do k = 2, size(temp, 3)
      temp_half(:, :) = 0.5*(temp(:, :, k) + temp(:, :, k - 1))
      rho_half(:, :) = press(:, :, k)/(rdgas*temp_half(:, :))
      rbop2d(:, :) = beta*temp_half(:, :)*sqrt(temp_half(:, :))/ &
                     (rho_half(:, :)*(temp_half(:, :) + rbop1))
      k_m(:, :, k) = rbop2d(:, :)
      k_t(:, :, k) = rbop2d(:, :)*rbop2
    end do

    k_m(:, :, 1) = 0.0
    k_t(:, :, 1) = 0.0

  end subroutine molecular_diff

!=======================================================================

end module diffusivity_mod
