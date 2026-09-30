!> Pressures and geopotential heights on the model levels.
!>
!> The vertical coordinate has half-level pressures `p_half(k) = pk(k) + bk(k)*ps`, where
!> `pk` and `bk` define the levels and `ps` is the surface pressure. The full-level values
!> follow Simmons and Burridge (1981):
!> `ln(p_full(k)) = ln(p_half(k+1)) - alpha`, with
!> `alpha = 1 - p_half(k)*(ln(p_half(k+1)) - ln(p_half(k)))/(p_half(k+1) - p_half(k))`.
!> If the top half level is at zero pressure, `ln(p_full(1)) = ln(p_half(2)) - 1`.
!> The geopotential is computed by integrating the hydrostatic and ideal-gas equations
!> exactly with the temperature (or virtual temperature) constant in each layer.
!>
!> References:
!>
!> * Simmons, A. J., and D. M. Burridge, 1981: An energy and angular-momentum conserving
!>   vertical finite-difference scheme and hybrid vertical coordinates. Mon. Wea. Rev.,
!>   109, 758-766.
module press_and_geopot_mod

  use fms_mod, only: mpp_pe, mpp_root_pe, error_mesg, FATAL, &
                     write_version_number

  use constants_mod, only: grav, rdgas, rvgas

  implicit none

  private

  public :: press_and_geopot_init, press_and_geopot_end, pressure_variables, half_level_pressures
  public :: compute_geopotential, compute_pressures_and_heights, compute_z_bot

  !> Returns the half-level pressures `pk + bk*surface_p` [Pa] for a surface pressure field
  !> or a single value.
  interface half_level_pressures
    module procedure half_level_pressures_1d, &
      half_level_pressures_3d
  end interface

  !> Computes the pressures and log pressures at the half and full levels [Pa] for a surface
  !> pressure field or a single value.
  interface pressure_variables
    module procedure pressure_variables_1d, &
      pressure_variables_3d
  end interface

!===============================================================================================

  character(len=128), parameter :: version = '$Id: press_and_geopot.f90,v 11.0 2004/09/28 19:29:51 fms Exp $'
  character(len=128), parameter :: tagname = '$Name: lima $'

!===============================================================================================

  real, allocatable, dimension(:) :: pk, bk
  real, allocatable, dimension(:, :) :: surf_geopotential

  real    :: ln_top_level_factor
  integer :: num_levels
  logical :: use_virtual_temperature
  character(len=64) :: vert_difference_option

  logical :: module_is_initialized = .false.

contains

!------------------------------------------------------------------------------

  !> Stores the vertical coordinate, the surface geopotential and the options.
  subroutine press_and_geopot_init(pk_in, bk_in, use_virtual_temperature_in, vert_difference_option_in, surf_geopotential_in)

    real, intent(in), dimension(:)   :: pk_in, bk_in
    !! `pk_in` [Pa], `bk_in`: coefficients of the half-level pressures `pk + bk*ps`
    logical, intent(in)                 :: use_virtual_temperature_in  !! use virtual temperature in the geopotential
    character(len=*), intent(in)        :: vert_difference_option_in  !! vertical differencing (`'simmons_and_burridge'`)
    real, intent(in), dimension(:, :) :: surf_geopotential_in  !! surface geopotential [m2/s2]

    integer :: k

    if (module_is_initialized) return

    call write_version_number(version, tagname)

    num_levels = size(pk_in, 1) - 1

    allocate (pk(num_levels + 1))
    allocate (bk(num_levels + 1))
    pk = pk_in
    bk = bk_in

    vert_difference_option = vert_difference_option_in
    use_virtual_temperature = use_virtual_temperature_in

    allocate (surf_geopotential(size(surf_geopotential_in, 1), size(surf_geopotential_in, 2)))
    surf_geopotential = surf_geopotential_in

    ln_top_level_factor = -1.0

    module_is_initialized = .true.

    return
  end subroutine press_and_geopot_init

!-----------------------------------------------------------------------

  function half_level_pressures_3d(surface_p) result(p_half)

    real, intent(in), dimension(:, :)     :: surface_p
    real, dimension(size(surface_p, 1), size(surface_p, 2), num_levels + 1) :: p_half

    integer :: k

    if (.not. module_is_initialized) then
      call error_mesg('half_level_pressures', 'press_and_geopot_init has not been called', FATAL)
    end if

    do k = 1, num_levels + 1
      p_half(:, :, k) = pk(k) + bk(k)*surface_p(:, :)
    end do

    return
  end function half_level_pressures_3d

!-----------------------------------------------------------------------

  function half_level_pressures_1d(surface_p) result(p_half)

    real, intent(in)             :: surface_p
    real, dimension(num_levels + 1) :: p_half

    integer :: k

    do k = 1, num_levels + 1
      p_half(k) = pk(k) + bk(k)*surface_p
    end do

    return
  end function half_level_pressures_1d

!-------------------------------------------------------------------------------------

  subroutine pressure_variables_3d(p_half, ln_p_half, p_full, ln_p_full, surface_p)

    real, intent(out), dimension(:, :, :) :: p_half, ln_p_half, p_full, ln_p_full
    real, intent(in), dimension(:, :)   :: surface_p

    real, dimension(size(p_half, 1), size(p_half, 2)) :: alpha

    integer :: k

    if (.not. module_is_initialized) then
      call error_mesg('pressure_variables', 'press_and_geopot_init has not been called', FATAL)
    end if

    p_half = half_level_pressures(surface_p)

    if (trim(vert_difference_option) == 'simmons_and_burridge') then

      if (pk(1) .eq. 0.0 .and. bk(1) .eq. 0.0) then

        do k = 2, size(p_half, 3)
          ln_p_half(:, :, k) = log(p_half(:, :, k))
        end do

        do k = 2, size(p_half, 3) - 1
          alpha = 1.0 - p_half(:, :, k)*(ln_p_half(:, :, k + 1) - ln_p_half(:, :, k))/(p_half(:, :, k + 1) - p_half(:, :, k))
          ln_p_full(:, :, k) = ln_p_half(:, :, k + 1) - alpha
        end do
        ln_p_full(:, :, 1) = ln_p_half(:, :, 2) + ln_top_level_factor
        ln_p_half(:, :, 1) = 0.0

      else

        do k = 1, size(p_half, 3)
          ln_p_half(:, :, k) = log(p_half(:, :, k))
        end do

        do k = 1, size(p_half, 3) - 1
          alpha = 1.0 - p_half(:, :, k)*(ln_p_half(:, :, k + 1) - ln_p_half(:, :, k))/(p_half(:, :, k + 1) - p_half(:, :, k))
          ln_p_full(:, :, k) = ln_p_half(:, :, k + 1) - alpha
        end do

      end if
      p_full = exp(ln_p_full)

    else

      call error_mesg('pressure_variables', '"'//trim(vert_difference_option)//'"'// &
                      ' is not a valid value for vert_difference_option', FATAL)

    end if

    return
  end subroutine pressure_variables_3d

!-------------------------------------------------------------------------------------
  !> Computes the height of the lowest full level above the surface.
  subroutine compute_z_bot(psg, tg, z_bot, qg)
    real, intent(in), dimension(:, :) :: psg, tg
    !! `psg`: surface pressure [Pa]; `tg`: temperature at the lowest level [K]
    real, intent(out), dimension(:, :) :: z_bot  !! height of the lowest full level above the surface [m]
    real, intent(in), optional, dimension(:, :) :: qg
    !! specific humidity at the lowest level [kg/kg] (required with virtual temperature)

    real, dimension(size(psg, 1), size(psg, 2)) ::    p_half_bot, p_half_nxt
    real, dimension(size(psg, 1), size(psg, 2)) :: ln_p_half_bot, ln_p_half_nxt
    real, dimension(size(psg, 1), size(psg, 2)) :: ln_p_full_bot, alpha, virtual_t

    num_levels = size(pk, 1) - 1
    p_half_bot = pk(num_levels + 1) + psg*bk(num_levels + 1)
    p_half_nxt = pk(num_levels) + psg*bk(num_levels)
    ln_p_half_bot = log(p_half_bot)

    if (trim(vert_difference_option) == 'simmons_and_burridge') then
      ln_p_half_nxt = log(p_half_nxt)
      alpha = 1.0 - p_half_nxt*(ln_p_half_bot - ln_p_half_nxt)/(p_half_bot - p_half_nxt)
      ln_p_full_bot = ln_p_half_bot - alpha
    end if

    if (use_virtual_temperature) then
      if (present(qg)) then
        virtual_t = tg*(1.+(rvgas/rdgas - 1.)*qg)
      else
        call error_mesg('compute_z_bot', 'qg must be present when use_virtual_temperature=.true.', FATAL)
      end if
    else
      virtual_t = tg
    end if

    z_bot = (rdgas*virtual_t*(ln_p_half_bot - ln_p_full_bot))/grav

    return
  end subroutine compute_z_bot
!-------------------------------------------------------------------------------------

  subroutine pressure_variables_1d(p_half, ln_p_half, p_full, ln_p_full, surface_p)

    real, intent(out), dimension(:) :: p_half, ln_p_half, p_full, ln_p_full
    real, intent(in)   :: surface_p

    real, dimension(1, 1)              :: surface_p_2d
    real, dimension(1, 1, num_levels + 1) :: p_half_3d
    real, dimension(1, 1, num_levels + 1) :: ln_p_half_3d
    real, dimension(1, 1, num_levels)   :: p_full_3d
    real, dimension(1, 1, num_levels)   :: ln_p_full_3d

    surface_p_2d(1, 1) = surface_p

    call pressure_variables_3d(p_half_3d, ln_p_half_3d, p_full_3d, ln_p_full_3d, surface_p_2d)

    p_half = p_half_3d(1, 1, :)
    ln_p_half = ln_p_half_3d(1, 1, :)
    p_full = p_full_3d(1, 1, :)
    ln_p_full = ln_p_full_3d(1, 1, :)

    return
  end subroutine pressure_variables_1d

!-----------------------------------------------------------------------

  !> Computes the geopotential at the full and half levels, upward from the surface
  !> geopotential.
  subroutine compute_geopotential(t_grid, ln_p_half, ln_p_full, geopot_full, geopot_half, q_grid)

    real, intent(in), dimension(:, :, :) :: t_grid, ln_p_half, ln_p_full
    !! `t_grid`: temperature [K]; `ln_p_half`, `ln_p_full`: log of the half- and full-level pressures
    real, intent(out), dimension(:, :, :) :: geopot_full, geopot_half
    !! `geopot_full`, `geopot_half`: geopotential at the full and half levels [m2/s2]
    real, intent(in), optional, dimension(:, :, :) :: q_grid
    !! specific humidity [kg/kg] (required with virtual temperature)

    real, dimension(size(t_grid, 1), size(t_grid, 2), size(t_grid, 3)) :: virtual_t

    integer :: ktop, num_levels, k

    if (.not. module_is_initialized) then
      call error_mesg('compute_geopotential', 'press_and_geopot_init has not been called', FATAL)
    end if

    num_levels = size(t_grid, 3)

    geopot_half(:, :, num_levels + 1) = surf_geopotential

    if (pk(1) .eq. 0.0) then
      ktop = 2
      geopot_half(:, :, 1) = 0.0
    else
      ktop = 1
    end if

    if (use_virtual_temperature) then
      if (present(q_grid)) then
        virtual_t = t_grid*(1.+(rvgas/rdgas - 1.)*q_grid)
      else
        call error_mesg('compute_geopotential', 'q_grid must be present when use_virtual_temperature=.true.', FATAL)
      end if
    else
      virtual_t = t_grid
    end if

    do k = num_levels, ktop, -1
      geopot_half(:, :, k) = geopot_half(:, :, k + 1) + rdgas*virtual_t(:, :, k)*(ln_p_half(:, :, k + 1) - ln_p_half(:, :, k))
    end do

    do k = 1, num_levels
      geopot_full(:, :, k) = geopot_half(:, :, k + 1) + rdgas*virtual_t(:, :, k)*(ln_p_half(:, :, k + 1) - ln_p_full(:, :, k))
    end do

    return
  end subroutine compute_geopotential

!-----------------------------------------------------------------------

  !> Computes the pressures and the geopotential heights at the full and half levels.
  subroutine compute_pressures_and_heights(t_grid, ps_grid, z_full, z_half, p_full, p_half, q_grid)

    real, intent(in), dimension(:, :, :) :: t_grid  !! temperature [K]
    real, intent(in), dimension(:, :) :: ps_grid  !! surface pressure [Pa]
    real, intent(in), optional, dimension(:, :, :) :: q_grid
    !! specific humidity [kg/kg] (required with virtual temperature)

    real, intent(out), dimension(size(t_grid, 1), size(t_grid, 2), size(t_grid, 3)) :: z_full, p_full
    !! `z_full`: geopotential height [m] and `p_full`: pressure [Pa] at the full levels
    real, intent(out), dimension(size(t_grid, 1), size(t_grid, 2), size(t_grid, 3) + 1) :: z_half, p_half
    !! `z_half`: geopotential height [m] and `p_half`: pressure [Pa] at the half levels

    real, dimension(size(t_grid, 1), size(t_grid, 2), size(t_grid, 3)) :: ln_p_full
    real, dimension(size(t_grid, 1), size(t_grid, 2), size(t_grid, 3) + 1) :: ln_p_half

    if (.not. module_is_initialized) then
      call error_mesg('compute_pressures_and_heights', 'press_and_geopot_init has not been called', FATAL)
    end if

    call pressure_variables(p_half, ln_p_half, p_full, ln_p_full, ps_grid)

    call compute_geopotential(t_grid, ln_p_half, ln_p_full, z_full, z_half, q_grid)

    z_full = z_full/grav
    z_half = z_half/grav

    return
  end subroutine compute_pressures_and_heights

!-----------------------------------------------------------------------
  !> Deallocates the module arrays.
  subroutine press_and_geopot_end

    if (.not. module_is_initialized) return

    deallocate (pk, bk, surf_geopotential)
    module_is_initialized = .false.

    return
  end subroutine press_and_geopot_end
!-----------------------------------------------------------------------

end module press_and_geopot_mod
