
!> Driver for the vertical diffusion coefficients and the surface-layer gustiness.
!>
!> Computes the diffusion coefficients with the non-local K scheme of `diffusivity_mod`
!> (`do_diffusivity`), optionally with molecular diffusion added
!> (`do_molecular_diffusion`), and the gustiness used by the surface fluxes: a constant
!> or, with `gust_scheme = 'beljaars'`, computed from `u_star` and `b_star` after Beljaars (1994) and
!> Beljaars and Viterbo (1999). Sends the boundary-layer diagnostics.
!>
!> Namelist: `vert_turb_driver_nml`
!> ([namelist reference](https://eddy-stanford.github.io/MiMA/Parameters/#vert_turb_driver_nml)).
module vert_turb_driver_mod

!---------------- modules ---------------------

  use diffusivity_mod, only: diffusivity, molecular_diff

  use diag_manager_mod, only: register_diag_field, send_data

  use time_manager_mod, only: time_type

  use fms_mod, only: mpp_pe, mpp_root_pe, stdlog, &
                     error_mesg, input_nml_file, &
                     check_nml_error, FATAL, &
                     write_version_number

  implicit none
  private

!---------------- interfaces ---------------------

  public vert_turb_driver_init, vert_turb_driver_end, vert_turb_driver

!-----------------------------------------------------------------------
!--------------------- version number ----------------------------------

  character(len=128) :: version = '$Id: vert_turb_driver.f90,v 11.0.8.1 2005/05/13 18:16:38 pjp Exp $'
  character(len=128) :: tagname = '$Name:  $'
  logical            :: module_is_initialized = .false.

!---------------- private data -------------------

  real :: gust_zi = 1000.   ! constant for computed gustiness (meters)

!-----------------------------------------------------------------------
!-------------------- namelist -----------------------------------------

  logical :: do_diffusivity = .true.
  !! compute diffusion coefficients with the non-local K scheme (`.false.`: no
  !! boundary-layer diffusion)
  logical :: do_molecular_diffusion = .false.  !! add molecular diffusion
  logical :: use_tau = .false.
  !! use the current time level (`.true.`) or the updated values (`.false.`)

  character(len=24) :: gust_scheme = 'constant'
  !! surface gustiness: `'constant'` (`constant_gust`) or `'beljaars'` (from u* and b*)
  real              :: constant_gust = 0.  !! [m/s] constant gustiness
  real              :: gust_factor = 1.0  !! factor for the `'beljaars'` gustiness

  namelist /vert_turb_driver_nml/ gust_scheme, constant_gust, use_tau, &
    do_molecular_diffusion, &
    do_diffusivity, gust_factor

!-------------------- diagnostics fields -------------------------------

  integer :: id_z_pbl, id_gust, id_diff_t, id_diff_m, id_z_full, id_z_half, &
             id_uwnd, id_vwnd

  real :: missing_value = -999.

  character(len=9) :: mod_name = 'vert_turb'

!-----------------------------------------------------------------------

contains

!#######################################################################

  !> Computes the diffusion coefficients, the boundary-layer depth and the gustiness, and
  !> sends their diagnostics.
  subroutine vert_turb_driver(is, js, Time_next, dt, &
                              p_half, p_full, z_half, z_full, u_star, &
                              b_star, u, v, t, q, um, vm, tm, qm, &
                              udt, vdt, tdt, qdt, diff_t, diff_m, &
                              gust, z_pbl, mask, kbot)

!-----------------------------------------------------------------------
    integer, intent(in)         :: is, js  !! starting i,j indices of the physics window
    type(time_type), intent(in)         :: Time_next  !! time at the end of the step (for the diagnostics)
    real, intent(in)         :: dt  !! time step [s]
    real, intent(in), dimension(:, :) :: u_star, b_star  !! friction velocity [m/s] and buoyancy scale [m/s2]
    real, intent(in), dimension(:, :, :) :: p_half, p_full, &
                                            z_half, z_full, &
                                            u, v, t, q, um, vm, tm, qm, &
                                            udt, vdt, tdt, qdt
    !! `p_half`, `p_full`: pressure at half and full levels [Pa]; `z_half`, `z_full`: height
    !! of half and full levels [m]; `u`, `v`, `t`, `q`: zonal and meridional wind [m/s],
    !! temperature [K] and specific humidity [kg/kg] at the current time level; `um`, `vm`,
    !! `tm`, `qm`: the same at the previous time level; `udt`, `vdt`, `tdt`, `qdt`: their
    !! tendencies [m/s2], [K/s], [kg/kg/s]
    real, intent(out), dimension(:, :, :) :: diff_t, diff_m
    !! diffusion coefficients for heat and moisture and for momentum at half levels [m2/s]
    real, intent(out), dimension(:, :)   :: gust, z_pbl
    !! `gust`: surface-layer gustiness [m/s]; `z_pbl`: boundary-layer depth (-999 if not
    !! `do_diffusivity`) [m]
    real, intent(in), optional, dimension(:, :, :) :: mask  !! mask for the diagnostics
    integer, intent(in), optional, dimension(:, :) :: kbot  !! index of the lowest model level
!-----------------------------------------------------------------------
    logical, dimension(size(t, 1), size(t, 2), size(t, 3) + 1) :: lmask
    real, dimension(size(t, 1), size(t, 2), size(t, 3) + 1) :: diag3
    real, dimension(size(t, 1), size(t, 2), size(t, 3))   :: tt, qq, uu, vv
    integer :: nlev
    logical :: used
!-----------------------------------------------------------------------
!----------------------- vertical turbulence ---------------------------
!-----------------------------------------------------------------------

    if (.not. module_is_initialized) call error_mesg &
      ('vert_turb_driver in vert_turb_driver_mod', &
       'initialization has not been called', FATAL)

    nlev = size(p_full, 3)

!-----------------------------------------------------------------------
!---- set up state variable used by this module ----

    if (use_tau) then
      !-- variables at time tau
      uu = u
      vv = v
      tt = t
      qq = q
    else
      !-- variables at time tau+1
      uu = um + dt*udt
      vv = vm + dt*vdt
      tt = tm + dt*tdt
      qq = qm + dt*qdt
    end if

!--------------------------------------------------------------------
! initialize output

    diff_t = 0.0
    diff_m = 0.0
    z_pbl = -999.0

!-----------------------------------------------------------------------
    if (do_diffusivity) then
!--------------------------------------------------------------------
!----------- compute molecular diffusion, if desired  ---------------

      if (do_molecular_diffusion) then
        call molecular_diff(tt, p_half, diff_m, diff_t)
      else
        diff_m = 0.0
        diff_t = 0.0
      end if

!---------------------------
!------------------- non-local K scheme --------------

      call diffusivity(tt, qq, uu, vv, p_full, p_half, z_full, z_half, &
                       u_star, b_star, z_pbl, diff_m, diff_t, &
                       kbot=kbot)

    end if

!-----------------------------------------------------------------------
!------------- define gustiness ------------

    if (trim(gust_scheme) == 'constant') then
      gust = constant_gust
    else if (trim(gust_scheme) == 'beljaars') then
!    --- from Beljaars (1994) and Beljaars and Viterbo (1999) ---
      where (b_star > 0.)
        gust = gust_factor*(u_star*b_star*gust_zi)**(1./3.)
      elsewhere
        gust = 0.
      end where
    end if

!-----------------------------------------------------------------------
!------------------------ diagnostics section --------------------------

!------- boundary layer depth -------
    if (id_z_pbl > 0) then
      used = send_data(id_z_pbl, z_pbl, Time_next, is, js)
    end if

!------- gustiness -------
    if (id_gust > 0) then
      used = send_data(id_gust, gust, Time_next, is, js)
    end if

!------- output diffusion coefficients ---------

    if (id_diff_t > 0 .or. id_diff_m > 0) then
!       --- set up local mask for fields without surface data ---
      if (present(mask)) then
        lmask(:, :, 1:nlev) = mask(:, :, 1:nlev) > 0.5
        lmask(:, :, nlev + 1) = .false.
      else
        lmask(:, :, 1:nlev) = .true.
        lmask(:, :, nlev + 1) = .false.
      end if
!       -- dummy data at surface --
      diag3(:, :, nlev + 1) = 0.0
    end if

!------- diffusion coefficient for heat/moisture -------
    if (id_diff_t > 0) then
      diag3(:, :, 1:nlev) = diff_t(:, :, 1:nlev)
      used = send_data(id_diff_t, diag3, Time_next, is, js, 1, mask=lmask)
    end if

!------- diffusion coefficient for momentum -------
    if (id_diff_m > 0) then
      diag3(:, :, 1:nlev) = diff_m(:, :, 1:nlev)
      used = send_data(id_diff_m, diag3, Time_next, is, js, 1, mask=lmask)
    end if

!--- geopotential height relative to the surface on full and half levels ----

    if (id_z_half > 0) then
      !--- set up local mask for fields with surface data ---
      if (present(mask)) then
        lmask(:, :, 1) = .true.
        lmask(:, :, 2:nlev + 1) = mask(:, :, 1:nlev) > 0.5
      else
        lmask = .true.
      end if
      used = send_data(id_z_half, z_half, Time_next, is, js, 1, mask=lmask)
    end if

    if (id_z_full > 0) then
      used = send_data(id_z_full, z_full, Time_next, is, js, 1, rmask=mask)
    end if

!--- zonal and meridional wind on mass grid -------

    if (id_uwnd > 0) then
      used = send_data(id_uwnd, uu, Time_next, is, js, 1, rmask=mask)
    end if

    if (id_vwnd > 0) then
      used = send_data(id_vwnd, vv, Time_next, is, js, 1, rmask=mask)
    end if

!-----------------------------------------------------------------------

  end subroutine vert_turb_driver

!#######################################################################

  !> Initializes the module: reads and checks `vert_turb_driver_nml` and registers the
  !> diagnostics.
  subroutine vert_turb_driver_init(axes, Time)

!-----------------------------------------------------------------------
    integer, intent(in) :: axes(4)  !! diagnostic axes (lon, lat, pfull, phalf)
    type(time_type), intent(in) :: Time  !! current time
!-----------------------------------------------------------------------
    integer, dimension(3) :: full = (/1, 2, 3/), half = (/1, 2, 4/)
    integer :: ierr, unit, io

    if (module_is_initialized) &
      call error_mesg &
      ('vert_turb_driver_init in vert_turb_driver_mod', &
       'attempting to call initialization twice', FATAL)

!-----------------------------------------------------------------------
!--------------- read namelist ------------------

    read (input_nml_file, nml=vert_turb_driver_nml, iostat=io)
    ierr = check_nml_error(io, 'vert_turb_driver_nml')

!---------- output namelist --------------------------------------------

    if (mpp_pe() == mpp_root_pe()) then
      call write_version_number(version, tagname)
      write (stdlog(), nml=vert_turb_driver_nml)
    end if

!     --- check namelist option ---
    if (trim(gust_scheme) /= 'constant' .and. &
        trim(gust_scheme) /= 'beljaars') call error_mesg &
      ('vert_turb_driver_mod', 'invalid value for namelist '// &
       'variable GUST_SCHEME', FATAL)

!-----------------------------------------------------------------------
!----- initialize diagnostic fields -----

    id_uwnd = register_diag_field(mod_name, 'uwnd', axes(full), Time, &
                                  'zonal wind on mass grid', 'm/s', &
                                  missing_value=missing_value)

    id_vwnd = register_diag_field(mod_name, 'vwnd', axes(full), Time, &
                                  'meridional wind on mass grid', 'm/s', &
                                  missing_value=missing_value)

    id_z_full = &
      register_diag_field(mod_name, 'z_full', axes(full), Time, &
                          'geopotential height relative to surface at full levels', &
                          'm', missing_value=missing_value)

    id_z_half = &
      register_diag_field(mod_name, 'z_half', axes(half), Time, &
                          'geopotential height relative to surface at half levels', &
                          'm', missing_value=missing_value)

    id_z_pbl = &
      register_diag_field(mod_name, 'z_pbl', axes(1:2), Time, &
                          'depth of planetary boundary layer', 'm')

    id_gust = &
      register_diag_field(mod_name, 'gust', axes(1:2), Time, &
                          'wind gustiness in surface layer', 'm/s')

    id_diff_t = &
      register_diag_field(mod_name, 'diff_t', axes(half), Time, &
                          'vert diff coeff for temp', 'm2/s', &
                          missing_value=missing_value)

    id_diff_m = &
      register_diag_field(mod_name, 'diff_m', axes(half), Time, &
                          'vert diff coeff for momentum', 'm2/s', &
                          missing_value=missing_value)

!-----------------------------------------------------------------------

    module_is_initialized = .true.

!-----------------------------------------------------------------------

  end subroutine vert_turb_driver_init

!#######################################################################

  !> Finalizes the module.
  subroutine vert_turb_driver_end

!-----------------------------------------------------------------------
    module_is_initialized = .false.

!-----------------------------------------------------------------------

  end subroutine vert_turb_driver_end

!#######################################################################

end module vert_turb_driver_mod

