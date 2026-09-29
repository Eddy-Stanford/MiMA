
module vert_turb_driver_mod

!-----------------------------------------------------------------------
!
!       driver for computing vertical diffusion coefficients
!
!         - non-local K scheme (diffusivity_mod), optionally with
!           molecular diffusion
!         - surface-layer gustiness
!
!-----------------------------------------------------------------------
!---------------- modules ---------------------


use    diffusivity_mod, only: diffusivity, molecular_diff

use   diag_manager_mod, only: register_diag_field, send_data

use   time_manager_mod, only: time_type

use       fms_mod,      only: mpp_pe, mpp_root_pe, stdlog, &
                              error_mesg, open_namelist_file, file_exist, &
                              check_nml_error, close_file, FATAL, &
                              write_version_number

implicit none
private

!---------------- interfaces ---------------------

public   vert_turb_driver_init, vert_turb_driver_end, vert_turb_driver


!-----------------------------------------------------------------------
!--------------------- version number ----------------------------------

character(len=128) :: version = '$Id: vert_turb_driver.f90,v 11.0.8.1 2005/05/13 18:16:38 pjp Exp $'
character(len=128) :: tagname = '$Name:  $'
logical            :: module_is_initialized = .false.

!---------------- private data -------------------

 real :: gust_zi = 1000.   ! constant for computed gustiness (meters)

!-----------------------------------------------------------------------
!-------------------- namelist -----------------------------------------

 logical :: do_diffusivity         = .false.
 logical :: do_molecular_diffusion = .false.
 logical :: use_tau          = .true.

 character(len=24) :: gust_scheme  = 'constant' ! valid schemes are:
                                                !   => 'constant'
                                                !   => 'beljaars'
 real              :: constant_gust = 1.0
 real              :: gust_factor   = 1.0

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

subroutine vert_turb_driver (is, js, Time_next, dt,                    &
                             p_half, p_full, z_half, z_full, u_star,   &
                             b_star, u, v, t, q, um, vm, tm, qm,       &
                             udt, vdt, tdt, qdt, diff_t, diff_m,       &
                             gust, z_pbl, mask, kbot                   )

!-----------------------------------------------------------------------
integer,         intent(in)         :: is, js
type(time_type), intent(in)         :: Time_next
   real,         intent(in)         :: dt
   real, intent(in), dimension(:,:) :: u_star, b_star
   real, intent(in), dimension(:,:,:) :: p_half, p_full, &
                                         z_half, z_full, &
                                         u, v, t, q, um, vm, tm, qm, &
                                         udt, vdt, tdt, qdt
   real, intent(out),   dimension(:,:,:) :: diff_t, diff_m
   real, intent(out),   dimension(:,:)   :: gust, z_pbl
   real, intent(in),optional, dimension(:,:,:) :: mask
integer, intent(in),optional, dimension(:,:) :: kbot
!-----------------------------------------------------------------------
logical, dimension(size(t,1),size(t,2),size(t,3)+1) :: lmask
real   , dimension(size(t,1),size(t,2),size(t,3)+1) :: diag3
real   , dimension(size(t,1),size(t,2),size(t,3))   :: tt, qq, uu, vv
integer :: nlev
logical :: used
!-----------------------------------------------------------------------
!----------------------- vertical turbulence ---------------------------
!-----------------------------------------------------------------------

      if (.not. module_is_initialized)  call error_mesg  &
                     ('vert_turb_driver in vert_turb_driver_mod',  &
                      'initialization has not been called', FATAL)

     nlev = size(p_full,3)

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
      endif

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
      call molecular_diff (tt, p_half, diff_m, diff_t)
    else
      diff_m = 0.0
      diff_t = 0.0
    endif

!---------------------------
!------------------- non-local K scheme --------------


    call diffusivity ( tt, qq, uu, vv, p_full, p_half, z_full, z_half,   &
                       u_star, b_star, z_pbl, diff_m, diff_t, &
                       kbot = kbot)

 endif

!-----------------------------------------------------------------------
!------------- define gustiness ------------

     if ( trim(gust_scheme) == 'constant' ) then
          gust = constant_gust
     else if ( trim(gust_scheme) == 'beljaars' ) then
!    --- from Beljaars (1994) and Beljaars and Viterbo (1999) ---
          where (b_star > 0.)
             gust = gust_factor * (u_star*b_star*gust_zi)**(1./3.)
          elsewhere
             gust = 0.
          endwhere
     endif

!-----------------------------------------------------------------------
!------------------------ diagnostics section --------------------------

!------- boundary layer depth -------
      if ( id_z_pbl > 0 ) then
         used = send_data ( id_z_pbl, z_pbl, Time_next, is, js )
      endif

!------- gustiness -------
      if ( id_gust > 0 ) then
         used = send_data ( id_gust, gust, Time_next, is, js )
      endif


!------- output diffusion coefficients ---------

  if ( id_diff_t > 0 .or. id_diff_m > 0 ) then
!       --- set up local mask for fields without surface data ---
        if (present(mask)) then
            lmask(:,:,1:nlev) = mask(:,:,1:nlev) > 0.5
            lmask(:,:,nlev+1) = .false.
        else
            lmask(:,:,1:nlev) = .true.
            lmask(:,:,nlev+1) = .false.
        endif
!       -- dummy data at surface --
        diag3(:,:,nlev+1)=0.0
  endif

!------- diffusion coefficient for heat/moisture -------
   if ( id_diff_t > 0 ) then
      diag3(:,:,1:nlev) = diff_t(:,:,1:nlev)
      used = send_data ( id_diff_t, diag3, Time_next, is, js, 1, mask=lmask )
   endif

!------- diffusion coefficient for momentum -------
   if ( id_diff_m > 0 ) then
      diag3(:,:,1:nlev) = diff_m(:,:,1:nlev)
      used = send_data ( id_diff_m, diag3, Time_next, is, js, 1, mask=lmask )
   endif

!--- geopotential height relative to the surface on full and half levels ----

   if ( id_z_half > 0 ) then
      !--- set up local mask for fields with surface data ---
      if ( present(mask) ) then
         lmask(:,:,1)        = .true.
         lmask(:,:,2:nlev+1) = mask(:,:,1:nlev) > 0.5
      else
         lmask = .true.
      endif
      used = send_data ( id_z_half, z_half, Time_next, is, js, 1, mask=lmask )
   endif

   if ( id_z_full > 0 ) then
      used = send_data ( id_z_full, z_full, Time_next, is, js, 1, rmask=mask)
   endif

!--- zonal and meridional wind on mass grid -------

   if ( id_uwnd > 0 ) then
      used = send_data ( id_uwnd, uu, Time_next, is, js, 1, rmask=mask)
   endif

   if ( id_vwnd > 0 ) then
      used = send_data ( id_vwnd, vv, Time_next, is, js, 1, rmask=mask)
   endif

!-----------------------------------------------------------------------

end subroutine vert_turb_driver

!#######################################################################

subroutine vert_turb_driver_init (axes, Time)

!-----------------------------------------------------------------------
   integer,         intent(in) :: axes(4)
   type(time_type), intent(in) :: Time
!-----------------------------------------------------------------------
   integer, dimension(3) :: full = (/1,2,3/), half = (/1,2,4/)
   integer :: ierr, unit, io

      if (module_is_initialized)  &
          call error_mesg  &
                   ('vert_turb_driver_init in vert_turb_driver_mod',  &
                    'attempting to call initialization twice', FATAL)

!-----------------------------------------------------------------------
!--------------- read namelist ------------------

      if (file_exist('input.nml')) then
         unit = open_namelist_file (file='input.nml')
         ierr=1; do while (ierr /= 0)
            read  (unit, nml=vert_turb_driver_nml, iostat=io, end=10)
            ierr = check_nml_error (io, 'vert_turb_driver_nml')
         enddo
  10     call close_file (unit)
      endif

!---------- output namelist --------------------------------------------

      if ( mpp_pe() == mpp_root_pe() ) then
           call write_version_number(version, tagname)
           write (stdlog(),nml=vert_turb_driver_nml)
      endif

!     --- check namelist option ---
      if ( trim(gust_scheme) /= 'constant' .and. &
           trim(gust_scheme) /= 'beljaars' ) call error_mesg &
         ('vert_turb_driver_mod', 'invalid value for namelist '//&
          'variable GUST_SCHEME', FATAL)

!-----------------------------------------------------------------------
!----- initialize diagnostic fields -----

   id_uwnd = register_diag_field ( mod_name, 'uwnd', axes(full), Time, &
        'zonal wind on mass grid', 'meters/second' ,                   &
         missing_value=missing_value    )

   id_vwnd = register_diag_field ( mod_name, 'vwnd', axes(full), Time, &
        'meridional wind on mass grid', 'meters/second' ,              &
        missing_value=missing_value    )

   id_z_full = &
   register_diag_field ( mod_name, 'z_full', axes(full), Time,    &
        'geopotential height relative to surface at full levels', &
         'meters' , missing_value=missing_value    )

   id_z_half = &
   register_diag_field ( mod_name, 'z_half', axes(half), Time,    &
        'geopotential height relative to surface at half levels', &
        'meters' , missing_value=missing_value    )

   id_z_pbl = &
   register_diag_field ( mod_name, 'z_pbl', axes(1:2), Time,       &
                        'depth of planetary boundary layer',  'm'  )

   id_gust = &
   register_diag_field ( mod_name, 'gust', axes(1:2), Time,        &
                        'wind gustiness in surface layer',  'm/s'  )

   id_diff_t = &
   register_diag_field ( mod_name, 'diff_t', axes(half), Time,    &
                        'vert diff coeff for temp',  'm2/s'   ,   &
                        missing_value=missing_value               )

   id_diff_m = &
   register_diag_field ( mod_name, 'diff_m', axes(half), Time,      &
                        'vert diff coeff for momentum',  'm2/s'   , &
                        missing_value=missing_value               )

!-----------------------------------------------------------------------

   module_is_initialized =.true.

!-----------------------------------------------------------------------

end subroutine vert_turb_driver_init


!#######################################################################

subroutine vert_turb_driver_end

!-----------------------------------------------------------------------
      module_is_initialized =.false.

!-----------------------------------------------------------------------

end subroutine vert_turb_driver_end

!#######################################################################

end module vert_turb_driver_mod

