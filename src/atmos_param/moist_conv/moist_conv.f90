
module moist_conv_mod

!-----------------------------------------------------------------------

  use mpp_mod, only: mpp_pe, &
                     mpp_root_pe, &
                     stdlog
  use time_manager_mod, only: time_type
  use Diag_Manager_Mod, only: register_diag_field, send_data
  use sat_vapor_pres_mod, only: EsComp, DEsComp
  use fms_mod, only: error_mesg, input_nml_file, &
                     check_nml_error, &
                     FATAL, WARNING, NOTE, mpp_pe, mpp_root_pe, &
                     write_version_number, stdlog
  use constants_mod, only: HLv, HLs, cp_air, grav, rdgas, rvgas

  use fms_mod, only: write_version_number, ERROR_MESG, FATAL
  use field_manager_mod, only: MODEL_ATMOS
  use tracer_manager_mod, only: get_number_tracers, &
                                get_tracer_names, &
                                get_tracer_indices, &
                                query_method, &
                                NO_TRACER

  implicit none
  private

!------- interfaces in this module ------------

  public :: moist_conv, moist_conv_Init, moist_conv_end

!-----------------------------------------------------------------------
!---- namelist ----

  real :: HC = 1.00
  real :: TOLmin = .02, TOLmax = .10
  integer :: ITSMOD = 30

  namelist /moist_conv_nml/ HC, TOLmin, TOLmax, ITSMOD

!-----------------------------------------------------------------------
!---- VERSION NUMBER -----

  character(len=128) :: version = '$Id: moist_conv.f90,v 11.0.2.1 2005/05/13 18:16:37 pjp Exp $'
  character(len=128) :: tagname = '$Name:  $'
  logical            :: module_is_initialized = .false.

!---------- initialize constants used by this module -------------------

  real, parameter :: d622 = rdgas/rvgas
  real, parameter :: d378 = 1.0 - d622
  real, parameter :: grav_inv = 1.0/grav
  real, parameter :: rocp = rdgas/cp_air

  real :: missing_value = -999.
!integer               :: num_tracers
!nteger, allocatable, dimension(:) :: id_tracer_conv, id_tracer_conv_col

  integer :: id_tdt_conv, id_qdt_conv, id_prec_conv, id_snow_conv, &
             id_q_conv_col, id_t_conv_col

  character(len=3) :: mod_name = 'mca'
  logical :: used

  logical :: do_mca_tracer = .false.
  integer :: num_mca_tracers = 0
  integer               :: num_tracers
  integer, allocatable, dimension(:) :: id_tracer_conv, id_tracer_conv_col

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

contains

!#######################################################################

  subroutine moist_conv(Tin, Qin, Pfull, Phalf, coldT, &
                        Tdel, Qdel, Rain, Snow, Lbot, &
                        dtinv, Time, mask, is, js, Conv, &
                        tracers, qtrmca)

!-----------------------------------------------------------------------
!
!                       MOIST CONVECTIVE ADJUSTMENT
!
!-----------------------------------------------------------------------
!
!   INPUT:   Tin     temperature at full model levels
!            Qin     specific humidity of water vapor at full
!                      model levels
!            Pfull   pressure at full model levels
!            Phalf   pressure at half model levels
!            coldT   Should MCA produce snow in this column?
!
!   OUTPUT:  Tdel    temperature adjustment at full model levels (deg k)
!            Qdel    specific humidity adjustment of water vapor at
!                       full model levels
!            Rain    liquid precipitiation (in Kg m-2)
!            Snow    ice phase precipitation (kg m-2)
!  OPTIONAL
!
!   INPUT:   Lbot    integer index of the lowest model level,
!                      Lbot is always <= size(Tin,3)
!
!  OUTPUT:   Conv    logical flag; TRUE then moist convective
!                       adjustment was performed at that model level.
!
!-----------------------------------------------------------------------
!----------------------PUBLIC INTERFACE ARRAYS--------------------------
    real, intent(INOUT), dimension(:, :, :)           :: Tin, Qin
    real, intent(IN), dimension(:, :, :)           :: Pfull, Phalf
    logical, intent(IN), dimension(:, :)             :: coldT
    real, intent(OUT), dimension(:, :, :)           :: Tdel, Qdel
    real, intent(OUT), dimension(:, :)             :: Rain, Snow
    integer, intent(IN), dimension(:, :), optional :: Lbot
    logical, intent(OUT), dimension(:, :, :), optional :: Conv
    real, intent(IN)                                :: dtinv
    integer, intent(IN)                                :: is, js
    real, dimension(:, :, :, :), intent(in), optional :: tracers
    real, dimension(:, :, :, :), intent(out), optional :: qtrmca
    type(time_type), intent(in)                         :: Time
    real, intent(in), dimension(:, :, :), optional :: mask

!-----------------------------------------------------------------------
!----------------------PRIVATE (LOCAL) ARRAYS---------------------------
! logical, dimension(size(Tin,1),size(Tin,2),size(Tin,3)) :: DO_ADJUST
!    real, dimension(size(Tin,1),size(Tin,2),size(Tin,3)) ::  &
!-----------------------------------------------------------------------
!----------------------PRIVATE (LOCAL) ARRAYS---------------------------
    integer, dimension(size(Tin, 1), size(Tin, 2)) :: ISMVD, ISMVF
    integer, dimension(size(Tin, 1), size(Tin, 2), size(Tin, 3)) :: IVF

    real, dimension(size(Tin, 1), size(Tin, 2), size(Tin, 3)) :: &
      Qdif, Temp, Qmix, Esat, Qsat, Test1, Test2

    real, dimension(size(Tin, 1), size(Tin, 2), size(Tin, 3) - 1) :: &
      Thalf, DelPoP, Esm, Esd, ALRM

    real, dimension(size(Tin, 3)) :: C, Ta, Qa
    real, dimension(size(Tin, 1), size(Tin, 2)) :: HL

    integer :: i, j, k, kk, KX, ITER, MXLEV, MXLEV1, kstart, KTOP, KBOT, KBOTM1
    real    :: ALTOL, Sum0, Sum1, Sum2, EsDiff, EsVal, Thaf, Pdelta

    real, dimension(size(Phalf, 1), size(Phalf, 2), size(Phalf, 3)) :: pmass
    real, dimension(size(Phalf, 1), size(Phalf, 2)) :: tempdiag
    integer  :: tr
!-----------------------------------------------------------------------

    if (.not. module_is_initialized) call ERROR_MESG('MCA', &
                                                     'moist_conv_init has not been called', FATAL)
    !moist_conv_init ( )

    do k = 1, size(Phalf, 3)
      pmass(:, :, k) = (Phalf(:, :, k + 1) - Phalf(:, :, k))/GRAV
    end do

    KX = size(Tin, 3)

!------ compute Proper HL
    HL = HLv

!------ convert spec hum to mixing ratio ------
    Temp(:, :, :) = Tin(:, :, :)
    Qmix(:, :, :) = Qin(:, :, :)

    do k = 1, KX - 1
      DelPoP(:, :, k) = (Pfull(:, :, k + 1) - Pfull(:, :, k))/Phalf(:, :, k + 1)
    end do

!-------------SATURATION VAPOR PRESSURE FROM ETABL----------------------

    call EsComp(Temp, Esat)

    Esat(:, :, :) = Esat(:, :, :)*HC
    Qsat(:, :, :) = Pfull(:, :, :)
    Qsat(:, :, :) = max(0.0, d622*Esat(:, :, :)/Qsat(:, :, :))
    Qdif(:, :, :) = max(0.0, Qmix(:, :, :) - Qsat(:, :, :))

!-----------------------------------------------------------------------
!                  MOIST CONVECTIVE ADJUSTMENT
!-----------------------------------------------------------------------

!  *** Set initial tolerance ***

    ALTOL = TOLmin

    do k = 1, KX - 1
      Thalf(:, :, k) = 0.50*(Temp(:, :, k) + Temp(:, :, k + 1))
    end do

    call EsComp(Thalf, Esm)
    call DEsComp(Thalf, Esd)

    do k = 1, KX - 1
      ALRM(:, :, k) = rocp*DelPoP(:, :, k)*Thalf(:, :, k) &
                      *(Phalf(:, :, k + 1) + d622*HL(:, :)*Esm(:, :, k)/Thalf(:, :, k)/rdgas) &
                      /(Phalf(:, :, k + 1) + d622*HL(:, :)*Esd(:, :, k)/cp_air)
    end do

    IVF(:, :, KX) = 0
    Test1(:, :, KX) = 0.0
    Test2(:, :, KX) = 0.0

    do k = 1, KX - 1
      Test1(:, :, k) = Temp(:, :, k + 1) - Temp(:, :, k)
      Test2(:, :, k) = ALRM(:, :, k) + ALTOL - Test1(:, :, k)
    end do

!!!!! Test1(:,:,:)=0.0-Qdif(:,:,:)
    Test1(:, :, :) = (0.0 - Qdif(:, :, :))*Qsat(:, :, :)

!-------IVF=1 in unstable layers where both levels are saturated--------

    do k = 1, KX - 1
      where (Test1(:, :, k) < 0.0 .and. Test1(:, :, k + 1) < 0.0 .and. &
             Test2(:, :, k) < 0.0)
        IVF(:, :, k) = 1
      elsewhere
        IVF(:, :, k) = 0
      end where
    end do

!  ------ Set convection flag (for optional output only) --------

    if (present(Conv)) then
      Conv(:, :, 1) = (IVF(:, :, 1) == 1)
      do k = 1, KX - 1
        Conv(:, :, k + 1) = (IVF(:, :, k) == 1 .or. IVF(:, :, k + 1) == 1)
      end do
    end if

!  ----- Set counter for each column -----

    ISMVF(:, :) = 0
    do k = 1, KX - 1
      ISMVF(:, :) = ISMVF(:, :) + IVF(:, :, k)
    end do

!-----------------------------------------------------------------------
!---------------LOOP OVER EACH VERTICAL COLUMN--------------------------
    do j = 1, size(Tin, 2)
      OUTER_LOOP: do i = 1, size(Tin, 1)
!-----------------------------------------------------------------------
        if (ISMVF(i, j) == 0) cycle

        if (present(Lbot)) then
          MXLEV = Lbot(i, j)
        else
          MXLEV = KX
        end if
        MXLEV1 = MXLEV - 1

!  *** Re-set initial tolerance ***
        ALTOL = TOLmin

!----------(return here after increasing tolerance)--------------------
1450    continue

!--------------Iterations at the same tolerance-------------------------
        do 1740 ITER = 1, ITSMOD
!-----------------------------------------------------------------------
          kstart = 1
1500      continue
!-------------TEST TO DETERMINE UNSTABLE LAYER BLOCKS-------------------
!-------Find top (KTOP) and bottom (KBOT) of unstable layers------------
          do k = kstart, MXLEV1
            if (IVF(i, j, k) == 1) GO TO 1505
          end do
          cycle OUTER_LOOP
1505      KTOP = k

          do k = KTOP, MXLEV1
            if (IVF(i, j, k + 1) == 0) then
              KBOT = k + 1
              GO TO 1510
            end if
          end do
          KBOT = MXLEV
1510      continue
!-----------------------------------------------------------------------

          KBOTM1 = KBOT - 1
          Sum1 = 0.0
          Sum2 = 0.0
!-----------------------------------------------------------------------
          do 1630 k = KTOP, KBOT
!-----------------------------------------------------------------------
            call DEsComp(Temp(i, j, k), EsDiff)
            C(k) = d622*HC*EsDiff/Pfull(i, j, k)

            Sum0 = 0.0
            if (k == KBOT) GO TO 1625
            kk = k
1620        if (kk > KBOTM1) GO TO 1625
            Sum0 = Sum0 + ALRM(i, j, kk)
            kk = kk + 1
            GO TO 1620

1625        continue
            Pdelta = Phalf(i, j, k + 1) - Phalf(i, j, k)
            Sum1 = Sum1 + Pdelta*((cp_air + HL(i, j)*C(k))*(Temp(i, j, k) + Sum0) + &
                                  HL(i, j)*(Qmix(i, j, k) - Qsat(i, j, k)))
            Sum2 = Sum2 + Pdelta*(cp_air + HL(i, j)*C(k))
!-----------------------------------------------------------------------
1630        continue
!-----------------------------------------------------------------------
            Ta(KBOT) = Sum1/Sum2
            k = KTOP
1645        if (k > KBOTM1) GO TO 1641
            Sum0 = 0.0
            kk = k
1640        if (kk > KBOTM1) GO TO 1642
            Sum0 = Sum0 + ALRM(i, j, kk)
            kk = kk + 1
            GO TO 1640
1642        Ta(k) = Ta(KBOT) - Sum0
            k = k + 1
            GO TO 1645

!---------UPDATE T,R,ES,Esm,Esd & Qsat FOR THE ADJUSTED POINTS----------

1641        do k = KTOP, KBOT
              Qa(k) = Qsat(i, j, k) + C(k)*(Ta(k) - Temp(i, j, k))
              Temp(i, j, k) = Ta(k)
              Qmix(i, j, k) = Qa(k)
!DIR$ INLINE
              call EsComp(Temp(i, j, k), EsVal)
!DIR$ NOINLINE
              Esat(i, j, k) = HC*EsVal
              Qsat(i, j, k) = Pfull(i, j, k)
              Qsat(i, j, k) = max(0.0, d622*Esat(i, j, k)/Qsat(i, j, k))
              Qdif(i, j, k) = max(0.0, Qmix(i, j, k) - Qsat(i, j, k))
            end do

            do k = KTOP, KBOTM1
              Thaf = 0.50*(Temp(i, j, k) + Temp(i, j, k + 1))
!DIR$ INLINE
              call EsComp(Thaf, EsVal)
              call DEsComp(Thaf, EsDiff)
!DIR$ NOINLINE
              Esm(i, j, k) = HC*EsVal
              Esd(i, j, k) = HC*EsDiff
              ALRM(i, j, k) = rocp*DelPoP(i, j, k)*Thaf* &
                              (Phalf(i, j, k + 1) + d622*HL(i, j)*Esm(i, j, k)/Thaf/rdgas)/ &
                              (Phalf(i, j, k + 1) + d622*HL(i, j)*Esd(i, j, k)/cp_air)
            end do

!------------Is this the bottom of the current column ???---------------
            kstart = KBOT + 1
            if (kstart <= MXLEV1) GO TO 1500
!-----------------------------------------------------------------------
            if (ITER == ITSMOD) GO TO 1740
!-----------------------------------------------------------------------

            do k = 1, MXLEV1
              Thaf = 0.50*(Temp(i, j, k) + Temp(i, j, k + 1))
!DIR$ INLINE
              call EsComp(Thaf, EsVal)
              call DEsComp(Thaf, EsDiff)
!DIR$ NOINLINE
              Esm(i, j, k) = HC*EsVal
              Esd(i, j, k) = HC*EsDiff
              ALRM(i, j, k) = rocp*DelPoP(i, j, k)*Thaf* &
                              (Phalf(i, j, k + 1) + d622*HL(i, j)*Esm(i, j, k)/Thaf/rdgas)/ &
                              (Phalf(i, j, k + 1) + d622*HL(i, j)*Esd(i, j, k)/cp_air)
            end do

            do k = 1, MXLEV1
              IVF(i, j, k) = 0
!!!!    if (Qdif(i,j,k) > 0.0 .and. Qdif(i,j,k+1) > 0.0 .and.  &
              if (Qdif(i, j, k)*Qsat(i, j, k) > 0.0 .and. &
                  Qdif(i, j, k + 1)*Qsat(i, j, k + 1) > 0.0 .and. &
                  (Temp(i, j, k + 1) - Temp(i, j, k)) > (ALRM(i, j, k) + ALTOL)) then
                IVF(i, j, k) = 1
              end if
            end do

!   ------ reset optional convection flag ------

            if (present(Conv)) then
              Conv(i, j, 1) = (IVF(i, j, 1) == 1)
              do k = 1, MXLEV1
                Conv(i, j, k + 1) = (IVF(i, j, k) == 1 .or. IVF(i, j, k + 1) == 1)
              end do
            end if

!   ------ Are all layers sufficiently stable ??? ------

            ISMVF(i, j) = 0
            do k = 1, MXLEV1
              ISMVF(i, j) = ISMVF(i, j) + IVF(i, j, k)
            end do
            if (ISMVF(i, j) == 0) cycle OUTER_LOOP

!-----------------------------------------------------------------------
1740        continue
!-----------------------------------------------------------------------

!---------Maximum iterations reached: Increase tolerance (ALTOL)--------
            ALTOL = 2.0*ALTOL
!del  WRITE (*,9902) I,ALTOL
            call error_mesg('moist_conv', 'Tolerence (ALTOL) doubled', NOTE)
            if (ALTOL <= TOLmax) GO TO 1450

!     WRITE (*,9903)
!     WRITE (*,9904) (k,Temp(i,j,k),Qmix(i,j,k),Qsat(i,j,k),  &
!                       Qdif(i,j,k),ALRM(i,j,k),k=1,MXLEV1)
!     WRITE (*,9904) (k,Temp(i,j,k),Qmix(i,j,k),Qsat(i,j,k),  &
!                       Qdif(i,j,k)            ,k=MXLEV,MXLEV)

            call error_mesg('moist_conv', 'maximum iterations reached', WARNING)
!-----------------------------------------------------------------------
          end do OUTER_LOOP
        end do
!-----------------------------------------------------------------------
!---------------------- END OF i,j LOOP --------------------------------
!-----------------------------------------------------------------------

!----- compute adjustments to temp and spec hum ----

        Tdel(:, :, :) = Temp(:, :, :) - Tin(:, :, :)
        Qdel(:, :, :) = Qmix(:, :, :) - Qin(:, :, :)

!----- integrate precip -----

        Rain(:, :) = 0.0
        Snow(:, :) = 0.0
        do k = 1, KX

          Rain(:, :) = Rain(:, :) + (Phalf(:, :, k) - Phalf(:, :, k + 1))* &
                       Qdel(:, :, k)*grav_inv

        end do
        Rain(:, :) = max(Rain(:, :), 0.0)
        Snow(:, :) = max(Snow(:, :), 0.0)
!-----------------------------------------------------------------------
!-----------------   PRINT FORMATS   -----------------------------------

9902    format(' *** ALTOL DOUBLED IN CONVAD AT I=', &
               I5, ' ,ALTOL=', F10.4)
9903    format(/, ' *** DIVERGENCE IN MOIST CONVECTIVE ADJUSTMENT ', /, &
                4x, 'K', 14x, 'T', 14x, 'R', 13x, 'Qsat', 14x, 'Qdif', 12x, 'ALRM',/)
9904    format(I5, 5e15.7)
!-----------------------------------------------------------------------

!------- update input values and compute tendency -------

        Tin = Tin + Tdel; Qin = Qin + Qdel

        Tdel = Tdel*dtinv; Qdel = Qdel*dtinv
        Rain = Rain*dtinv; Snow = Snow*dtinv
!---------------------------------------------------------------------
!   define the effect of moist convective adjustment on the tracer
!   fields. code to do so does not currently exist.
!---------------------------------------------------------------------
        if (present(qtrmca)) then
          qtrmca = 0.
        end if

!------- diagnostics for dt/dt_ras -------
        if (id_tdt_conv > 0) then
          used = send_data(id_tdt_conv, Tdel, Time, is, js, 1, &
                           rmask=mask)
        end if
!------- diagnostics for dq/dt_ras -------
        if (id_qdt_conv > 0) then
          used = send_data(id_qdt_conv, Qdel, Time, is, js, 1, &
                           rmask=mask)
        end if
!------- diagnostics for precip_ras -------
        if (id_prec_conv > 0) then
          used = send_data(id_prec_conv, Rain + Snow, Time, is, js)
        end if
!------- diagnostics for snow_ras -------
        if (id_snow_conv > 0) then
          used = send_data(id_snow_conv, Snow, Time, is, js)
        end if

!------- diagnostics for water vapor path tendency ----------
        if (id_q_conv_col > 0) then
          tempdiag(:, :) = 0.
          do k = 1, kx
            tempdiag(:, :) = tempdiag(:, :) + Qdel(:, :, k)*pmass(:, :, k)
          end do
          used = send_data(id_q_conv_col, tempdiag, Time, is, js)
        end if

!------- diagnostics for dry static energy tendency ---------
        if (id_t_conv_col > 0) then
          tempdiag(:, :) = 0.
          do k = 1, kx
            tempdiag(:, :) = tempdiag(:, :) + Tdel(:, :, k)*cp_air*pmass(:, :, k)
          end do
          used = send_data(id_t_conv_col, tempdiag, Time, is, js)
        end if

        do tr = 1, num_mca_tracers
!------- diagnostics for dtracer/dt from RAS -------------
          if (id_tracer_conv(tr) > 0) then
            used = send_data(id_tracer_conv(tr), qtrmca(:, :, :, tr), Time, is, js, 1, &
                             rmask=mask)
          end if

!------- diagnostics for column tracer path tendency -----
          if (id_tracer_conv_col(tr) > 0) then
            tempdiag(:, :) = 0.
            do k = 1, kx
              tempdiag(:, :) = tempdiag(:, :) + qtrmca(:, :, k, tr)*pmass(:, :, k)
            end do
            used = send_data(id_tracer_conv_col(tr), tempdiag, Time, is, js)
          end if

        end do

        end subroutine moist_conv

!#######################################################################

!#######################################################################

        subroutine moist_conv_init(axes, Time, tracers_in_mca)

          integer, intent(in) :: axes(4)
          type(time_type), intent(in) :: Time
          logical, dimension(:), intent(in), optional :: tracers_in_mca

!-----------------------------------------------------------------------

          integer :: unit, io, ierr
          integer :: nn, tr
          character(len=128) :: diagname, diaglname, tendunits, name, units

!-----------------------------------------------------------------------

          read (input_nml_file, nml=moist_conv_nml, iostat=io)
          ierr = check_nml_error(io, 'moist_conv_nml')

!---------- output namelist --------------------------------------------

          if (mpp_pe() == mpp_root_pe()) then
            call write_version_number(version, tagname)
            write (stdlog(), nml=moist_conv_nml)
          end if

          id_tdt_conv = register_diag_field(mod_name, &
                                            'tdt_conv', axes(1:3), Time, &
                                            'Temperature tendency from moist conv adj', 'K/s', &
                                            missing_value=missing_value)

          id_qdt_conv = register_diag_field(mod_name, &
                                            'qdt_conv', axes(1:3), Time, &
                                            'Spec humidity tendency from moist conv adj', 'kg/kg/s', &
                                            missing_value=missing_value)

          id_prec_conv = register_diag_field(mod_name, &
                                             'prec_conv', axes(1:2), Time, &
                                             'Precipitation rate from moist conv adj', 'kg/m2/s')

          id_snow_conv = register_diag_field(mod_name, &
                                             'snow_conv', axes(1:2), Time, &
                                             'Frozen precip rate from moist conv adj', 'kg/m2/s')

          id_q_conv_col = register_diag_field(mod_name, &
                                              'q_conv_col', axes(1:2), Time, &
                                              'Water vapor path tendency from moist conv adj', 'kg/m2/s')

          id_t_conv_col = register_diag_field(mod_name, &
                                              't_conv_col', axes(1:2), Time, &
                                              'Column static energy tendency from moist conv adj', 'W/m2')

!---------------------------------------------------------------------
! --- Find the tracer indices
!---------------------------------------------------------------------
          call get_number_tracers(MODEL_ATMOS, num_tracers)
          if (num_tracers .gt. 0) then
          else
            call error_mesg('moist_conv_init', 'No atmospheric tracers found', FATAL)
          end if

!----------------------------------------------------------------------
!    determine how many tracers are to be transported by moist_conv_mod.
!----------------------------------------------------------------------
          num_mca_tracers = count(tracers_in_mca)
          if (num_mca_tracers > 0) then
            do_mca_tracer = .true.
          else
            do_mca_tracer = .false.
          end if

!---------------------------------------------------------------------
!    allocate the arrays to hold the diagnostics for the moist_conv
!    tracers.
!---------------------------------------------------------------------
          allocate (id_tracer_conv(num_mca_tracers)); id_tracer_conv = 0
          allocate (id_tracer_conv_col(num_mca_tracers)); id_tracer_conv_col = 0
          nn = 1
          do tr = 1, num_tracers
            if (tracers_in_mca(tr)) then
              call get_tracer_names(MODEL_ATMOS, tr, name=name, units=units)

!----------------------------------------------------------------------
!    for the column tendencies, the name for the diagnostic will be
!    the name of the tracer followed by 'dt_MCA'. the longname will be
!    the name of the tracer followed by ' tendency from MCA'. units are
!    the supplied units of the tracer divided by seconds.
!----------------------------------------------------------------------
              diagname = trim(name)//'dt_MCA'
              diaglname = trim(name)//' tendency from MCA'
              tendunits = trim(units)//'/s'
              id_tracer_conv(nn) = register_diag_field(mod_name, &
                                                       trim(diagname), axes(1:3), Time, &
                                                       trim(diaglname), trim(tendunits), &
                                                       missing_value=missing_value)

!----------------------------------------------------------------------
!    for the column integral  tendencies, the name for the diagnostic
!    will be the name of the tracer followed by 'dt_MCA_col'. the long-
!    name will be the name of the tracer followed by ' path tendency
!    from MCA'. units are the supplied units of the tracer multiplied
!    by kg/m2 divided by seconds.
!----------------------------------------------------------------------
              diagname = trim(name)//'dt_MCA_col'
              diaglname = trim(name)//' path tendency from MCA'
              tendunits = trim(units)//' kg/m2/s'
              id_tracer_conv_col(nn) = register_diag_field(mod_name, &
                                                           trim(diagname), axes(1:2), Time, &
                                                           trim(diaglname), trim(tendunits), &
                                                           missing_value=missing_value)
              nn = nn + 1
            end if
          end do

          module_is_initialized = .true.

!-----------------------------------------------------------------------

        end subroutine moist_conv_init

!#######################################################################
        subroutine moist_conv_end

          integer :: log_unit

          log_unit = stdlog()
          if (mpp_pe() == mpp_root_pe()) then
            write (log_unit, '(/,(a))') 'Exiting moist_conv.'
          end if

          module_is_initialized = .false.

        end subroutine moist_conv_end

!#######################################################################

        end module moist_conv_mod
