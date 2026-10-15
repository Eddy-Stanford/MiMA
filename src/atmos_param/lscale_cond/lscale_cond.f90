
!> Large-scale condensation.
!>
!> Where the specific humidity exceeds `hc` times the saturation value, the temperature
!> and humidity are adjusted to saturation, releasing latent heat; the condensate falls
!> out as precipitation. With `do_evap`, falling precipitation re-evaporates in
!> sub-saturated layers below. All precipitation is returned as rain. Used with
!> `do_lsc = .true.` in `moist_processes_nml`.
!>
!> Namelist: `lscale_cond_nml`
!> ([namelist reference](https://eddy-stanford.github.io/MiMA/Parameters/#lscale_cond_nml)).
module lscale_cond_mod

!-----------------------------------------------------------------------
  use fms_mod, only: error_mesg, input_nml_file, &
                     check_nml_error, mpp_pe, mpp_root_pe, FATAL, &
                     write_version_number, stdlog
  use sat_vapor_pres_mod, only: escomp, descomp
  use constants_mod, only: HLv, HLs, Cp_Air, Grav, rdgas, rvgas

  implicit none
  private
!-----------------------------------------------------------------------
!  ---- public interfaces ----

  public lscale_cond, lscale_cond_init, lscale_cond_end

!-----------------------------------------------------------------------
!   ---- version number ----

  character(len=128) :: version = '$Id: lscale_cond.f90,v 10.0.6.1 2005/05/13 18:16:37 pjp Exp $'
  character(len=128) :: tagname = '$Name:  $'
  logical            :: module_is_initialized = .false.

!-----------------------------------------------------------------------
!   ---- local/private data ----

  real, parameter :: d622 = rdgas/rvgas
  real, parameter :: d378 = 1.-d622

!-----------------------------------------------------------------------
!   --- namelist ----

  real    :: hc = 1.00  !! relative humidity at which condensation occurs (0 <= `hc` <= 1)
  logical :: do_evap = .true.  !! re-evaporate falling precipitation in sub-saturated layers below

  namelist /lscale_cond_nml/ hc, do_evap

contains

!#######################################################################

  !> Computes the large-scale condensation: the adjustments of temperature and specific
  !> humidity and the resulting precipitation.
  subroutine lscale_cond(tin, qin, pfull, phalf, coldT, &
                         rain, snow, tdel, qdel, mask, conv)

!--------------------- interface arguments -----------------------------

    real, intent(in), dimension(:, :, :) :: tin, qin, pfull, phalf
    !! `tin`: temperature at full levels [K]; `qin`: specific humidity at full levels [kg/kg];
    !! `pfull`, `phalf`: pressure at full and half levels [Pa]
    logical, intent(in), dimension(:, :):: coldT  !! whether precipitation should be snow (not used)
    real, intent(out), dimension(:, :)   :: rain, snow
    !! liquid and frozen precipitation [kg/m2]; `snow` is always zero
    real, intent(out), dimension(:, :, :) :: tdel, qdel
    !! changes of temperature [K] and specific humidity [kg/kg] at full levels
    real, intent(in), dimension(:, :, :), optional :: mask  !! mask (0 or 1); no adjustment where it is 0
    logical, intent(in), dimension(:, :, :), optional :: conv
    !! no large-scale adjustment where true (e.g. where convection occurred)
!-----------------------------------------------------------------------
!---------------------- local data -------------------------------------

    logical, dimension(size(tin, 1), size(tin, 2), size(tin, 3)) :: do_adjust
    real, dimension(size(tin, 1), size(tin, 2), size(tin, 3)) :: &
      esat, qsat, desat, dqsat, pmes, pmass
    real, dimension(size(tin, 1), size(tin, 2))             :: hlcp, precip
    integer :: k, kx, i, j
!-----------------------------------------------------------------------
!     computation of precipitation by condensation processes
!-----------------------------------------------------------------------

    if (.not. module_is_initialized) call error_mesg('lscale_cond', &
                                                     'lscale_cond_init has not been called.', FATAL)

    kx = size(tin, 3)

!----- compute proper latent heat --------------------------------------
    hlcp = HLv/Cp_Air

!----- saturation vapor pressure (esat) & specific humidity (qsat) -----

    call escomp(tin, esat)
    call descomp(tin, desat)

    esat(:, :, :) = esat(:, :, :)*hc

    do k = 1, kx
    do j = 1, size(tin, 2)
    do i = 1, size(tin, 1)
      if (pfull(i, j, k) > d378*esat(i, j, k)) then
        pmes(i, j, k) = 1.0/pfull(i, j, k)
        qsat(i, j, k) = d622*esat(i, j, k)*pmes(i, j, k)
        qsat(i, j, k) = max(0.0, qsat(i, j, k))
        dqsat(i, j, k) = d622*pfull(i, j, k)*desat(i, j, k)*pmes(i, j, k)*pmes(i, j, k)
      else
        pmes(i, j, k) = 0.0
        qsat(i, j, k) = 0.0
        dqsat(i, j, k) = 0.0
      end if
    end do
    end do
    end do

!--------- do adjustment where greater than saturated value ------------

    if (present(conv)) then
!     do_adjust(:,:,:)=(.not.conv(:,:,:) .and. qin(:,:,:) > qsat(:,:,:))
      do_adjust(:, :, :) = (.not. conv(:, :, :) .and. &
                            (qin(:, :, :) - qsat(:, :, :))*qsat(:, :, :) > 0.0)
    else
!     do_adjust(:,:,:)=(qin(:,:,:) > qsat(:,:,:))
      do_adjust(:, :, :) = ((qin(:, :, :) - qsat(:, :, :))*qsat(:, :, :) > 0.0)
    end if

    if (present(mask)) then
      do_adjust(:, :, :) = do_adjust(:, :, :) .and. (mask(:, :, :) > 0.5)
    end if

!----------- compute adjustments to temp and spec humidity -------------
    do k = 1, kx
      where (do_adjust(:, :, k))
        qdel(:, :, k) = (qsat(:, :, k) - qin(:, :, k))/(1.0 + hlcp(:, :)*dqsat(:, :, k))
        tdel(:, :, k) = -hlcp(:, :)*qdel(:, :, k)
      elsewhere
        qdel(:, :, k) = 0.0
        tdel(:, :, k) = 0.0
      end where
    end do
!------------ pressure mass of each layer ------------------------------

    do k = 1, kx
      pmass(:, :, k) = (phalf(:, :, k + 1) - phalf(:, :, k))/Grav
    end do

!------------ re-evaporation of precipitation in dry layer below -------

    if (do_evap) then
      if (present(mask)) then
        call precip_evap(pmass, tin, qin, qsat, dqsat, hlcp, tdel, qdel, mask)
      else
        call precip_evap(pmass, tin, qin, qsat, dqsat, hlcp, tdel, qdel)
      end if
    end if

!------------ integrate precip -----------------------------------------

    precip(:, :) = 0.0
    do k = 1, kx
      precip(:, :) = precip(:, :) - pmass(:, :, k)*qdel(:, :, k)
    end do
    precip(:, :) = max(precip(:, :), 0.0)

    !assign precip to snow or rain
    rain = precip
    snow = 0.

!-----------------------------------------------------------------------

  end subroutine lscale_cond

!#######################################################################

  !> Re-evaporates falling precipitation in sub-saturated layers below.
  subroutine precip_evap(pmass, tin, qin, qsat, dqsat, hlcp, &
                         tdel, qdel, mask)

!-----------------------------------------------------------------------
    real, intent(in), dimension(:, :, :) :: pmass, tin, qin, qsat, dqsat
    real, intent(in), dimension(:, :)   :: hlcp
    real, intent(inout), dimension(:, :, :) :: tdel, qdel
    real, intent(in), dimension(:, :, :), optional :: mask
!-----------------------------------------------------------------------
    real, dimension(size(tin, 1), size(tin, 2)) :: exq, def

    integer k
!-----------------------------------------------------------------------
    exq(:, :) = 0.0

    do k = 1, size(tin, 3)

      where (qdel(:, :, k) < 0.0) exq(:, :) = exq(:, :) - &
                                              qdel(:, :, k)*pmass(:, :, k)

      if (present(mask)) exq(:, :) = exq(:, :)*mask(:, :, k)

!  ---- evaporate precip where needed ------

      where ((qdel(:, :, k) >= 0.0) .and. (exq(:, :) > 0.0))
        exq(:, :) = exq(:, :)/pmass(:, :, k)
        def(:, :) = (qsat(:, :, k) - qin(:, :, k))/(1.+hlcp(:, :)*dqsat(:, :, k))
        def(:, :) = min(max(def(:, :), 0.0), exq(:, :))
        qdel(:, :, k) = qdel(:, :, k) + def(:, :)
        tdel(:, :, k) = tdel(:, :, k) - def(:, :)*hlcp(:, :)
        exq(:, :) = (exq(:, :) - def(:, :))*pmass(:, :, k)
      end where

    end do

!-----------------------------------------------------------------------

  end subroutine precip_evap

!#######################################################################

  !> Initializes the module: reads `lscale_cond_nml` and writes it to the log file.
  subroutine lscale_cond_init()


    integer unit, io, ierr

!----------- read namelist ---------------------------------------------

    read (input_nml_file, nml=lscale_cond_nml, iostat=io)
    ierr = check_nml_error(io, 'lscale_cond_nml')

!---------- output namelist --------------------------------------------

    if (mpp_pe() == mpp_root_pe()) then
      call write_version_number(version, tagname)
      write (stdlog(), nml=lscale_cond_nml)
    end if

    module_is_initialized = .true.

  end subroutine lscale_cond_init

!#######################################################################
  !> Marks the module as not initialized.
  subroutine lscale_cond_end

    module_is_initialized = .false.

!---------------------------------------------------------------------

  end subroutine lscale_cond_end

!#######################################################################

end module lscale_cond_mod
