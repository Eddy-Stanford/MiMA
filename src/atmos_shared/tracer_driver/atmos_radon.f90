!> An extremely simplified radon tracer.
!>
!> A very simple tracer which bears some characteristics of radon (Rn222): a surface source
!> over land and radioactive decay. Up to `ncopies_radon` copies (`radon`, `radon_2`, ...)
!> are used if they are in the `field_table`.
!>
!> Namelist: `atmos_radon_nml`
!> ([namelist reference](https://eddy-stanford.github.io/MiMA/Parameters/#atmos_radon_nml)).
!>
!> Original authors: William Cooke.
module atmos_radon_mod

!-----------------------------------------------------------------------

  use fms_mod, only: &
    write_version_number, &
    mpp_pe, &
    mpp_root_pe, &
    error_mesg, &
    FATAL, WARNING, NOTE, &
    stdlog
  use time_manager_mod, only: time_type
  use diag_manager_mod, only: send_data
  use tracer_manager_mod, only: get_tracer_index
  use field_manager_mod, only: MODEL_ATMOS
  use atmos_tracer_utilities_mod, only: wet_deposition, &
                                        dry_deposition

  implicit none
  private
!-----------------------------------------------------------------------
!----- interfaces -------

  public atmos_radon_sourcesink, atmos_radon_init, atmos_radon_end

!-----------------------------------------------------------------------
!----------- namelist -------------------
!-----------------------------------------------------------------------
  integer  :: ncopies_radon = 9  !! number of copies of the radon tracer (`radon`, `radon_2`, ...) looked for in the
                                !! field table; at most 9

  namelist /atmos_radon_nml/ &
    ncopies_radon

!--- Arrays to help calculate tracer sources/sinks ---

  character(len=6), parameter :: module_name = 'tracer'

  logical :: module_is_initialized = .false.

!---- version number -----
  character(len=128) :: version = '$Id: atmos_radon.f90,v 11.0 2004/09/28 19:26:41 fms Exp $'
  character(len=128) :: tagname = '$Name: lima $'
!-----------------------------------------------------------------------

contains

!#######################################################################
  !> Computes the tendency of radon due to its sources and sinks.
  !>
  !> This is a very rudimentary implementation of radon. The Rn222 flux is assumed to be
  !> 3.69e-21 kg/m2/s over land between 60S and 60N, half of that between 60N and 70N
  !> (without `kbot`: except between 300E and 336E), and zero elsewhere; it is put into the
  !> lowest model level. The mixing ratio is scaled by 1e21. Rn222 has a half-life of
  !> 3.83 days, which corresponds to an e-folding time of 5.52 days.
  subroutine atmos_radon_sourcesink(lon, lat, land, pwt, radon, radon_dt, &
                                    Time, kbot)

!-----------------------------------------------------------------------
    real, intent(in), dimension(:, :)   :: lon, lat
    !! longitude and latitude of the centres of the grid cells [rad]
    real, intent(in), dimension(:, :)   :: land  !! land fraction (land where > 0.5)
    real, intent(in), dimension(:, :, :) :: pwt, radon
    !! `pwt`: pressure weight dp/grav [kg/m2]; `radon`: radon mixing ratio
    real, intent(out), dimension(:, :, :) :: radon_dt  !! tendency of the radon mixing ratio [1/s]
    type(time_type), intent(in) :: Time  !! model time
    integer, intent(in), dimension(:, :), optional :: kbot  !! index of the lowest model level above the surface
!-----------------------------------------------------------------------
    real, dimension(size(radon, 1), size(radon, 2), size(radon, 3)) :: &
      source, sink
    logical, dimension(size(radon, 1), size(radon, 2)) ::  maskeq, masknh
    real radon_flux, dtr, deg60, deg70, deg300, deg336
    integer i, j, kb, id, jd, kd, lat1
!-----------------------------------------------------------------------

    id = size(radon, 1); jd = size(radon, 2); kd = size(radon, 3)

    dtr = acos(0.0)/90.
    deg60 = 60.*dtr; deg70 = 70.*dtr; deg300 = 300.*dtr; deg336 = 336.*dtr

!----------- compute radon source ------------
!
!  rn222 flux is 3.69e-21 kg/m*m/sec over land for latitudes lt 60n
!   between 60n and 70n the source  = source * .5
!
!  molecular wt. of air is 28.9644 gm/mole
!  molecular wt. of radon is 222 gm/mole
!  scaling facter to get reasonable mixing ratio is 1.e+21
!
!  source = 3.69e-21 * g * 28.9644 * 1.e+21/(pwt * 222.) or
!
!  source = g * .4814353 / pwt
!
!  must initialize all rn to .001
!

    radon_flux = 3.69e-21*28.9644*1.e+21/222.
    source = 0.0
    maskeq = (land > 0.5) .and. lat > -deg60 .and. lat < deg60
    masknh = (land > 0.5) .and. lat >= deg60 .and. lat < deg70

    if (present(kbot)) then
      do j = 1, jd
      do i = 1, id
        kb = kbot(i, j)
        if (maskeq(i, j)) source(i, j, kb) = radon_flux/pwt(i, j, kb)
        if (masknh(i, j)) source(i, j, kb) = 0.5*radon_flux/pwt(i, j, kb)
      end do
      end do
    else
      where (maskeq) source(:, :, kd) = radon_flux/pwt(:, :, kd)
      where (masknh) source(:, :, kd) = 0.5*radon_flux/pwt(:, :, kd)
      where (masknh .and. lon > deg300 .and. lon < deg336) &
        source(:, :, kd) = 0.0
    end if

!------- compute radon sink --------------
!
!  rn222 has a half-life time of 3.83days
!   (corresponds to an e-folding time of 5.52 days)
!
!  sink = 1./(86400.*5.52) = 2.09675e-6
!

    where (radon(:, :, :) >= 0.0)
      sink(:, :, :) = -2.09675e-6*radon(:, :, :)
    elsewhere
      sink(:, :, :) = 0.0
    end where

!------- tendency ------------------

    radon_dt = source + sink

!-----------------------------------------------------------------------

  end subroutine atmos_radon_sourcesink

!#######################################################################

  !> Initializes the radon module: finds the radon tracers in the field table.
  subroutine atmos_radon_init(r, axes, Time, nradon, mask)

!-----------------------------------------------------------------------
    real, intent(inout), dimension(:, :, :, :) :: r  !! tracer fields (nlon, nlat, nlev, ntrace)
    type(time_type), intent(in)                        :: Time  !! model time
    integer, intent(in)                        :: axes(4)  !! diagnostic axes (lon, lat, pfull, phalf)
    integer, dimension(:), pointer                         :: nradon
    !! allocated here: tracer indices of the `ncopies_radon` radon copies (-1 if not in the field table)
    real, intent(in), dimension(:, :, :), optional        :: mask
    !! 1. above the ground, 0. below (nlon, nlat, nlev)

    logical :: flag
    integer :: n
!
!-----------------------------------------------------------------------
!
    integer log_unit, unit, io, index, ntr, nt
    character(len=16) ::  fld
    character(len=64) ::  search_name
    character(len=4) ::  chname
    integer :: nn

    if (module_is_initialized) return

!---- write namelist ------------------

    call write_version_number(version, tagname)
    if (mpp_pe() == mpp_root_pe()) &
      write (stdlog(), nml=atmos_radon_nml)

    if (ncopies_radon > 9) then
      call error_mesg('atmos_radonm_mod', &
                      'currently no more than 9 copies of the radon tracer '// &
                      'are allowed', FATAL)
    end if
    allocate (nradon(ncopies_radon))
    nradon = -1

    do nn = 1, ncopies_radon
      write (chname, '(i1)') nn
      if (nn > 1) then
        search_name = 'radon_'//trim(chname)
      else
        search_name = 'radon'
      end if
!----- set initial value of radon ------------

      n = get_tracer_index(MODEL_ATMOS, search_name)
      if (n > 0) then
        nradon(nn) = n
        if (nradon(nn) > 0 .and. mpp_pe() == mpp_root_pe()) write (*, 30) trim(search_name), nradon(nn)
        if (nradon(nn) > 0 .and. mpp_pe() == mpp_root_pe()) write (stdlog(), 30) trim(search_name), nradon(nn)
      end if

    end do

30  format(A, ' was initialized as tracer number ', i2)
!

    module_is_initialized = .true.

!-----------------------------------------------------------------------

  end subroutine atmos_radon_init

!#######################################################################

  !> Terminates the radon module.
  subroutine atmos_radon_end

    module_is_initialized = .false.

  end subroutine atmos_radon_end

end module atmos_radon_mod

