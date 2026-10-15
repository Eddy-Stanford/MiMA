!> An arbitrarily specified tracer for testing the convective transport of tracers.
!>
!> The tracer has no sources or sinks. Up to `ncopies_cnvct_trcr` copies (`cnvct_trcr`,
!> `cnvct_trcr_2`, ...) are used if they are in the `field_table`; each is initialized with a
!> profile that decreases exponentially from 1 at the lowest level to exp(-1) at the top,
!> unless `INPUT/tracer_<name>.res` exists.
!>
!> Namelist: `atmos_convection_tracer_nml`
!> ([namelist reference](https://eddy-stanford.github.io/MiMA/Parameters/#atmos_convection_tracer_nml)).
!>
!> Original authors: Richard Hemler.
module atmos_convection_tracer_mod

!-----------------------------------------------------------------------

  use fms_mod, only: &
    write_version_number, &
    error_mesg, &
    FATAL, WARNING, NOTE, &
    mpp_pe, mpp_root_pe, stdlog
  use fms2_io_mod, only: file_exists
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

  public atmos_cnvct_tracer_sourcesink, &
    atmos_convection_tracer_init, &
    atmos_convection_tracer_end

!-----------------------------------------------------------------------
!----------- namelist -------------------

  integer  :: ncopies_cnvct_trcr = 9  !! number of copies of the convection tracer looked for in the field table;
                                     !! at most 9

  namelist /atmos_convection_tracer_nml/ &
    ncopies_cnvct_trcr

!-----------------------------------------------------------------------

!--- Arrays to help calculate tracer sources/sinks ---

  character(len=6), parameter :: module_name = 'tracer'

  logical :: module_is_initialized = .false.

!---- version number -----
  character(len=128) :: version = '$Id: atmos_convection_tracer.f90,v 11.0 2004/09/28 19:26:35 fms Exp $'
  character(len=128) :: tagname = '$Name: lima $'
!-----------------------------------------------------------------------

contains

!#######################################################################
  !> Returns the tendency of the convection tracer due to its sources and sinks, which is
  !> zero: the tracer is assumed to have no source or sink.
  subroutine atmos_cnvct_tracer_sourcesink(lon, lat, land, pwt, &
                                           convtr, convtr_dt, &
                                           Time, is, ie, js, je, &
                                           kbot)

!-----------------------------------------------------------------------
    real, intent(in), dimension(:, :)   :: lon, lat
    !! longitude and latitude of the centres of the grid cells [rad]
    real, intent(in), dimension(:, :)   :: land  !! land fraction
    real, intent(in), dimension(:, :, :) :: pwt, convtr
    !! `pwt`: pressure weight dp/grav [kg/m2]; `convtr`: convection tracer mixing ratio
    real, intent(out), dimension(:, :, :) :: convtr_dt  !! tendency of the convection tracer mixing ratio [1/s]
    type(time_type), intent(in) :: Time  !! model time
    integer, intent(in)       :: is, ie, js, je  !! local domain boundaries
    integer, intent(in), dimension(:, :), optional :: kbot  !! index of the lowest model level above the surface
!-----------------------------------------------------------------------
    real, dimension(size(convtr, 1), size(convtr, 2), size(convtr, 3)) :: &
      source, sink
!-----------------------------------------------------------------------

!------  define source and sink of convection_tracer -------
!
!   it is currently assumed that the convection tracer has no source
!   or sink

    source = 0.
    sink = 0.

!------- tendency ------------------

    convtr_dt = source + sink

!-----------------------------------------------------------------------

  end subroutine atmos_cnvct_tracer_sourcesink

!#######################################################################

  !> Initializes the convection tracer module: finds the convection tracers in the field
  !> table and sets their initial profile if there is no restart file for them.
  subroutine atmos_convection_tracer_init(r, phalf, axes, Time, &
                                          nconvect, mask)

!-----------------------------------------------------------------------
    real, intent(inout), dimension(:, :, :, :) :: r  !! tracer fields (nlon, nlat, nlev, ntrace)
    real, intent(in), dimension(:, :, :)   :: phalf  !! pressure at half levels [Pa]
    type(time_type), intent(in)                        :: Time  !! model time
    integer, intent(in)                        :: axes(4)  !! diagnostic axes (lon, lat, pfull, phalf)
    integer, dimension(:), pointer                         :: nconvect
    !! allocated here: tracer indices of the `ncopies_cnvct_trcr` copies (-1 if not in the field table)
    real, intent(in), dimension(:, :, :), optional        :: mask
    !! 1. above the ground, 0. below (nlon, nlat, nlev)

    logical :: flag
    integer :: n
    character(len=64) ::  search_name(10)
    character(len=4) ::  chname
    integer :: nn
!
!-----------------------------------------------------------------------
!
    real, dimension(size(r, 1), size(r, 2), size(r, 3)) :: xgcm, pfull
    integer log_unit, unit, io, index, ntr, nt
    character(len=16) ::  fld

    real :: xba = 1.0
    integer :: nlev, k
    character(len=64) :: filename

    nlev = size(r, 3)

!---------------------------------------------------------------------
    if (module_is_initialized) return

!---- write namelist ------------------

    call write_version_number(version, tagname)
    if (mpp_pe() == mpp_root_pe()) &
      write (stdlog(), nml=atmos_convection_tracer_nml)

!----- set initial value of convection tracer ------------

    if (ncopies_cnvct_trcr > 9) then
      call error_mesg('atmos_convection_tracer_mod', &
                      'currently no more than 9 copies of the convection tracer '// &
                      'are allowed', FATAL)
    end if
    allocate (nconvect(ncopies_cnvct_trcr))
    nconvect = -1

    do nn = 1, ncopies_cnvct_trcr
      write (chname, '(i1)') nn
      if (nn > 1) then
        search_name(nn) = 'cnvct_trcr_'//trim(chname)
      else
        search_name(nn) = 'cnvct_trcr'
      end if

      n = get_tracer_index(MODEL_ATMOS, search_name(nn))
      if (n > 0) then
        nconvect(nn) = n
        if (nconvect(nn) > 0 .and. mpp_pe() == mpp_root_pe()) write (*, 30) trim(search_name(nn)), nconvect(nn)
        if (nconvect(nn) > 0 .and. mpp_pe() == mpp_root_pe()) write (stdlog(), 30) trim(search_name(nn)), nconvect(nn)
      end if

    end do

30  format(A, ' was initialized as tracer number ', i2)
!

!---------------------------------------------------------------------
!    if a convection_tracer.res file exists, it will have been prev-
!    iously processed. there is no need to do anything here.
!---------------------------------------------------------------------
    do nn = 1, ncopies_cnvct_trcr
      if (nconvect(nn) > 0) then
        filename = 'INPUT/tracer_'//trim(search_name(nn))//'.res'
        if (file_exists(filename)) then

!--------------------------------------------------------------------
!    if a .res file does not exist, initialize the convection_tracer.
!--------------------------------------------------------------------
        else
          do k = 1, nlev
            pfull(:, :, k) = 0.5*(phalf(:, :, k) + phalf(:, :, k + 1))
          end do
          do k = 1, nlev
            xgcm(:, :, k) = xba* &
                            exp((pfull(:, :, k) - pfull(:, :, 1))/ &
                                (pfull(:, :, 1) - pfull(:, :, nlev)))
          end do
          do k = 1, nlev
            r(:, :, nlev + 1 - k, nconvect(nn)) = xgcm(:, :, k)
          end do
        end if  ! (file_exist)
      end if
    end do

    module_is_initialized = .true.

!-----------------------------------------------------------------------

  end subroutine atmos_convection_tracer_init

!#######################################################################

  !> Terminates the convection tracer module.
  subroutine atmos_convection_tracer_end

    module_is_initialized = .false.

  end subroutine atmos_convection_tracer_end

end module atmos_convection_tracer_mod

