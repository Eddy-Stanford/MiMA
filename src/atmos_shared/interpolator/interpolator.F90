!> Interpolation of climatology files (e.g. ozone) to the model grid, time and levels.
!>
!> `interpolator_init` opens the netCDF file `INPUT/file_name` and sets up an
!> `interpolate_type` for some or all of its fields; `interpolator` then returns a field
!> interpolated conservatively to the model grid, to the model time and, for 3-D and 4-D
!> output, to the model pressure levels. The file's latitudes and longitudes may be in
!> degrees or radians, and its levels pressure (`pfull`, `phalf`, in Pa, hPa or mb) or
!> sigma (`sigma_full`, `sigma_half`). The time units are "days since" or "months since"
!> a date; a year of 0 marks a climatology for all years, and julian and noleap file
!> calendars are converted to the model calendar.
!>
!> A file with one record is used at all times, files with up to 12 records are
!> interpolated cyclically over the year, and longer files linearly in time. Files of
!> monthly data with gaps between years (e.g. only years ending in 0) are interpolated
!> between months and between the years with data. Files with 1 to 4 records are read
!> at initialization, others record by record as needed (or all at initialization with
!> `read_all_on_init`). Fields in column-integral units (kg/m2) are converted to mixing
!> ratios for the vertical interpolation and back. `init_clim_diag` registers diagnostics
!> of the column integrals.
!>
!> Namelist: `interpolator_nml`
!> ([namelist reference](https://eddy-stanford.github.io/MiMA/Parameters/#interpolator_nml)).
!>
!> Original authors: William Cooke.
module mima_interpolator_mod

  use mpp_mod, only: mpp_error, &
                     FATAL, &
                     mpp_pe, &
                     mpp_init, &
                     mpp_exit, &
                     mpp_npes, &
                     WARNING, &
                     NOTE
  use netcdf, only: nf90_open, nf90_close, nf90_inquire, &
                    nf90_inquire_dimension, &
                    nf90_inquire_variable, nf90_inq_varid, &
                    nf90_inq_dimid, nf90_inq_attname, &
                    nf90_inquire_attribute, nf90_get_att, &
                    nf90_get_var, nf90_strerror, &
                    NF90_NOWRITE, NF90_NOERR, NF90_CHAR, &
                    NF90_MAX_VAR_DIMS
  use mpp_domains_mod, only: mpp_domains_init, &
                             mpp_update_domains, &
                             mpp_define_domains, &
                             mpp_global_field, &
                             domain2d, &
                             mpp_define_layout, &
                             mpp_get_compute_domain
  use diag_manager_mod, only: diag_manager_init, get_base_time, &
                              register_diag_field, send_data, &
                              diag_axis_init
  use fms_mod, only: lowercase, write_version_number, &
                     fms_init, &
                     mpp_root_pe, stdlog
  use fms2_io_mod, only: file_exists
  use horiz_interp_mod, only: horiz_interp_type, &
                              horiz_interp_new, &
                              horiz_interp, &
                              horiz_interp_del
  use time_manager_mod, only: time_type, &
                              set_time, &
                              set_date, &
                              get_date, &
                              get_calendar_type, &
                              JULIAN, NOLEAP, &
                              THIRTY_DAY_MONTHS, & !mj
                              get_date_julian, set_date_no_leap, &
                              set_date_julian, get_date_no_leap, &
                              print_date, &
                              operator(+), &
                              operator(-), &
                              operator(*), &
                              operator(>), &
                              operator(<), &
                              decrement_time
  use time_interp_mod, only: time_interp, YEAR
  use constants_mod, only: grav, PI

  implicit none
  private

  public interpolator_init, &
    interpolator, &
    interpolator_end, &
    init_clim_diag, &
    query_interpolator

  !> Returns a climatology field interpolated to the model grid and time and, for 3-D
  !> and 4-D output, to the model pressure levels.
  !>
  !> `call interpolator(clim_type, Time, phalf, interp_data, field_name, is, js,
  !> clim_units)`; the 2-D version has no `phalf`. The 2-D and 3-D versions return the
  !> field `field_name`; the 4-D version returns all the fields of `clim_type` (last
  !> dimension of `interp_data`), which must then all have the same `vert_interp` and
  !> `data_out_of_bounds`. Stops if `field_name` is not a field of `clim_type`.
  interface interpolator
    module procedure interpolator_4D
    module procedure interpolator_3D
    module procedure interpolator_2D
  end interface

  !> Interpolates columns from one set of pressure layers to another, as the
  !> pressure-weighted average of the overlapping layers.
  interface interp_weighted_scalar
    module procedure interp_weighted_scalar_1D
    module procedure interp_weighted_scalar_2D
  end interface interp_weighted_scalar
  character(len=128) :: version = &
                        '$Id: interpolator.F90,v 12.0 2005/04/14 15:51:43 fms Exp $'
  character(len=128) :: tagname = '$Name: lima $'
  logical            :: module_is_initialized = .false.
  logical            :: clim_diag_initialized = .false.

! Metadata of a file dimension and its coordinate variable.
  type axis_info_type
    character(len=128) :: name
    integer            :: len
    character(len=128) :: units = 'nounits'
    character(len=128) :: calendar = 'unspecified'
    integer            :: sense = 0     ! 'positive' attribute: 1 up, -1 down
    real, allocatable  :: data(:)          ! zero for the time axis
  end type axis_info_type

! Metadata of a file variable that is not a coordinate variable.
  type field_info_type
    character(len=128)   :: name = 'noname'
    character(len=128)   :: units = 'nounits'
    integer              :: varid = -1
    real                 :: scale = 1.0    ! scale_factor
    real                 :: add = 0.0    ! add_offset
    integer, allocatable :: count(:)       ! read shape, 1 for the time dimension
    integer              :: tdim = 0      ! position of the time dimension, if any
  end type field_info_type

  !> A climatology file and the fields read from it, set up by `interpolator_init` (all
  !> components are private).
  type, public  :: interpolate_type
    private
!Redundant data between fields
!All climatology data
    real, pointer            :: lat(:) => null()  !! latitudes of the climatology grid [rad]
    real, pointer            :: lon(:) => null()  !! longitudes of the climatology grid [rad]
    real, pointer            :: latb(:) => null()  !! latitude boundaries of the climatology grid [rad]
    real, pointer            :: lonb(:) => null()  !! longitude boundaries of the climatology grid [rad]
    real, pointer            :: levs(:) => null()  !! climatology full levels [Pa], or sigma
    real, pointer            :: halflevs(:) => null()  !! climatology half levels [Pa], or sigma
    type(horiz_interp_type)  :: interph  !! horizontal interpolation to the model grid
    type(time_type), pointer :: time_slice(:) => null() !! times of the records in the file
    integer                  :: ncid          !! netCDF id of the open climatology file
    character(len=64)        :: file_name     !! climatology file name
    integer                  :: TIME_FLAG     !! time interpolation: `LINEAR`, `SEASONAL` or `BILINEAR`
    integer                  :: level_type    !! `PRESSURE` or `SIGMA` levels
    integer                  :: is, ie, js, je
    !! compute domain of the climatology grid, for the `climo` diagnostics
    integer                  :: vertical_indices !! direction of the vertical axis in the file:
    !! `INCREASING_DOWNWARD` or `INCREASING_UPWARD`

!Field specific data  for nfields
    type(field_info_type), pointer :: field_type(:) => null()   !! netCDF information of each field
    character(len=64), pointer :: field_name(:) => null()   !! name of each field
    integer, pointer :: time_init(:, :) => null()  !! records held in `data` (second index: 2 or ntime)
    integer, pointer :: mr(:) => null()           !! conversion to mixing ratio: `NO_CONV` or `KG_M2`
    integer, pointer :: out_of_bounds(:) => null()!! `CONSTANT` or `ZERO` (see `interpolator_init`)
!++lwh
    integer, pointer :: vert_interp(:) => null()  !! type of vertical interpolation
!--lwh
    real, pointer :: data(:, :, :, :, :) => null() !! horizontally interpolated data
    !! (lon, lat, climatology level, record, field)

    real, pointer :: pmon_pyear(:, :, :, :) => null()
    !! for `BILINEAR`: data of the month before the model time, in the data year before it
    real, pointer :: pmon_nyear(:, :, :, :) => null()
    !! for `BILINEAR`: data of the month before the model time, in the data year after it
    real, pointer :: nmon_nyear(:, :, :, :) => null()
    !! for `BILINEAR`: data of the month after the model time, in the data year after it
    real, pointer :: nmon_pyear(:, :, :, :) => null()
    !! for `BILINEAR`: data of the month after the model time, in the data year before it
!integer                    :: indexm, indexp, climatology
    integer, dimension(:), pointer :: indexm => null()  !! for `BILINEAR`: month before the model time
    integer, dimension(:), pointer :: indexp => null()  !! for `BILINEAR`: month after the model time
    integer, dimension(:), pointer :: climatology => null()
    !! for `BILINEAR`: data year before the model time

    type(time_type), pointer :: clim_times(:, :) => null()  !! times of the records, by (month, year)
  end type interpolate_type

  integer :: ndim, nvar, ntime
  integer :: nlat, nlatb, nlon, nlonb, nlev, nlevh
  integer ::          len, ntime_in, num_fields
  type(axis_info_type), allocatable :: axes(:)
  type(field_info_type), allocatable :: varfields(:)

! pletzer real, allocatable :: time_in(:)
! sjs real, allocatable :: climdata(:,:,:), climdata2(:,:,:)

  character(len=32) :: name, units
  integer           :: sense

  integer, parameter :: max_diag_fields = 30

! flags to indicate direction of vertical axis in  data file
  integer, parameter :: INCREASING_DOWNWARD = 1, INCREASING_UPWARD = -1
!++lwh
! Flags to indicate whether the time interpolation should be linear or some other scheme for seasonal data.
  integer, parameter :: LINEAR = 1, SEASONAL = 2, BILINEAR = 3

! Flags to indicate where climatology pressure levels are pressure or sigma levels
  integer, parameter :: PRESSURE = 1, SIGMA = 2

! Flags to indicate whether the climatology units are mixing ratio (kg/kg) or column integral (kg/m2).
! Vertical interpolation scheme requires mixing ratio at this time.
  integer, parameter :: NO_CONV = 1, KG_M2 = 2

! Flags to indicate what to do when the model surface pressure exceeds the  climatology surface pressure level.
  integer, parameter, public :: CONSTANT = 1, ZERO = 2
  !! values of `data_out_of_bounds`: with `CONSTANT` the climatology's top and bottom layers
  !! are extended to cover the model column; with `ZERO` they are not (with
  !! `INTERP_WEIGHTED_P` the field is then zero beyond the climatology's levels)

! Flags to indicate the type of vertical interpolation
  integer, parameter, public :: INTERP_WEIGHTED_P = 10, INTERP_LINEAR_P = 20, INTERP_LOG_P = 30
  !! values of `vert_interp`: `INTERP_WEIGHTED_P`, pressure-weighted average of the
  !! overlapping climatology layers; `INTERP_LINEAR_P`, linear in pressure;
  !! `INTERP_LOG_P` is not implemented (rejected by `interpolator_init`)
!--lwh

  integer :: num_clim_diag = 0
  character(len=64) :: climo_diag_name(max_diag_fields)
  integer :: climo_diag_id(max_diag_fields), hinterp_id(max_diag_fields)
  real ::  missing_value = -1.e10
! sjs integer :: itaum, itaup

  logical :: read_all_on_init = .false.  !! read all time levels of a file at initialization
  integer :: verbose = 0  !! amount of diagnostic printout

  namelist /interpolator_nml/ &
    read_all_on_init, verbose

contains
!
!#######################################################################
!
  !> Opens a climatology file and sets up an interpolate type for its fields.
  !>
  !> Reads the axes, times and field information of `INPUT/file_name`, sets up the
  !> horizontal interpolation to the model grid and allocates the data; seasonal data (1 to
  !> 4 records), and all data with `read_all_on_init`, are read here. Stops if the file or
  !> a field in `data_names` does not exist.
  subroutine interpolator_init(clim_type, file_name, lonb_mod, latb_mod, &
                               data_names, data_out_of_bounds, &
                               vert_interp, clim_units)
    type(interpolate_type), intent(inout) :: clim_type
    !! set up with the file and field data needed by `interpolator`
    character(len=*), intent(in)            :: file_name  !! name of the climatology file in `INPUT/`
    real, intent(in)            :: lonb_mod(:), latb_mod(:)
    !! `lonb_mod`, `latb_mod`: longitudes and latitudes of the model grid-box boundaries [rad]
    character(len=*), intent(in), optional :: data_names(:)
    !! names of the fields in the file to read (default: all fields)
!++lwh
    integer, intent(in)            :: data_out_of_bounds(:)
    !! what to do where the model column extends beyond the climatology's levels, `CONSTANT`
    !! or `ZERO`: one value for all fields, or one per field
    integer, intent(in), optional  :: vert_interp(:)
    !! type of vertical interpolation, `INTERP_WEIGHTED_P` (default) or `INTERP_LINEAR_P`:
    !! one value for all fields, or one per field
!--lwh
    character(len=*), intent(out), optional :: clim_units(:)  !! units of the fields read

    integer                      :: ncid, log_unit
    character(len=64)            :: src_file, cart, namelev
!++lwh
    integer                      :: num_files
    real                         :: dlat, dlon
!--lwh
    type(horiz_interp_type)      :: interph
    type(time_type), allocatable :: time_slice(:)
    type(time_type)              :: base_time
    type(time_type)              :: last_time
    logical                      :: NAME_PRESENT
    real                         :: dtr, tpi
    integer                      :: fileday, filemon, fileyr, filehr, filemin, filesec, m, m1
    character(len=20)           :: fileunits
    real, dimension(:), allocatable  :: alpha
    integer   :: j, k, ii, i
    logical :: non_monthly
    character(len=24) :: file_calendar
    integer :: model_calendar
    integer :: yr, mo, dy, hr, mn, sc
    integer :: n
    type(time_type) :: Julian_time, Noleap_time
    real, allocatable :: time_in(:), time_values(:)

    if (.not. module_is_initialized) then
      call fms_init
      call diag_manager_init
    end if

    tpi = 2.0*PI ! 4.*acos(0.)
    dtr = tpi/360.

    num_fields = 0

!--------------------------------------------------------------------
! open source file containing fields to be interpolated
!--------------------------------------------------------------------
    src_file = 'INPUT/'//trim(file_name)

    if (file_exists(trim(src_file))) then
      call nc_check(nf90_open(trim(src_file), NF90_NOWRITE, ncid), src_file)
    else
!Climatology file doesn't exist, so exit
      call mpp_error(FATAL, 'Interpolator_init : Data file '//trim(src_file)//' does not exist')
    end if

!Read the axes, the fields (nvar of them) and the times in this file
    call get_file_info(ncid, src_file, axes, varfields, time_values, ntime)
    if (ntime == 0) call mpp_error(FATAL, 'interpolator_init : '//trim(src_file)//' has no time records')
    ndim = size(axes)
    nvar = size(varfields)
    clim_type%ncid = ncid
    clim_type%file_name = trim(file_name)

    num_fields = nvar
    if (present(data_names)) num_fields = size(data_names(:))

    nlon = 0 ! Number of longitudes (center-points) in the climatology.
    nlat = 0 ! Number of latitudes (center-points) in the climatology.
    nlev = 0 ! Number of levels (center-points) in the climatology.
    nlatb = 0 ! Number of longitudes (boundaries) in the climatology.
    nlonb = 0 ! Number of latitudes (boundaries) in the climatology.
    nlevh = 0 ! Number of levels (boundaries) in the climatology.

    clim_type%level_type = 0 ! Default value

!++lwh
! -------------------------------------------------------------------
! For 2-D fields, set a default value of nlev=nlevh=1
! -------------------------------------------------------------------
    nlev = 1
    nlevh = 1
!--lwh
    clim_type%vertical_indices = 0  ! initial value

    do i = 1, ndim
      name = axes(i)%name
      len = axes(i)%len
      units = axes(i)%units
      file_calendar = axes(i)%calendar
      sense = axes(i)%sense
      !mj if we want to use a previous output file as an input file
      ! (eg to force the model), the calendar might be "360_day" (CF) or
      ! the older "360", both of which mean "thirty_day_months"
      if (trim(file_calendar) .eq. '360_day' .or. trim(file_calendar) .eq. '360') &
        file_calendar = 'thirty_day_months'
      select case (name)
      case ('lat')
        nlat = len
        allocate (clim_type%lat(nlat))
        clim_type%lat = axes(i)%data
        select case (units(1:6))
        case ('degree')
          clim_type%lat = clim_type%lat*dtr
        case ('radian')
        case default
          call mpp_error(FATAL, "interpolator_init : Units for lat not recognised in file "//file_name)
        end select
      case ('lon')
        nlon = len
        allocate (clim_type%lon(nlon))
        clim_type%lon = axes(i)%data
        select case (units(1:6))
        case ('degree')
          clim_type%lon = clim_type%lon*dtr
        case ('radian')
        case default
          call mpp_error(FATAL, "interpolator_init : Units for lon not recognised in file "//file_name)
        end select
      case ('latb')
        nlatb = len
        allocate (clim_type%latb(nlatb))
        clim_type%latb = axes(i)%data
        select case (units(1:6))
        case ('degree')
          clim_type%latb = clim_type%latb*dtr
        case ('radian')
        case default
          call mpp_error(FATAL, "interpolator_init : Units for latb not recognised in file "//file_name)
        end select
      case ('lonb')
        nlonb = len
        allocate (clim_type%lonb(nlonb))
        clim_type%lonb = axes(i)%data
        select case (units(1:6))
        case ('degree')
          clim_type%lonb = clim_type%lonb*dtr
        case ('radian')
        case default
          call mpp_error(FATAL, "interpolator_init : Units for lonb not recognised in file "//file_name)
        end select
      case ('pfull')
        nlev = len
        allocate (clim_type%levs(nlev))
        clim_type%levs = axes(i)%data
        clim_type%level_type = PRESSURE
        ! Convert to Pa
        if (chomp(units) == "mb" .or. chomp(units) == "hPa") then
          clim_type%levs = clim_type%levs*100.
        end if
! define the direction of the vertical data axis
! switch index order if necessary so that indx 1 is at lowest pressure,
! index nlev at highest pressure.
        if (sense == 1) then
          clim_type%vertical_indices = INCREASING_UPWARD
          allocate (alpha(nlev))
          do n = 1, nlev
            alpha(n) = clim_type%levs(nlev - n + 1)
          end do
          do n = 1, nlev
            clim_type%levs(n) = alpha(n)
          end do
          deallocate (alpha)
        else
          clim_type%vertical_indices = INCREASING_DOWNWARD
        end if

      case ('phalf')
        nlevh = len
        allocate (clim_type%halflevs(nlevh))
        clim_type%halflevs = axes(i)%data
        clim_type%level_type = PRESSURE
        ! Convert to Pa
        if (chomp(units) == "mb" .or. chomp(units) == "hPa") then
          clim_type%halflevs = clim_type%halflevs*100.
        end if
! define the direction of the vertical data axis
! switch index order if necessary so that indx 1 is at lowest pressure,
! index nlev at highest pressure.
        if (sense == 1) then
          clim_type%vertical_indices = INCREASING_UPWARD
          allocate (alpha(nlevh))
          do n = 1, nlevh
            alpha(n) = clim_type%halflevs(nlevh - n + 1)
          end do
          do n = 1, nlevh
            clim_type%halflevs(n) = alpha(n)
          end do
          deallocate (alpha)
        else
          clim_type%vertical_indices = INCREASING_DOWNWARD
        end if
      case ('sigma_full')
        nlev = len
        allocate (clim_type%levs(nlev))
        clim_type%levs = axes(i)%data
        clim_type%level_type = SIGMA
      case ('sigma_half')
        nlevh = len
        allocate (clim_type%halflevs(nlevh))
        clim_type%halflevs = axes(i)%data
        clim_type%level_type = SIGMA

      case ('time')
        model_calendar = get_calendar_type()
        fileday = 0
        filemon = 0
        fileyr = 0
        filehr = 0
        filemin = 0
        filesec = 0
        select case (units(:3))
        case ('day')
          fileunits = units(12:) !Assuming "days since YYYY-MM-DD HH:MM:SS"
          read (fileunits(1:4), *) fileyr
          read (fileunits(6:7), *) filemon
          read (fileunits(9:10), *) fileday
          read (fileunits(12:13), *) filehr
          read (fileunits(15:16), *) filemin
          read (fileunits(18:19), *) filesec
        case ('mon')
          fileunits = units(14:) !Assuming "months since YYYY-MM-DD HH:MM:SS"
          read (fileunits(1:4), *) fileyr
          read (fileunits(6:7), *) filemon
          read (fileunits(9:10), *) fileday
          read (fileunits(12:13), *) filehr
          read (fileunits(15:16), *) filemin
          read (fileunits(18:19), *) filesec
        case default
          call mpp_error(FATAL, 'Interpolator_init : Time units not recognised in file '//file_name)
        end select

        if (fileyr /= 0) then

!----------------------------------------------------------------------
!    if file date has a non-zero year in the base time, determine that
!    base_time based on the netcdf info.
!----------------------------------------------------------------------
          if ((model_calendar == JULIAN .and. &
               trim(file_calendar) == 'julian') .or. &
              (model_calendar == NOLEAP .and. &
               trim(file_calendar) == 'noleap') .or. &
              (model_calendar == THIRTY_DAY_MONTHS .and. & !mj
               trim(file_calendar) == 'thirty_day_months')) then
            call mpp_error(NOTE, 'mima_interpolator_mod: Model and file&
                      & calendars are the same ( '//trim(file_calendar)//' ) for file '//   &
                      & trim(file_name)//'; no calendar conversion  &
                      &needed')
            base_time = set_date(fileyr, filemon, fileday, filehr, &
                                 filemin, filesec)
          else if ((model_calendar == JULIAN .and. &
                    trim(file_calendar) == 'noleap')) then
            call mpp_error(NOTE, 'mima_interpolator_mod: Using julian &
                              &model calendar and noleap file calendar&
                              & for file '//trim(file_name)//   &
                              &'; calendar conversion needed')
            base_time = set_date_no_leap(fileyr, filemon, fileday, &
                                         filehr, filemin, filesec)
          else if ((model_calendar == NOLEAP .and. &
                    trim(file_calendar) == 'julian')) then
            call mpp_error(NOTE, 'mima_interpolator_mod: Using noleap &
                              &model calendar and julian file calendar&
                              & for file '//trim(file_name)//  &
                              &'; calendar conversion needed')
            base_time = set_date_julian(fileyr, filemon, fileday, &
                                        filehr, filemin, filesec)
          else
            call mpp_error(FATAL, 'mima_interpolator_mod: Model and file&
                 & calendars ( '//trim(file_calendar)//' ) differ  &
                 &for file '//trim(file_name)//';  this calendar  &
                 &conversion not currently available')
          end if

        else

! if the year is specified as '0000', then the file is intended to
! apply to all years -- the time variables within the file refer to
! the displacement from the start of each year to the time of the
! associated data. Time interpolation is to be done with interface
! time_interp_list, with the optional argument modtime=YEAR. base_time
! is set to an arbitrary value here; it's only use will be as a
! timestamp for optionally generated diagnostics.

          base_time = get_base_time()
        end if

        ntime_in = 1
        if (ntime > 0) then
          allocate (time_in(ntime), clim_type%time_slice(ntime))
          allocate (clim_type%clim_times(12, (ntime + 11)/12))
          time_in = 0.0
          clim_type%time_slice = set_time(0, 0) + base_time
          clim_type%clim_times = set_time(0, 0) + base_time
          time_in = time_values
          ntime_in = ntime
! determine whether the data is a continuous set of monthly values or
! a series of annual cycles spread throughout the period of data
          non_monthly = .false.
          do n = 1, ntime - 1
!  Assume that the times in the data file correspond to days only.
            if (time_in(n + 1) > (time_in(n) + 32)) then
              non_monthly = .true.
              exit
            end if
          end do
          if (fileyr == 0) then
            call mpp_error(NOTE, 'mima_interpolator_mod :'// &
                           trim(file_name)//' is a year-independent climatology file')
          else
            call mpp_error(NOTE, 'mima_interpolator_mod :'// &
                           trim(file_name)//' is a timeseries file')
          end if

          do n = 1, ntime
!Assume that the times in the data file correspond to days only.

            if (fileyr == 0) then
! RSH NOTE:
! for this case, do not add base_time. time_slice will be sent to
! time_interp_list with the optional argument modtime=YEAR, so that
! the time that is needed in time_slice is the displacement into the
! year, not the displacement from a base_time.
              clim_type%time_slice(n) = set_time(0, int(time_in(n)))
            else

!--------------------------------------------------------------------
!    if fileyr /= 0, then define the times associated with each time-
!    slice. if calendar conversion between data file and model calendar
!    is needed, do it so that data from the file is associated with the
!    same calendar time in the model. here the time_slice needs to
!    include the base_time; values will be generated relative to the
!    "real" time.
!--------------------------------------------------------------------
              if ((model_calendar == JULIAN .and. &
                   trim(file_calendar) == 'julian') .or. &
                  (model_calendar == NOLEAP .and. &
                   trim(file_calendar) == 'noleap') .or. &
                  (model_calendar == THIRTY_DAY_MONTHS .and. & !mj
                   trim(file_calendar) == 'thirty_day_months')) then

!---------------------------------------------------------------------
!    no calendar conversion needed.
!---------------------------------------------------------------------
                clim_type%time_slice(n) = set_time(0, int(time_in(n))) + &
                                          base_time

!---------------------------------------------------------------------
!    convert file times from noleap to julian.
!---------------------------------------------------------------------
              else if ((model_calendar == JULIAN .and. &
                        trim(file_calendar) == 'noleap')) then
                Noleap_time = set_time(0, int(time_in(n))) + base_time
                call get_date_no_leap(Noleap_time, yr, mo, dy, hr, &
                                      mn, sc)
                clim_type%time_slice(n) = set_date_julian(yr, mo, dy, &
                                                          hr, mn, sc)
                if (n == 1) then
                  call print_date(clim_type%time_slice(1), &
                          str='for file '//trim(file_name)//', the &
                                &first time slice is mapped to :')
                end if
                if (n == ntime) then
                  call print_date(clim_type%time_slice(ntime), &
                           str='for file '//trim(file_name)//', the &
                                 &last time slice is mapped to:')
                end if

!---------------------------------------------------------------------
!    convert file times from julian to noleap.
!---------------------------------------------------------------------
              else if ((model_calendar == NOLEAP .and. &
                        trim(file_calendar) == 'julian')) then
                Julian_time = set_time(0, int(time_in(n))) + base_time
                call get_date_julian(Julian_time, yr, mo, dy, hr, mn, sc)
                clim_type%time_slice(n) = set_date_no_leap(yr, mo, dy, &
                                                           hr, mn, sc)
                if (n == 1) then
                  call print_date(clim_type%time_slice(1), &
                           str='for file '//trim(file_name)//', the &
                                 &first time slice is mapped to :')
                end if
                if (n == ntime) then
                  call print_date(clim_type%time_slice(ntime), &
                           str='for file '//trim(file_name)//', the &
                                 &last time slice is mapped to:')
                end if

!---------------------------------------------------------------------
!    any other calendar combinations would have caused a fatal error
!    above.
!---------------------------------------------------------------------
              end if
            end if

            m = (n - 1)/12 + 1; m1 = n - (m - 1)*12
            clim_type%clim_times(m1, m) = clim_type%time_slice(n)
          end do
          deallocate (time_in)
        end if
      end select ! case(name)
    end do

    if (.not. associated(clim_type%time_slice)) then
      if (ntime > 0) call mpp_error(FATAL, 'interpolator_init : the time dimension of '// &
                                    trim(file_name)//' must be called time')
! No time dimension: the file holds one time-independent record.
      ntime = 1
      allocate (clim_type%time_slice(1), clim_type%clim_times(1, 1))
      base_time = get_base_time()
      clim_type%time_slice = base_time
      clim_type%clim_times = base_time
    end if

! -------------------------------------------------------------------
! For 2-D fields, allocate levs and halflevs here
!  code is still needed for case when only halflevs are in data file.
! -------------------------------------------------------------------
    if (.not. associated(clim_type%levs)) then
      allocate (clim_type%levs(nlev))
      clim_type%levs = 0.0
    end if
    if (.not. associated(clim_type%halflevs)) then
      allocate (clim_type%halflevs(nlev + 1))
      clim_type%halflevs(1) = 0.0
      if (clim_type%level_type == PRESSURE) then
        clim_type%halflevs(nlev + 1) = 1013.25*100.0   ! MKS
      else if (clim_type%level_type == SIGMA) then
        clim_type%halflevs(nlev + 1) = 1.0
      end if
      do n = 2, nlev
        clim_type%halflevs(n) = 0.5*(clim_type%levs(n) + &
                                     clim_type%levs(n - 1))
      end do
    end if
    deallocate (axes, time_values)

! In the case where only the midpoints of the longitudes are defined we force the definition
! of the boundaries to be half-way between the midpoints.
    if (.not. associated(clim_type%lon) .and. .not. associated(clim_type%lonb)) &
      call mpp_error(FATAL, 'Interpolator_init : There appears to be no longitude axis in file '//file_name)

    if (.not. associated(clim_type%lonb)) then

      if (size(clim_type%lon(:)) /= 1) then
        allocate (clim_type%lonb(size(clim_type%lon(:)) + 1))
        dlon = (clim_type%lon(2) - clim_type%lon(1))/2.0
        clim_type%lonb(1) = clim_type%lon(1) - dlon
        clim_type%lonb(2:) = clim_type%lon(1:) + dlon
      else

! this is the case for zonal mean data, lon = 1, lonb not present
! in file.

        allocate (clim_type%lonb(2))
        clim_type%lonb(1) = -360.*dtr
        clim_type%lonb(2) = 360.0*dtr
        clim_type%lon(1) = 0.0
      end if
    end if

!clim_type%lonb=clim_type%lonb*dtr
! This assumes the lonb are in degrees in the NetCDF file!

    if (.not. associated(clim_type%lat) .and. .not. associated(clim_type%latb)) &
      call mpp_error(FATAL, 'Interpolator_init : There appears to be no latitude axis in file '//file_name)
! In the case where only the grid midpoints of the latitudes are defined we force the
! definition of the boundaries to be half-way between the midpoints.
    if (.not. associated(clim_type%latb)) then
      allocate (clim_type%latb(nlat + 1))
      dlat = (clim_type%lat(2) - clim_type%lat(1))*0.5
!  clim_type%latb(1) = min( 90., max(-90., clim_type%lat(1) - dlat) )
      clim_type%latb(1) = min(PI/2., max(-PI/2., clim_type%lat(1) - dlat))
      clim_type%latb(2:nlat) = (clim_type%lat(1:nlat - 1) + clim_type%lat(2:nlat))*0.5
      dlat = (clim_type%lat(nlat) - clim_type%lat(nlat - 1))*0.5
!  clim_type%latb(nlat+1) = min( 90., max(-90., clim_type%lat(nlat) + dlat) )
      clim_type%latb(nlat + 1) = min(PI/2., max(-PI/2., clim_type%lat(nlat) + dlat))
    end if
!clim_type%latb=clim_type%latb*dtr

!Assume that the horizontal interpolation within a file is the same for each variable.

    call horiz_interp_new(clim_type%interph, &
                          clim_type%lonb, clim_type%latb, &
                          lonb_mod, latb_mod, interp_method='conservative')

!--------------------------------------------------------------------
!  allocate the variable clim_type%data . This will be the climatology
!  data horizontally interpolated, so it will be on the model horizontal
!  grid, but it will still be on the climatology vertical grid.
!--------------------------------------------------------------------

    select case (ntime)
    case (13:)
! This may  be data that does not have a continous time-line
! i.e. IPCC data where decadal data is present but we wish to retain
! the seasonal nature of the data.
! RSH: the following test will not always work; instead use the
! RSH: non-monthly variable to test on.
!RSHlast_time = clim_type%time_slice(1) + ( ntime -1 ) * &
!RSH        ( clim_type%time_slice(2) - clim_type%time_slice(1) )

!RSHif ( last_time < clim_type%time_slice(ntime)) then

      if (non_monthly) then
! We have a broken time-line. e.g. We have monthly data but only for years ending in 0. 1960,1970 etc.
        if (mod(ntime, 12) /= 0) call mpp_error(FATAL, 'interpolator_init : '//trim(file_name)// &
                                                ' has gaps in time but is not a set of whole years of monthly data')
!   allocate(clim_type%data(size(lonb_mod(:))-1, size(latb_mod(:))-1, nlev, 2, num_fields))
        allocate (clim_type%pmon_pyear(size(lonb_mod(:)) - 1, size(latb_mod(:)) - 1, nlev, num_fields))
        allocate (clim_type%pmon_nyear(size(lonb_mod(:)) - 1, size(latb_mod(:)) - 1, nlev, num_fields))
        allocate (clim_type%nmon_nyear(size(lonb_mod(:)) - 1, size(latb_mod(:)) - 1, nlev, num_fields))
        allocate (clim_type%nmon_pyear(size(lonb_mod(:)) - 1, size(latb_mod(:)) - 1, nlev, num_fields))
        clim_type%pmon_pyear = 0.0
        clim_type%pmon_nyear = 0.0
        clim_type%nmon_nyear = 0.0
        clim_type%nmon_pyear = 0.0
        clim_type%TIME_FLAG = BILINEAR
      else
! We have a continuous time-line so treat as for 5-12 timelevels as below.
        if (.not. read_all_on_init) then
          allocate (clim_type%data(size(lonb_mod(:)) - 1, size(latb_mod(:)) - 1, nlev, 2, num_fields))
        else
          allocate (clim_type%data(size(lonb_mod(:)) - 1, size(latb_mod(:)) - 1, nlev, &
                                   ntime, num_fields))
        end if
        clim_type%data = 0.0
        clim_type%TIME_FLAG = LINEAR
      end if

    case (5:12)
! We have more than 4 timelevels
! Assume we have monthly or higher time resolution datasets (climatology or time series)
! So we only need to read 2 datasets and apply linear temporal interpolation.
      if (.not. read_all_on_init) then
        allocate (clim_type%data(size(lonb_mod(:)) - 1, size(latb_mod(:)) - 1, nlev, 2, num_fields))
      else
        allocate (clim_type%data(size(lonb_mod(:)) - 1, size(latb_mod(:)) - 1, nlev, &
                                 ntime, num_fields))
      end if
      clim_type%data = 0.0
      clim_type%TIME_FLAG = LINEAR
    case (1:4)
! Assume we have seasonal data and read in all the data.
! They are interpolated linearly in time, cyclically over the year.

      allocate (clim_type%data(size(lonb_mod(:)) - 1, size(latb_mod(:)) - 1, nlev, ntime, num_fields))
      clim_type%data = 0.0
      clim_type%TIME_FLAG = SEASONAL
! case (default)

    end select

!------------------------------------------------------------------
!    Allocate space for the single time level of the climatology on its
!    grid size.
!----------------------------------------------------------------------

    if (clim_type%TIME_FLAG .eq. LINEAR) then
      allocate (clim_type%time_init(num_fields, 2))
    else
      allocate (clim_type%time_init(num_fields, ntime))
    end if
    allocate (clim_type%indexm(num_fields), &
              clim_type%indexp(num_fields), &
              clim_type%climatology(num_fields))
    clim_type%time_init(:, :) = 0
    clim_type%indexm(:) = 0
    clim_type%indexp(:) = 0
    clim_type%climatology(:) = 0

    allocate (clim_type%field_name(num_fields))
    allocate (clim_type%field_type(num_fields))
    allocate (clim_type%mr(num_fields))
    allocate (clim_type%out_of_bounds(num_fields))
    clim_type%out_of_bounds(:) = 0
    allocate (clim_type%vert_interp(num_fields))
    clim_type%vert_interp(:) = 0
    if (present(data_names)) then

!++lwh
      if (size(data_out_of_bounds(:)) /= size(data_names(:)) .and. size(data_out_of_bounds(:)) /= 1) &
         call mpp_error(FATAL, 'interpolator_init : The size of the data_out_of_bounds array must be 1&
                               & or size(data_names)')
      if (present(vert_interp) .and. &
          size(vert_interp(:)) /= size(data_names(:)) .and. size(vert_interp(:)) /= 1) &
         call mpp_error(FATAL, 'interpolator_init : The size of the vert_interp array must be 1&
                               & or size(data_names)')
! Only read the fields named in data_names
      do j = 1, size(data_names(:))
        NAME_PRESENT = .false.
        do i = 1, nvar
          name = varfields(i)%name
          units = varfields(i)%units
          if (name == data_names(j)) then
            units = chomp(units)
            if (mpp_pe() == 0) write (*, *) 'Initializing src field : ', trim(name)
            clim_type%field_name(j) = name
            clim_type%field_type(j) = varfields(i)
            clim_type%mr(j) = check_climo_units(units)
            NAME_PRESENT = .true.
            if (present(clim_units)) clim_units(j) = units
            clim_type%out_of_bounds(j) = data_out_of_bounds(min(j, size(data_out_of_bounds(:))))
            if (clim_type%out_of_bounds(j) /= CONSTANT .and. &
                clim_type%out_of_bounds(j) /= ZERO) &
               call mpp_error(FATAL, "Interpolator_init: data_out_of_bounds must be&
                                    & set to ZERO or CONSTANT")
            if (present(vert_interp)) then
              clim_type%vert_interp(j) = vert_interp(min(j, size(vert_interp(:))))
              if (clim_type%vert_interp(j) /= INTERP_WEIGHTED_P .and. &
                  clim_type%vert_interp(j) /= INTERP_LINEAR_P) &
                 call mpp_error(FATAL, "Interpolator_init: vert_interp must be&
                                      & set to INTERP_WEIGHTED_P or INTERP_LINEAR_P")
            else
              clim_type%vert_interp(j) = INTERP_WEIGHTED_P
            end if
          end if
        end do
        if (.not. NAME_PRESENT) &
          call mpp_error(FATAL, 'interpolator_init : Check names of fields being passed. ' &
                         //trim(data_names(j))//' does not exist.')
      end do
    else

      if (size(data_out_of_bounds(:)) /= nvar .and. size(data_out_of_bounds(:)) /= 1) &
         call mpp_error(FATAL, 'interpolator_init : The size of the out of bounds array must be 1&
                              & or the number of fields in the climatology dataset')
      if (present(vert_interp)) then
        if (size(vert_interp(:)) /= nvar .and. size(vert_interp(:)) /= 1) &
        call mpp_error(FATAL, 'interpolator_init : The size of the vert_interp array must be 1&
                             & or the number of fields in the climatology dataset')
      end if

! Read all the fields within the climatology data file.
      do i = 1, nvar
        name = varfields(i)%name
        units = varfields(i)%units
        if (mpp_pe() == 0) write (*, *) 'Initializing src field : ', trim(name)
        clim_type%field_name(i) = lowercase(trim(name))
        clim_type%field_type(i) = varfields(i)
        clim_type%mr(i) = check_climo_units(units)
        if (present(clim_units)) clim_units(i) = units
        clim_type%out_of_bounds(i) = data_out_of_bounds(min(i, size(data_out_of_bounds(:))))
        if (clim_type%out_of_bounds(i) /= CONSTANT .and. &
            clim_type%out_of_bounds(i) /= ZERO) &
           call mpp_error(FATAL, "Interpolator_init: data_out_of_bounds must be&
                                & set to ZERO or CONSTANT")
        if (present(vert_interp)) then
          clim_type%vert_interp(i) = vert_interp(min(i, size(vert_interp(:))))
          if (clim_type%vert_interp(i) /= INTERP_WEIGHTED_P .and. &
              clim_type%vert_interp(i) /= INTERP_LINEAR_P) &
             call mpp_error(FATAL, "Interpolator_init: vert_interp must be&
                                  & set to INTERP_WEIGHTED_P or INTERP_LINEAR_P")
        else
          clim_type%vert_interp(i) = INTERP_WEIGHTED_P
        end if
      end do
!--lwh
    end if

    deallocate (varfields)

    if (clim_type%TIME_FLAG .eq. SEASONAL) then
! Read all the data at this point.
      do i = 1, num_fields
        do n = 1, ntime
          call read_data(clim_type, clim_type%field_type(i), &
                         clim_type%data(:, :, :, n, i), n, i, base_time)
        end do
      end do
    end if

    if (clim_type%TIME_FLAG .eq. LINEAR .and. read_all_on_init) then
! Read all the data at this point.
      do i = 1, num_fields
        do n = 1, ntime
          call read_data(clim_type, clim_type%field_type(i), &
                         clim_type%data(:, :, :, n, i), n, i, base_time)
        end do
      end do

      call nc_check(nf90_close(ncid), src_file)
    end if

    module_is_initialized = .true.

    call write_version_number(version, tagname)

  end subroutine interpolator_init
!
!#######################################################################
!
  !> Returns `KG_M2` if the units are those of a column integral (kg/m2), else `NO_CONV`.
  function check_climo_units(units)
! Function to check the units that the climatology data is using.
! This is needed to allow for conversion of datasets to mixing ratios which is what the
! vertical interpolation scheme requires
! The default is to assume no conversion is needed.
! If the units are those of a column burden (kg/m2) then conversion to mixing ratio is flagged.
!
    character(len=*), intent(in) :: units

    integer :: check_climo_units

    check_climo_units = NO_CONV
    select case (chomp(units))
    case ('kg/m2')
      check_climo_units = KG_M2
    case ('kg/m^2')
      check_climo_units = KG_M2
    case ('kg/m**2')
      check_climo_units = KG_M2
    case ('kg m^-2')
      check_climo_units = KG_M2
    case ('kg m**-2')
      check_climo_units = KG_M2
    end select

  end function check_climo_units
!
!#######################################################################
!
  !> Registers diagnostics of the column integrals of the climatology fields [kg/m2].
  !>
  !> For each field, module `climo` has the column integral on the climatology grid, as
  !> read from the file, and module `hinterp` the column integral after the horizontal
  !> interpolation to the model grid; both diagnose the interpolation. Also sets up a
  !> domain decomposition of the climatology grid for sending the `climo` fields.
  subroutine init_clim_diag(clim_type, mod_axes, init_time)
    type(interpolate_type), intent(inout)  :: clim_type  !! set up by `interpolator_init`
    integer, intent(in)     :: mod_axes(:)  !! diagnostic axes of the model grid (lon, lat, ...)
    type(time_type), intent(in)     :: init_time  !! model initial time

    integer :: axes(2), nxd, nyd, ndivs, i
    type(domain2d) :: domain
    integer :: domain_layout(2), iscomp, iecomp, jscomp, jecomp

    if (.not. module_is_initialized .or. .not. associated(clim_type%lon)) &
      call mpp_error(FATAL, "init_clim_diag : You must call interpolator_init before calling init_clim_diag")

    ndivs = mpp_npes()
    nxd = size(clim_type%lon(:))
    nyd = size(clim_type%lat(:))

! Define the domain decomposition of the climatology file. This may be (probably is) different from the model domain.
    call mpp_define_layout((/1, nxd, 1, nyd/), ndivs, domain_layout)
    call mpp_define_domains((/1, nxd, 1, nyd/), domain_layout, domain, xhalo=0, yhalo=0)
    call mpp_get_compute_domain(domain, iscomp, iecomp, jscomp, jecomp)
    axes(1) = diag_axis_init(clim_type%file_name(1:5)//'x', clim_type%lon, units='degrees', cart_name='x', domain2=domain)
    axes(2) = diag_axis_init(clim_type%file_name(1:5)//'y', clim_type%lat, units='degrees', cart_name='y', domain2=domain)
    clim_type%is = iscomp
    clim_type%ie = iecomp
    clim_type%js = jscomp
    clim_type%je = jecomp

!init_time = set_date(1980,1,1,0,0,0)

    if ((num_clim_diag + size(clim_type%field_name(:))) .gt. max_diag_fields) &
      call mpp_error(FATAL, "init_clim_diag : Trying to set up too many diagnostic fields for the climatology data")
    do i = 1, size(clim_type%field_name(:))
      climo_diag_name(i + num_clim_diag) = clim_type%field_name(i)
      climo_diag_id(i + num_clim_diag) = register_diag_field('climo', clim_type%field_name(i), axes(1:2), init_time, &
                                                             'column integral of '//trim(clim_type%field_name(i))// &
                                                             ' (climatology grid)', &
                                                             'kg/m2', missing_value)
      hinterp_id(i + num_clim_diag) = register_diag_field('hinterp', clim_type%field_name(i), mod_axes(1:2), init_time, &
                                                          'column integral of '//trim(clim_type%field_name(i))// &
                                                          ' (interpolated to the model grid)', &
                                                          'kg/m2', missing_value)
    end do
! Total number of climatology diagnostics (num_clim_diag). This can be from multiple climatology fields with different spatial axes.
! It is simply a holder for the diagnostic indices.
    num_clim_diag = num_clim_diag + size(clim_type%field_name(:))

    clim_diag_initialized = .true.

  end subroutine init_clim_diag
!

!---------------------------------------------------------------------

  !> Returns all the fields of an interpolate type, interpolated to the model grid, time
  !> and levels (see `interpolator`).
  subroutine interpolator_4D(clim_type, Time, phalf, interp_data, &
                             field_name, is, js, clim_units)
    type(interpolate_type), intent(inout)  :: clim_type  !! set up by `interpolator_init`
    character(len=*), intent(in)  :: field_name
    !! name of one of the fields (all the fields are interpolated)
    type(time_type), intent(in)  :: Time  !! model time to interpolate to
    real, dimension(:, :, :), intent(in)  :: phalf  !! model pressure at half levels [Pa]
    real, dimension(:, :, :, :), intent(out) :: interp_data
    !! the interpolated fields (last index: field)
    integer, intent(in), optional :: is, js
    !! `is`, `js`: indices of the first point of the physics window (default 1)
    character(len=*), intent(out), optional :: clim_units  !! units of the first field
    real :: tweight, tweight1, tweight2, tweight3
    integer :: taum, taup, ilon
    real :: hinterp_data(size(interp_data, 1), size(interp_data, 2), size(clim_type%levs(:)), size(clim_type%field_name(:)))
    real :: p_fact(size(interp_data, 1), size(interp_data, 2))
    real :: col_data(size(interp_data, 1), size(interp_data, 2), &
                     size(clim_type%field_name(:)))
    real :: pclim(size(clim_type%halflevs(:)))
    integer :: istart, iend, jstart, jend
    logical :: result, found
    logical :: found_field = .false.
    integer :: modyear, modmonth, modday, modhour, modminute, modsecond
    integer :: climyear, climmonth, climday, climhour, climminute, climsecond
    integer :: year1, month1, day, hour, minute, second
    integer :: taum1, taup1, taum2, taup2, climatology, m
    type(time_type) :: clim_datem, clim_datep, mod_time, prev_clim_time, t_prev, t_next
    type(time_type), dimension(2) :: month
    integer :: indexm, indexp, yearm, yearp
    integer :: i, j, k, n, itaum, itaup

    if (.not. module_is_initialized .or. .not. associated(clim_type%lon)) &
      call mpp_error(FATAL, "interpolator_3D : You must call interpolator_init before calling interpolator")

    do n = 2, size(clim_type%field_name(:))
      if (clim_type%vert_interp(n) /= clim_type%vert_interp(n - 1) .or. &
          clim_type%out_of_bounds(n) /= clim_type%out_of_bounds(n - 1)) then
        if (mpp_pe() == mpp_root_pe()) then
          print *, 'processing file '//trim(clim_type%file_name)
        end if
        call mpp_error(FATAL, 'mima_interpolator_mod: &
                &cannot use 4D interface to interpolator for this file')
      end if
    end do

    istart = 1
    if (present(is)) istart = is
    iend = istart - 1 + size(interp_data, 1)

    jstart = 1
    if (present(js)) jstart = js
    jend = jstart - 1 + size(interp_data, 2)

    do i = 1, size(clim_type%field_name(:))
!++lwh
      if (field_name == clim_type%field_name(i)) then
!--lwh
        found_field = .true.
        exit
      end if
    end do
    i = 1

    if (present(clim_units)) then
      clim_units = clim_type%field_type(i)%units
      clim_units = chomp(clim_units)
    end if
    if (size(clim_type%time_slice(:)) == 1) then
      taum = 1; taup = 1; tweight = 0.0   ! a single record, used at all times
    else if (size(clim_type%time_slice(:)) .le. 12) then
      call time_interp(Time, clim_type%time_slice, tweight, taum, taup, modtime=YEAR)
    else
      call time_interp(Time, clim_type%time_slice, tweight, taum, taup)
    end if

    if (clim_type%TIME_FLAG .ne. LINEAR .or. read_all_on_init) then
      itaum = taum
      itaup = taup
    end if

    if (clim_type%TIME_FLAG .eq. BILINEAR) then
      ! Check if delta-time is greater than delta of first two climatology time-slices.
      if ((Time - clim_type%time_slice(taum)) > (clim_type%time_slice(2) - clim_type%time_slice(1)) .or. &
          (clim_type%time_slice(taup) - Time) > (clim_type%time_slice(2) - clim_type%time_slice(1))) then
        ! The difference between the model time and the last climatology time-slice previous to the model time.
        ! We need 2 time levels.
        itaum = 0
        itaup = 0
        ! Assume this is monthly data. So we need to get the data applicable to the model date but substitute
        ! the climatology year into the appropriate place.

        ! We need to get the previous months data for the climatology year before
        ! and after the model year.
        call get_date(Time, modyear, modmonth, modday, modhour, modminute, modsecond)
        call get_date(clim_type%time_slice(taum), climyear, climmonth, climday, climhour, climminute, climsecond)

        climatology = 1
        do m = 1, size(clim_type%clim_times(:, :), 2)
          !Assume here that a climatology is for 1 year and consists of 12 months starting in January.
          call get_date(clim_type%clim_times(1, m), year1, month1, day, hour, minute, second)
          if (year1 == climyear) climatology = m
        end do
        do m = 1, 12
          !Find which month we are trying to look at and set clim_date[mp] to the dates spanning that.
          call get_date(clim_type%clim_times(m, climatology), year1, month1, day, hour, minute, second)
          if (month1 == modmonth) then
!RSHBUGFX   if ( modday <= day ) then
            if (modday < day) then
              indexm = m - 1; indexp = m
            else
              indexm = m; indexp = m + 1
            end if
          end if

        end do
        if (indexm == 0) then
          indexm = 12
          yearm = modyear - 1
        else
          yearm = modyear
        end if
        call get_date(clim_type%time_slice(indexm + (climatology - 1)*12), &
                      climyear, climmonth, climday, climhour, climminute, climsecond)
        month(1) = set_date(yearm, indexm, climday, climhour, climminute, climsecond)
        if (indexp == 13) then
          indexp = 1
          yearp = modyear + 1
        else
          yearp = modyear
        end if
        call get_date(clim_type%time_slice(indexp + (climatology - 1)*12), &
                      climyear, climmonth, climday, climhour, climminute, climsecond)
        month(2) = set_date(yearp, indexp, climday, climhour, climminute, climsecond)

        call time_interp(Time, month, tweight3, taum, taup) ! tweight3 is the time weight between the months.

        month(1) = clim_type%time_slice(indexm + (climatology - 1)*12)
        month(2) = clim_type%time_slice(indexm + climatology*12)
        call get_date(month(1), climyear, climmonth, climday, climhour, climminute, climsecond)
        t_prev = set_date(yearm, climmonth, climday, climhour, climminute, climsecond)
        call time_interp(t_prev, month, tweight1, taum, taup) ! tweight1 is the time weight between the climatology years.
        month(1) = clim_type%time_slice(indexp + (climatology - 1)*12)
        month(2) = clim_type%time_slice(indexp + climatology*12)
        call get_date(month(1), climyear, climmonth, climday, climhour, climminute, climsecond)
        t_next = set_date(yearp, climmonth, climday, climhour, climminute, climsecond)
        call time_interp(t_next, month, tweight2, taum, taup) ! tweight1 is the time weight between the climatology years.

        if (indexm == clim_type%indexm(1) .and. &
            indexp == clim_type%indexp(1) .and. &
            climatology == clim_type%climatology(1)) then
        else
          clim_type%indexm(:) = indexm
          clim_type%indexp(:) = indexp
          clim_type%climatology(:) = climatology
          do i = 1, size(clim_type%field_name(:))
            call read_data(clim_type, clim_type%field_type(i), &
                           clim_type%pmon_pyear(:, :, :, i), &
                           clim_type%indexm(i) + (clim_type%climatology(i) - 1)*12, i, Time)
! Read the data for the next month in the previous climatology.
            call read_data(clim_type, clim_type%field_type(i), &
                           clim_type%nmon_pyear(:, :, :, i), &
                           clim_type%indexp(i) + (clim_type%climatology(i) - 1)*12, i, Time)
            call read_data(clim_type, clim_type%field_type(i), &
                           clim_type%pmon_nyear(:, :, :, i), &
                           clim_type%indexm(i) + clim_type%climatology(i)*12, i, Time)
            call read_data(clim_type, clim_type%field_type(i), &
                           clim_type%nmon_nyear(:, :, :, i), &
                           clim_type%indexp(i) + clim_type%climatology(i)*12, i, Time)
          end do
        end if

      else ! We are within a climatology data set

        do i = 1, size(clim_type%field_name(:))
          if (taum /= clim_type%time_init(i, 1) .or. &
              taup /= clim_type%time_init(i, 2)) then

            call read_data(clim_type, clim_type%field_type(i), &
                           clim_type%pmon_pyear(:, :, :, i), taum, i, Time)
! Read the data for the next month in the previous climatology.
            call read_data(clim_type, clim_type%field_type(i), &
                           clim_type%nmon_pyear(:, :, :, i), taup, i, Time)
            clim_type%time_init(i, 1) = taum
            clim_type%time_init(i, 2) = taup
          end if
        end do

!       clim_type%pmon_nyear = 0.0
!       clim_type%nmon_nyear = 0.0

! set to zero so when next return to bilinear section will be sure to
! have proper data (relevant when running fixed_year case for more than
! one year in a single job)
        clim_type%indexm(:) = 0
        clim_type%indexp(:) = 0
        clim_type%climatology(:) = 0

!       tweight3 = 0.0 ! This makes [pn]mon_nyear irrelevant. Set them to 0 to test.
        tweight1 = 0.0
        tweight2 = 0.0
        tweight3 = tweight
      end if
    end if   !(BILINEAR)

    if (clim_type%TIME_FLAG .eq. LINEAR .and. &
        (.not. read_all_on_init)) then
! We need 2 time levels. Check we have the correct data.
      itaum = 0
      itaup = 0
      do n = 1, size(clim_type%time_init, 2)
        if (clim_type%time_init(1, n) .eq. taum) itaum = n
        if (clim_type%time_init(1, n) .eq. taup) itaup = n
      end do

      if (itaum .eq. 0 .and. itaup .eq. 0) then
!Neither time is set so we need to read 2 time slices.
!Set up
! field(:,:,:,1) as the previous time slice.
! field(:,:,:,2) as the next time slice.
        do i = 1, size(clim_type%field_name(:))
          call read_data(clim_type, clim_type%field_type(i), clim_type%data(:, :, :, 1, i), taum, i, Time)
          clim_type%time_init(i, 1) = taum
          itaum = 1
          call read_data(clim_type, clim_type%field_type(i), clim_type%data(:, :, :, 2, i), taup, i, Time)
          clim_type%time_init(i, 2) = taup
          itaup = 2
        end do
      end if ! itaum.eq.itaup.eq.0
      if (itaum .eq. 0 .and. itaup .ne. 0) then
! Can't think of a situation where we would have the next time level but not the previous.
        call mpp_error(FATAL, 'interpolator_3D : No data from the previous climatology time &
                                & but we have the next time. How did this happen?')
      end if
      if (itaum .ne. 0 .and. itaup .eq. 0) then
!We have the previous time step but not the next time step data
        itaup = 1
        if (itaum .eq. 1) itaup = 2
        do i = 1, size(clim_type%field_name(:))
          call read_data(clim_type, clim_type%field_type(i), clim_type%data(:, :, :, itaup, i), taup, i, Time)
          clim_type%time_init(i, itaup) = taup
        end do
      end if

    end if! TIME_FLAG

    select case (clim_type%TIME_FLAG)
    case (LINEAR, SEASONAL)
      do n = 1, size(clim_type%field_name(:))
        hinterp_data(:, :, :, n) = (1 - tweight)* &
                                   clim_type%data(istart:iend, jstart:jend, :, itaum, n) + &
                                   tweight* &
                                   clim_type%data(istart:iend, jstart:jend, :, itaup, n)
      end do
    case (BILINEAR)
      do n = 1, size(clim_type%field_name(:))
        hinterp_data(:, :, :, n) = (1 - tweight1)*(1 - tweight3)* &
                                   clim_type%pmon_pyear(istart:iend, jstart:jend, :, n) + &
                                   (1 - tweight2)*tweight3* &
                                   clim_type%nmon_pyear(istart:iend, jstart:jend, :, n) + &
                                   tweight1*(1 - tweight3)* &
                                   clim_type%pmon_nyear(istart:iend, jstart:jend, :, n) + &
                                   tweight2*tweight3* &
                                   clim_type%nmon_nyear(istart:iend, jstart:jend, :, n)

      end do

    end select

    select case (clim_type%level_type)
    case (PRESSURE)
      p_fact = 1.0
    case (SIGMA)
      p_fact = maxval(phalf, 3)! max pressure in the column !(:,:,size(phalf,3))
    end select

    col_data(:, :, :) = 0.0
    do i = 1, size(clim_type%field_name(:))

      select case (clim_type%mr(i))
      case (NO_CONV)
        do k = 1, size(hinterp_data, 3)
          col_data(:, :, i) = col_data(:, :, i) + hinterp_data(:, :, k, i)* &
                              (clim_type%halflevs(k + 1) - clim_type%halflevs(k))/grav
        end do

      case (KG_M2)
        do k = 1, size(hinterp_data, 3)
          col_data(:, :, i) = col_data(:, :, i) + hinterp_data(:, :, k, i)
          hinterp_data(:, :, k, i) = hinterp_data(:, :, k, i)/ &
                                     ((clim_type%halflevs(k + 1) - clim_type%halflevs(k))*p_fact)
        end do
      end select
    end do

    do i = 1, size(clim_type%field_name(:))
      found = .false.
      do j = 1, size(climo_diag_name(:))
        if (climo_diag_name(j) .eq. clim_type%field_name(i)) then
          found = .true.
          exit
        end if
      end do

      if (found) then
        if (hinterp_id(j) > 0) then
          result = send_data(hinterp_id(j), col_data(:, :, i), Time)
        end if
      end if

    end do

    i = 1

!++lwh
    do j = 1, size(phalf, 2)
      do ilon = 1, size(phalf, 1)
        pclim = p_fact(ilon, j)*clim_type%halflevs
        if (maxval(phalf(ilon, j, :)) > maxval(pclim)) then
          if (verbose > 3) then
            call mpp_error(NOTE, "Interpolator: model surface pressure&
                                & is greater than climatology surface pressure for " &
                                //trim(clim_type%file_name))
          end if
          select case (clim_type%out_of_bounds(i))
          case (CONSTANT)
            pclim(maxloc(pclim)) = maxval(phalf(ilon, j, :))
!           case(ZERO)
!              pclim( maxloc(pclim)) = 0
          end select
        end if
        if (minval(phalf(ilon, j, :)) < minval(pclim)) then
          if (verbose > 3) then
            call mpp_error(NOTE, "Interpolator: model top pressure&
                                & is less than climatology top pressure for " &
                                //trim(clim_type%file_name))
          end if
          select case (clim_type%out_of_bounds(i))
          case (CONSTANT)
            pclim(minloc(pclim)) = minval(phalf(ilon, j, :))
!           case(ZERO)
!              pclim( maxloc(pclim)) = 0
          end select
        end if
        select case (clim_type%vert_interp(i))
        case (INTERP_WEIGHTED_P)
          call interp_weighted_scalar(pclim, phalf(ilon, j, :), hinterp_data(ilon, j, :, :), interp_data(ilon, j, :, :))
        case (INTERP_LINEAR_P)
          do n = 1, size(clim_type%field_name(:))
            call interp_linear(pclim, phalf(ilon, j, :), hinterp_data(ilon, j, :, n), interp_data(ilon, j, :, n))
          end do
!        case(INTERP_LOG)
        end select
      end do
    end do

!--lwh
    do i = 1, size(clim_type%field_name(:))

      select case (clim_type%mr(i))
      case (KG_M2)
        do k = 1, size(interp_data, 3)
          interp_data(:, :, k, i) = interp_data(:, :, k, i)*(phalf(:, :, k + 1) - phalf(:, :, k))
        end do
      end select

    end do

    if (.not. found_field) then !field name is not in interpolator file.ERROR.
      call mpp_error(FATAL, "Interpolator: the field name is not contained in this &
                       &intepolate_type: "//trim(field_name))
    end if
  end subroutine interpolator_4D
!
!#######################################################################
!#######################################################################
!
  !> Returns a 3-D field interpolated to the model grid, time and levels (see
  !> `interpolator`).
  subroutine interpolator_3D(clim_type, Time, phalf, interp_data, field_name, is, js, clim_units)
    type(interpolate_type), intent(inout)  :: clim_type  !! set up by `interpolator_init`
    character(len=*), intent(in)  :: field_name  !! name of the field to interpolate
    type(time_type), intent(in)  :: Time  !! model time to interpolate to
    real, dimension(:, :, :), intent(in)  :: phalf  !! model pressure at half levels [Pa]
    real, dimension(:, :, :), intent(out) :: interp_data  !! the interpolated field
    integer, intent(in), optional :: is, js
    !! `is`, `js`: indices of the first point of the physics window (default 1)
    character(len=*), intent(out), optional :: clim_units  !! units of the field
    real :: tweight, tweight1, tweight2, tweight3
    integer :: taum, taup, ilon
    real :: hinterp_data(size(interp_data, 1), size(interp_data, 2), size(clim_type%levs(:)))
    real :: p_fact(size(interp_data, 1), size(interp_data, 2))
    real :: col_data(size(interp_data, 1), size(interp_data, 2))
    real :: pclim(size(clim_type%halflevs(:)))
    integer :: istart, iend, jstart, jend
    logical :: result, found
    logical :: found_field = .false.
    integer :: modyear, modmonth, modday, modhour, modminute, modsecond
    integer :: climyear, climmonth, climday, climhour, climminute, climsecond
    integer :: year1, month1, day, hour, minute, second
    integer :: taum1, taup1, taum2, taup2, climatology, m
    type(time_type) :: clim_datem, clim_datep, mod_time, prev_clim_time, t_prev, t_next
    type(time_type), dimension(2) :: month
    integer :: indexm, indexp, yearm, yearp
    integer :: i, j, k, itaum, itaup, n

    if (.not. module_is_initialized .or. .not. associated(clim_type%lon)) &
      call mpp_error(FATAL, "interpolator_3D : You must call interpolator_init before calling interpolator")

    istart = 1
    if (present(is)) istart = is
    iend = istart - 1 + size(interp_data, 1)

    jstart = 1
    if (present(js)) jstart = js
    jend = jstart - 1 + size(interp_data, 2)

    do i = 1, size(clim_type%field_name(:))
!++lwh
      if (field_name == clim_type%field_name(i)) then
!--lwh
        found_field = .true.
        if (present(clim_units)) then
          clim_units = clim_type%field_type(i)%units
          clim_units = chomp(clim_units)
        end if
        if (size(clim_type%time_slice(:)) == 1) then
          taum = 1; taup = 1; tweight = 0.0   ! a single record, used at all times
        else if (size(clim_type%time_slice(:)) .le. 12) then
          call time_interp(Time, clim_type%time_slice, tweight, taum, taup, modtime=YEAR)
        else
          call time_interp(Time, clim_type%time_slice, tweight, taum, taup)
        end if

!   if(clim_type%TIME_FLAG .ne. LINEAR ) then
        if (clim_type%TIME_FLAG .ne. LINEAR .or. read_all_on_init) then
          itaum = taum
          itaup = taup
        end if

        if (clim_type%TIME_FLAG .eq. BILINEAR) then
          ! Check if delta-time is greater than delta of first two climatology time-slices.
          if ((Time - clim_type%time_slice(taum)) > (clim_type%time_slice(2) - clim_type%time_slice(1)) .or. &
              (clim_type%time_slice(taup) - Time) > (clim_type%time_slice(2) - clim_type%time_slice(1))) then
            ! The difference between the model time and the last climatology time-slice previous to the model time.
            ! We need 2 time levels.
            itaum = 0
            itaup = 0
            ! Assume this is monthly data. So we need to get the data applicable to the model date but substitute
            ! the climatology year into the appropriate place.

            ! We need to get the previous months data for the climatology year before
            ! and after the model year.
            call get_date(Time, modyear, modmonth, modday, modhour, modminute, modsecond)
            call get_date(clim_type%time_slice(taum), climyear, climmonth, climday, climhour, climminute, climsecond)

            climatology = 1
            do m = 1, size(clim_type%clim_times(:, :), 2)
              !Assume here that a climatology is for 1 year and consists of 12 months starting in January.
              call get_date(clim_type%clim_times(1, m), year1, month1, day, hour, minute, second)
              if (year1 == climyear) climatology = m
            end do
            do m = 1, 12
              !Find which month we are trying to look at and set clim_date[mp] to the dates spanning that.
              call get_date(clim_type%clim_times(m, climatology), year1, month1, day, hour, minute, second)
              if (month1 == modmonth) then
!RSHBUGFX   if ( modday <= day ) then
                if (modday < day) then
                  indexm = m - 1; indexp = m
                else
                  indexm = m; indexp = m + 1
                end if
              end if

            end do
            if (indexm == 0) then
              indexm = 12
              yearm = modyear - 1
            else
              yearm = modyear
            end if
            call get_date(clim_type%time_slice(indexm + (climatology - 1)*12), &
                          climyear, climmonth, climday, climhour, climminute, climsecond)
            month(1) = set_date(yearm, indexm, climday, climhour, climminute, climsecond)
            if (indexp == 13) then
              indexp = 1
              yearp = modyear + 1
            else
              yearp = modyear
            end if
            call get_date(clim_type%time_slice(indexp + (climatology - 1)*12), &
                          climyear, climmonth, climday, climhour, climminute, climsecond)
            month(2) = set_date(yearp, indexp, climday, climhour, climminute, climsecond)

            call time_interp(Time, month, tweight3, taum, taup) ! tweight3 is the time weight between the months.

            month(1) = clim_type%time_slice(indexm + (climatology - 1)*12)
            month(2) = clim_type%time_slice(indexm + climatology*12)
            call get_date(month(1), climyear, climmonth, climday, climhour, climminute, climsecond)
            t_prev = set_date(yearm, climmonth, climday, climhour, climminute, climsecond)
            call time_interp(t_prev, month, tweight1, taum, taup) ! tweight1 is the time weight between the climatology years.

            month(1) = clim_type%time_slice(indexp + (climatology - 1)*12)
            month(2) = clim_type%time_slice(indexp + climatology*12)
            call get_date(month(1), climyear, climmonth, climday, climhour, climminute, climsecond)
            t_next = set_date(yearp, climmonth, climday, climhour, climminute, climsecond)
            call time_interp(t_next, month, tweight2, taum, taup) ! tweight1 is the time weight between the climatology years.

            if (indexm == clim_type%indexm(i) .and. &
                indexp == clim_type%indexp(i) .and. &
                climatology == clim_type%climatology(i)) then
            else
              clim_type%indexm(i) = indexm
              clim_type%indexp(i) = indexp
              clim_type%climatology(i) = climatology
              call read_data(clim_type, clim_type%field_type(i), &
                             clim_type%pmon_pyear(:, :, :, i), &
                             clim_type%indexm(i) + (clim_type%climatology(i) - 1)*12, i, Time)
! Read the data for the next month in the previous climatology.
              call read_data(clim_type, clim_type%field_type(i), &
                             clim_type%nmon_pyear(:, :, :, i), &
                             clim_type%indexp(i) + (clim_type%climatology(i) - 1)*12, i, Time)
              call read_data(clim_type, clim_type%field_type(i), &
                             clim_type%pmon_nyear(:, :, :, i), &
                             clim_type%indexm(i) + clim_type%climatology(i)*12, i, Time)
              call read_data(clim_type, clim_type%field_type(i), &
                             clim_type%nmon_nyear(:, :, :, i), &
                             clim_type%indexp(i) + clim_type%climatology(i)*12, i, Time)
            end if

          else ! We are within a climatology data set

            if (taum /= clim_type%time_init(i, 1) .or. &
                taup /= clim_type%time_init(i, 2)) then

              call read_data(clim_type, clim_type%field_type(i), clim_type%pmon_pyear(:, :, :, i), taum, i, Time)
! Read the data for the next month in the previous climatology.
              call read_data(clim_type, clim_type%field_type(i), clim_type%nmon_pyear(:, :, :, i), taup, i, Time)
!RSHbug   clim_type%pmon_nyear = 0.0
!RSHbug   clim_type%nmon_nyear = 0.0

!         clim_type%pmon_nyear(:,:,:,i) = 0.0
!         clim_type%nmon_nyear(:,:,:,i) = 0.0

! set to zero so when next return to bilinear section will be sure to
! have proper data (relevant when running fixed_year case for more than
! one year in a single job)
              clim_type%indexm(i) = 0
              clim_type%indexp(i) = 0
              clim_type%climatology(i) = 0

              clim_type%time_init(i, 1) = taum
              clim_type%time_init(i, 2) = taup
            end if
!       tweight3 = 0.0 ! This makes [pn]mon_nyear irrelevant. Set them to 0 to test.
            tweight1 = 0.0; tweight2 = 0.0
            tweight3 = tweight
          end if

        end if ! (BILINEAR)

        if (clim_type%TIME_FLAG .eq. LINEAR .and. &
            (.not. read_all_on_init)) then
! We need 2 time levels. Check we have the correct data.
          itaum = 0
          itaup = 0
          do n = 1, size(clim_type%time_init, 2)
            if (clim_type%time_init(i, n) .eq. taum) itaum = n
            if (clim_type%time_init(i, n) .eq. taup) itaup = n
          end do

          if (itaum .eq. 0 .and. itaup .eq. 0) then
!Neither time is set so we need to read 2 time slices.
!Set up
! field(:,:,:,1) as the previous time slice.
! field(:,:,:,2) as the next time slice.
            call read_data(clim_type, clim_type%field_type(i), clim_type%data(:, :, :, 1, i), taum, i, Time)
            clim_type%time_init(i, 1) = taum
            itaum = 1
            call read_data(clim_type, clim_type%field_type(i), clim_type%data(:, :, :, 2, i), taup, i, Time)
            clim_type%time_init(i, 2) = taup
            itaup = 2
          end if ! itaum.eq.itaup.eq.0
          if (itaum .eq. 0 .and. itaup .ne. 0) then
! Can't think of a situation where we would have the next time level but not the previous.
            call mpp_error(FATAL, 'interpolator_3D : No data from the previous climatology time &
                                    & but we have the next time. How did this happen?')
          end if
          if (itaum .ne. 0 .and. itaup .eq. 0) then
!We have the previous time step but not the next time step data
            itaup = 1
            if (itaum .eq. 1) itaup = 2
            call read_data(clim_type, clim_type%field_type(i), clim_type%data(:, :, :, itaup, i), taup, i, Time)
            clim_type%time_init(i, itaup) = taup
          end if

        end if! TIME_FLAG

        select case (clim_type%TIME_FLAG)
        case (LINEAR, SEASONAL)
          hinterp_data = (1 - tweight)*clim_type%data(istart:iend, jstart:jend, :, itaum, i) + &
                         tweight*clim_type%data(istart:iend, jstart:jend, :, itaup, i)
        case (BILINEAR)
          hinterp_data = &
            (1 - tweight1)*(1 - tweight3)*clim_type%pmon_pyear(istart:iend, jstart:jend, :, i) + &
            (1 - tweight2)*tweight3*clim_type%nmon_pyear(istart:iend, jstart:jend, :, i) + &
            tweight1*(1 - tweight3)*clim_type%pmon_nyear(istart:iend, jstart:jend, :, i) + &
            tweight2*tweight3*clim_type%nmon_nyear(istart:iend, jstart:jend, :, i)

        end select

        select case (clim_type%level_type)
        case (PRESSURE)
          p_fact = 1.0
        case (SIGMA)
          p_fact = maxval(phalf, 3)! max pressure in the column !(:,:,size(phalf,3))
        end select

        col_data(:, :) = 0.0
        select case (clim_type%mr(i))
        case (NO_CONV)
          do k = 1, size(hinterp_data, 3)
            col_data(:, :) = col_data(:, :) + hinterp_data(:, :, k)* &
                             (clim_type%halflevs(k + 1) - clim_type%halflevs(k))/grav
          end do

        case (KG_M2)
          do k = 1, size(hinterp_data, 3)
            col_data(:, :) = col_data(:, :) + hinterp_data(:, :, k)
            hinterp_data(:, :, k) = hinterp_data(:, :, k)/ &
                                    ((clim_type%halflevs(k + 1) - clim_type%halflevs(k))*p_fact)
          end do
        end select

        found = .false.
        do j = 1, size(climo_diag_name(:))
          if (climo_diag_name(j) .eq. clim_type%field_name(i)) then
            found = .true.
            exit
          end if
        end do

        if (found) then
          if (hinterp_id(j) > 0) then
            result = send_data(hinterp_id(j), col_data, Time)
          end if
        end if

!++lwh
        do j = 1, size(phalf, 2)
          do ilon = 1, size(phalf, 1)
            pclim = p_fact(ilon, j)*clim_type%halflevs
            if (maxval(phalf(ilon, j, :)) > maxval(pclim)) then
              if (verbose > 3) then
                call mpp_error(NOTE, "Interpolator: model surface pressure&
                                    & is greater than climatology surface pressure for " &
                                    //trim(clim_type%file_name))
              end if
              select case (clim_type%out_of_bounds(i))
              case (CONSTANT)
                pclim(maxloc(pclim)) = maxval(phalf(ilon, j, :))
!           case(ZERO)
!              pclim( maxloc(pclim)) = 0
              end select
            end if
            if (minval(phalf(ilon, j, :)) < minval(pclim)) then
              if (verbose > 3) then
                call mpp_error(NOTE, "Interpolator: model top pressure&
                                    & is less than climatology top pressure for " &
                                    //trim(clim_type%file_name))
              end if
              select case (clim_type%out_of_bounds(i))
              case (CONSTANT)
                pclim(minloc(pclim)) = minval(phalf(ilon, j, :))
!           case(ZERO)
!              pclim( maxloc(pclim)) = 0
              end select
            end if
            select case (clim_type%vert_interp(i))
            case (INTERP_WEIGHTED_P)
              call interp_weighted_scalar(pclim, phalf(ilon, j, :), hinterp_data(ilon, j, :), interp_data(ilon, j, :))
            case (INTERP_LINEAR_P)
              call interp_linear(pclim, phalf(ilon, j, :), hinterp_data(ilon, j, :), interp_data(ilon, j, :))
!        case(INTERP_LOG)
            end select
          end do
        end do

!--lwh

        select case (clim_type%mr(i))
        case (KG_M2)
          do k = 1, size(interp_data, 3)
            interp_data(:, :, k) = interp_data(:, :, k)*(phalf(:, :, k + 1) - phalf(:, :, k))
          end do
        end select

      end if !field_name
    end do !End of i loop
    if (.not. found_field) then !field name is not in interpolator file.ERROR.
      call mpp_error(FATAL, "Interpolator: the field name is not contained in this &
                       &intepolate_type: "//trim(field_name))
    end if
  end subroutine interpolator_3D
!
!#######################################################################
!
!++lwh
  !> Returns a 2-D field interpolated to the model grid and time (see `interpolator`).
  subroutine interpolator_2D(clim_type, Time, interp_data, field_name, is, js, clim_units)

    type(interpolate_type), intent(inout)  :: clim_type  !! set up by `interpolator_init`
    character(len=*), intent(in)     :: field_name  !! name of the field to interpolate
    type(time_type), intent(in)     :: Time  !! model time to interpolate to
    real, dimension(:, :), intent(out)    :: interp_data  !! the interpolated field
    integer, intent(in), optional :: is, js
    !! `is`, `js`: indices of the first point of the physics window (default 1)
    character(len=*), intent(out), optional :: clim_units  !! units of the field
    real :: tweight, tweight1, tweight2
    integer :: taum, taup, ilon
    real :: hinterp_data(size(interp_data, 1), size(interp_data, 2), size(clim_type%levs(:)))
    real :: p_fact(size(interp_data, 1), size(interp_data, 2))
    real :: col_data(size(interp_data, 1), size(interp_data, 2))
    integer :: istart, iend, jstart, jend
    logical :: result, found
    logical :: found_field = .false.
    integer :: modyear, modmonth, modday, modhour, modminute, modsecond
    integer :: climyear, climmonth, climday, climhour, climminute, climsecond
    integer :: taum1, taup1, taum2, taup2
    type(time_type) :: clim_datem, clim_datep
    integer :: j, k, i, itaum, itaup, n

    if (.not. module_is_initialized .or. .not. associated(clim_type%lon)) &
      call mpp_error(FATAL, "interpolator_2D : You must call interpolator_init before calling interpolator")

    istart = 1
    if (present(is)) istart = is
    iend = istart - 1 + size(interp_data, 1)

    jstart = 1
    if (present(js)) jstart = js
    jend = jstart - 1 + size(interp_data, 2)

    do i = 1, size(clim_type%field_name(:))
!++lwh
      if (field_name == clim_type%field_name(i)) then
!--lwh

        found_field = .true.

        if (present(clim_units)) then
          clim_units = clim_type%field_type(i)%units
          clim_units = chomp(clim_units)
        end if
        if (size(clim_type%time_slice(:)) == 1) then
          taum = 1; taup = 1; tweight = 0.0   ! a single record, used at all times
        else if (size(clim_type%time_slice(:)) .le. 12) then
          call time_interp(Time, clim_type%time_slice, tweight, taum, taup, modtime=YEAR)
        else
          call time_interp(Time, clim_type%time_slice, tweight, taum, taup)
        end if

! If the climatology file has seasonal, a split time-line or has all the data
! read in then enter this loop.
!
        if (clim_type%TIME_FLAG .ne. LINEAR .or. read_all_on_init) then
          itaum = taum
          itaup = taup
        end if

!    if(clim_type%TIME_FLAG .eq. BILINEAR ) then
!      ! Check if delta-time is greater than delta of first two climatology time-slices.
!      if ( (Time - clim_type%time_slice(taum) ) > ( clim_type%time_slice(2)- clim_type%time_slice(1) ) .or. &
!           (clim_type%time_slice(taup)  - Time) > ( clim_type%time_slice(2)- clim_type%time_slice(1) ) ) then
!      ! The difference between the model time and the last climatology time-slice previous to the model time.
!      ! We need 2 time levels. Check we have the correct data.
!        itaum=0
!        itaup=0
!      ! Assume this is monthly data. So we need to get the data applicable to the model date but substitute
!      ! the climatology year into the appropriate place.
!
!        call get_date(Time, modyear, modmonth, modday, modhour, modminute, modsecond)
!        call get_date(clim_type%time_slice(taum), climyear, climmonth, climday, climhour, climminute, climsecond)
!        clim_datem = set_date(climyear, modmonth, modday, modhour, modminute, modsecond)
!        call time_interp(clim_datem, clim_type%time_slice, tweight1, taum1, taup1 )
!
!
!        call get_date(clim_type%time_slice(taup), climyear, climmonth, climday, climhour, climminute, climsecond)
!        clim_datep = set_date(climyear, modmonth, modday, modhour, modminute, modsecond)
!
!
!      endif
!
!    endif

        if (clim_type%TIME_FLAG .eq. LINEAR .and. &
            (.not. read_all_on_init)) then
! We need 2 time levels. Check we have the correct data.
          itaum = 0
          itaup = 0
          do n = 1, size(clim_type%time_init, 2)
            if (clim_type%time_init(i, n) .eq. taum) itaum = n
            if (clim_type%time_init(i, n) .eq. taup) itaup = n
          end do

          if (itaum .eq. 0 .and. itaup .eq. 0) then
            !Neither time is set so we need to read 2 time slices.
            !Set up
            ! field(:,:,:,1) as the previous time slice.
            ! field(:,:,:,2) as the next time slice.
            call read_data(clim_type, clim_type%field_type(i), clim_type%data(:, :, :, 1, i), taum, i, Time)
            clim_type%time_init(i, 1) = taum
            itaum = 1
            call read_data(clim_type, clim_type%field_type(i), clim_type%data(:, :, :, 2, i), taup, i, Time)
            clim_type%time_init(i, 2) = taup
            itaup = 2
          end if ! itaum.eq.itaup.eq.0
          if (itaum .eq. 0 .and. itaup .ne. 0) then
            ! Can't think of a situation where we would have the next time level but not the previous.
            call mpp_error(FATAL, 'interpolator_2D : No data from the previous climatology time but we have&
                                & the next time. How did this happen?')
          end if
          if (itaum .ne. 0 .and. itaup .eq. 0) then
            !We have the previous time step but not the next time step data
            itaup = 1
            if (itaum .eq. 1) itaup = 2
            call read_data(clim_type, clim_type%field_type(i), clim_type%data(:, :, :, itaup, i), taup, i, Time)
            clim_type%time_init(i, itaup) = taup
          end if
        end if! TIME_FLAG .eq. LINEAR .and. (.not. read_all_on_init)

        select case (clim_type%TIME_FLAG)
        case (LINEAR, SEASONAL)
          hinterp_data = (1 - tweight)*clim_type%data(istart:iend, jstart:jend, :, itaum, i) &
                         + tweight*clim_type%data(istart:iend, jstart:jend, :, itaup, i)
        case (BILINEAR)
          call mpp_error(FATAL, 'interpolator_2D : data with gaps in time (from '// &
                         trim(clim_type%file_name)//') need the 3-D or 4-D interface')

        end select

        found = .false.
        do j = 1, size(climo_diag_name(:))
          if (climo_diag_name(j) .eq. clim_type%field_name(i)) then
            found = .true.
            exit
          end if
        end do

        if (found) then
          if (hinterp_id(j) > 0) then
            result = send_data(hinterp_id(j), hinterp_data, Time)
          end if
        end if

        interp_data(:, :) = hinterp_data(:, :, 1)

      end if !field_name
    end do !End of i loop

    if (.not. found_field) then !field name is not in interpolator file.ERROR.
      call mpp_error(FATAL, "Interpolator: the field name is not contained in this &
                       &intepolate_type: "//trim(field_name))
    end if
  end subroutine interpolator_2D
!--lwh
!
!#######################################################################
!
  !> Frees the memory of an interpolate type and closes its file.
  subroutine interpolator_end(clim_type)
    type(interpolate_type), intent(inout) :: clim_type  !! the interpolate type to free
    integer :: log_unit

    if (mpp_pe() == mpp_root_pe()) then
      write (stdlog(), '(/,(a))') 'Exiting interpolator, have a nice day ...'
    end if

    deallocate (clim_type%lat)
    deallocate (clim_type%lon)
    deallocate (clim_type%latb)
    deallocate (clim_type%lonb)
    deallocate (clim_type%levs)
    deallocate (clim_type%halflevs)
    call horiz_interp_del(clim_type%interph)
    deallocate (clim_type%time_slice)
    deallocate (clim_type%field_type)
    deallocate (clim_type%field_name)
    deallocate (clim_type%time_init)
    deallocate (clim_type%mr)
    if (associated(clim_type%data)) then
      deallocate (clim_type%data)
    end if
    if (associated(clim_type%pmon_pyear)) then
      deallocate (clim_type%pmon_pyear)
      deallocate (clim_type%pmon_nyear)
      deallocate (clim_type%nmon_nyear)
      deallocate (clim_type%nmon_pyear)
    end if

! RSH mod
    if (.not. (clim_type%TIME_FLAG .eq. LINEAR .and. &
               !     read_all_on_init)) .or. clim_type%TIME_FLAG .eq. BILINEAR  ) then
               read_all_on_init)) then
      call nc_check(nf90_close(clim_type%ncid), clim_type%file_name)
    end if

    module_is_initialized = .false.

  end subroutine interpolator_end
!
!#######################################################################
!
  !> Reads one record of a field and interpolates it horizontally to the model grid.
  subroutine read_data(clim_type, src_field, hdata, nt, i, Time)
!
!  INTENT IN
!    clim_type : The interpolate type which contains the data
!    src_field : The field type
!    nt        : The index of the time slice of the climatology that you wish to read.
!    i         : The index of the field name that you are trying to read. (optional)
!    Time      : The model time. Used for diagnostic purposes only. (optional)
!
!  INTENT OUT
!
!    hdata     : The horizontally interpolated climatology field. This
!                field will still be on the climatology vertical grid.
!
    type(interpolate_type), intent(in)  :: clim_type
    type(field_info_type), intent(in)  :: src_field
    integer, intent(in)  :: nt
    real, intent(out) :: hdata(:, :, :)
    integer, optional, intent(in)  :: i
    type(time_type), optional, intent(in)  :: Time

    integer   :: k, km, n
    integer   :: start(size(src_field%count))
! sjs
    real, allocatable :: climdata(:, :, :), climdata2(:, :, :), buf(:)

    allocate (climdata(size(clim_type%lon(:)), size(clim_type%lat(:)), &
                       size(clim_type%levs(:))))

!  read time level nt of the field and unpack it, as mpp_read did
    n = product(src_field%count)
    if (n > size(climdata)) call mpp_error(FATAL, 'interpolator read_data : '// &
                                           trim(src_field%name)//' in '//trim(clim_type%file_name)// &
                                           ' is larger than the climatology grid')
    allocate (buf(size(climdata)))
    buf = 0.0
    start = 1
    if (src_field%tdim > 0) start(src_field%tdim) = nt
    call nc_check(nf90_get_var(clim_type%ncid, src_field%varid, buf(1:n), start, src_field%count), &
                  trim(clim_type%file_name)//' '//src_field%name)
    buf(1:n) = buf(1:n)*src_field%scale + src_field%add
    climdata = reshape(buf, shape(climdata))
    deallocate (buf)

!  if vertical index increases upward, flip the data so that lowest
!  pressure level data is at index 1, rather than the highest pressure
!  level data. the indices themselves were previously flipped.
    if (clim_type%vertical_indices == INCREASING_UPWARD) then
      allocate (climdata2(size(clim_type%lon(:)), &
                          size(clim_type%lat(:)), &
                          size(clim_type%levs(:))))
      km = size(clim_type%levs(:))
      do k = 1, km
        climdata2(:, :, k) = climdata(:, :, km + 1 - k)
      end do
      climdata = climdata2
      deallocate (climdata2)
    end if

    call horiz_interp(clim_type%interph, climdata, hdata)
    if (clim_diag_initialized) &
      call diag_read_data(clim_type, climdata, i, Time)
    deallocate (climdata)

  end subroutine read_data
!
!#######################################################################
!
  !> Sends the `climo` column-integral diagnostic of a field read by `read_data`.
  subroutine diag_read_data(clim_type, model_data, i, Time)
!
! A routine to diagnose the data read in by read_data
!
!  INTENT IN
!    clim_type  : The interpolate type.
!    model_data : The data read in from file that is being diagnosed.
!    i          : The index of the field name that you are diagnosing.
!    Time       : The model time
!
    type(interpolate_type), intent(in) :: clim_type
    real, intent(in) :: model_data(:, :, :)
    integer, intent(in) :: i
    type(time_type), intent(in) :: Time

    integer :: j, k
    real :: col_data(size(model_data, 1), size(model_data, 2))
    logical :: result, found

    found = .false.
    do j = 1, size(climo_diag_name(:))
      if (climo_diag_name(j) .eq. clim_type%field_name(i)) then
        found = .true.
        exit
      end if
    end do

    if (found) then
      if (climo_diag_id(j) > 0) then
        col_data(:, :) = 0.0
        do k = 1, size(model_data, 3)
          col_data(:, :) = col_data(:, :) + &
                           model_data(:, :, k)* &
                           (clim_type%halflevs(k + 1) - clim_type%halflevs(k))/grav
        end do
        result = send_data(climo_diag_id(j), col_data(clim_type%is:clim_type%ie, clim_type%js:clim_type%je), Time)
      end if
    end if

  end subroutine diag_read_data
!
!#######################################################################
!
!++lwh
  !> Returns the number of fields and the field names of an interpolate type.
  subroutine query_interpolator(clim_type, nfields, field_names)
    type(interpolate_type), intent(in)                    :: clim_type  !! set up by `interpolator_init`
    integer, intent(out), optional                        :: nfields  !! number of fields
    character(len=*), dimension(:), intent(out), optional :: field_names  !! names of the fields

    if (present(nfields)) nfields = size(clim_type%field_name(:))
    if (present(field_names)) field_names = clim_type%field_name

  end subroutine query_interpolator
!--lwh
!
!#######################################################################
!
  !> Returns a string without a trailing CHAR(0), as read from netCDF files.
  function chomp(string)
!
! A function to remove CHAR(0) from the end of strings read from NetCDF files.
!
    character(len=*), intent(in) :: string
    character(len=64) :: chomp

    integer :: len

    len = len_trim(string)
    if (string(len:len) == char(0)) len = len - 1

    chomp = string(:len)

  end function chomp
!
!#######################################################################
!
  !> Reads the metadata of an open netCDF file: its axes, fields and times.
  subroutine get_file_info(ncid, file_name, axes, fields, time_values, ntime)
!
! Read the metadata of an open netCDF file, as mpp_io's mpp_read_meta did.
!
! INTENT OUT
!  axes        :: Every dimension, with the data and attributes of its coordinate variable.
!  fields      :: Every variable that is not a coordinate variable, in file order.
!  time_values :: The values of the record (time) coordinate variable.
!  ntime       :: The length of the record dimension, or -1 if there is none.
!
    integer, intent(in)  :: ncid
    character(len=*), intent(in)  :: file_name
    type(axis_info_type), allocatable, intent(out) :: axes(:)
    type(field_info_type), allocatable, intent(out) :: fields(:)
    real, allocatable, intent(out) :: time_values(:)
    integer, intent(out) :: ntime

    integer :: ndims, nvars, recdim, nf, i, j, k, nvdims, nvatts, dimid
    integer :: dimids(NF90_MAX_VAR_DIMS)
    character(len=128) :: name, attname, positive
    logical, allocatable :: isdim(:)

    call nc_check(nf90_inquire(ncid, ndims, nvars, unlimitedDimId=recdim), file_name)
    allocate (axes(ndims))
    do i = 1, ndims
      call nc_check(nf90_inquire_dimension(ncid, i, axes(i)%name, axes(i)%len), file_name)
      allocate (axes(i)%data(axes(i)%len))
      axes(i)%data = 0.0
    end do
! Without an unlimited dimension, a fixed-size dimension called time is the record dimension.
    if (recdim <= 0) then
      if (nf90_inq_dimid(ncid, 'time', dimid) == NF90_NOERR) recdim = dimid
    end if
    ntime = -1
    if (recdim > 0) then
      ntime = axes(recdim)%len
      ! mpp_io required a coordinate variable for the record dimension
      call nc_check(nf90_inq_varid(ncid, trim(axes(recdim)%name), k), trim(file_name)//' '//axes(recdim)%name)
    end if
    allocate (time_values(max(ntime, 0)))

! A variable whose name matches a dimension name is a coordinate variable.
    allocate (isdim(nvars))
    do i = 1, nvars
      call nc_check(nf90_inquire_variable(ncid, i, name=name), file_name)
      isdim(i) = .false.
      do j = 1, ndims
        if (trim(lowercase(name)) == trim(lowercase(axes(j)%name))) isdim(i) = .true.
      end do
    end do
    allocate (fields(count(.not. isdim)))

    nf = 0
    do i = 1, nvars
      call nc_check(nf90_inquire_variable(ncid, i, name, ndims=nvdims, dimids=dimids, nAtts=nvatts), file_name)
      if (isdim(i)) then
        call nc_check(nf90_inq_dimid(ncid, trim(name), dimid), trim(file_name)//' '//name)
        if (dimid == recdim) then
          call nc_check(nf90_get_var(ncid, i, time_values), trim(file_name)//' '//name)
        else
          call nc_check(nf90_get_var(ncid, i, axes(dimid)%data), trim(file_name)//' '//name)
        end if
        do k = 1, nvatts
          call nc_check(nf90_inq_attname(ncid, i, k, attname), file_name)
          select case (trim(attname))
          case ('units')
            call get_text_att(ncid, i, attname, axes(dimid)%units)
          case ('calendar', 'calendar_type')
            call get_text_att(ncid, i, attname, axes(dimid)%calendar)
            j = index(axes(dimid)%calendar, achar(0))
            if (j > 0) axes(dimid)%calendar(j:j) = ' '
            axes(dimid)%calendar = lowercase(axes(dimid)%calendar)
            select case (trim(axes(dimid)%calendar))
            case ('none')
              axes(dimid)%calendar = 'no_calendar'
            case ('no_leap')
              axes(dimid)%calendar = 'noleap'
            case ('365_days')
              axes(dimid)%calendar = '365_day'
            case ('360_days')
              axes(dimid)%calendar = '360_day'
            end select
          case ('positive')
            positive = ''
            call get_text_att(ncid, i, attname, positive)
            if (positive == 'down') then
              axes(dimid)%sense = -1
            else if (positive == 'up') then
              axes(dimid)%sense = 1
            end if
          end select
        end do
      else
        nf = nf + 1
        fields(nf)%name = name
        fields(nf)%varid = i
        allocate (fields(nf)%count(nvdims))
        do j = 1, nvdims
          if (dimids(j) == recdim) then
            fields(nf)%count(j) = 1
            fields(nf)%tdim = j
          else
            fields(nf)%count(j) = axes(dimids(j))%len
          end if
        end do
        do k = 1, nvatts
          call nc_check(nf90_inq_attname(ncid, i, k, attname), file_name)
          select case (trim(attname))
          case ('units')
            call get_text_att(ncid, i, attname, fields(nf)%units)
          case ('scale_factor')
            call get_real_att(ncid, i, attname, fields(nf)%scale)
          case ('add_offset')
            call get_real_att(ncid, i, attname, fields(nf)%add)
          end select
        end do
      end if
    end do
    deallocate (isdim)

  end subroutine get_file_info
!
!#######################################################################
!
  subroutine get_text_att(ncid, varid, attname, value)
! Set value to a text attribute. Other attribute types leave it unchanged.
    integer, intent(in)    :: ncid, varid
    character(len=*), intent(in)    :: attname
    character(len=*), intent(inout) :: value

    integer :: xtype, attlen
    character(len=:), allocatable :: buf

    call nc_check(nf90_inquire_attribute(ncid, varid, trim(attname), xtype, attlen), attname)
    if (xtype /= NF90_CHAR) return
    allocate (character(len=attlen) :: buf)
    call nc_check(nf90_get_att(ncid, varid, trim(attname), buf), attname)
    value = buf

  end subroutine get_text_att
!
!#######################################################################
!
  subroutine get_real_att(ncid, varid, attname, value)
! Set value to the first element of a numeric attribute. Text attributes leave it unchanged.
    integer, intent(in)    :: ncid, varid
    character(len=*), intent(in)    :: attname
    real, intent(inout) :: value

    integer :: xtype, attlen
    real, allocatable :: buf(:)

    call nc_check(nf90_inquire_attribute(ncid, varid, trim(attname), xtype, attlen), attname)
    if (xtype == NF90_CHAR .or. attlen < 1) return
    allocate (buf(attlen))
    call nc_check(nf90_get_att(ncid, varid, trim(attname), buf), attname)
    value = buf(1)

  end subroutine get_real_att
!
!#######################################################################
!
  subroutine nc_check(status, context)
! FATAL error on a netCDF error.
    integer, intent(in) :: status
    character(len=*), intent(in) :: context

    if (status /= NF90_NOERR) call mpp_error(FATAL, 'mima_interpolator_mod: '//trim(context)// &
                                             ': '//trim(nf90_strerror(status)))

  end subroutine nc_check
!
!#################################################################
!
  !> Interpolates several columns (second index) by pressure-weighted averaging.
  subroutine interp_weighted_scalar_2D(grdin, grdout, datin, datout)
    real, intent(in), dimension(:) :: grdin, grdout
    real, intent(in), dimension(:, :) :: datin
    real, intent(out), dimension(:, :) :: datout

    integer :: j, k, n

    if (size(grdin(:)) .ne. (size(datin, 1) + 1)) &
      call mpp_error(FATAL, 'interp_weighted_scalar : input data and pressure do not have the same number of levels')
    if (size(grdout(:)) .ne. (size(datout, 1) + 1)) &
      call mpp_error(FATAL, 'interp_weighted_scalar : output data and pressure do not have the same number of levels')

    do k = 1, size(datout, 1)
      datout(k, :) = 0.0

      do j = 1, size(datin, 1)

        if (grdin(j) <= grdout(k) .and. &
            grdin(j + 1) >= grdout(k) .and. &
            grdin(j + 1) <= grdout(k + 1)) then

          do n = 1, size(datin, 2)
            datout(k, n) = datout(k, n) + datin(j, n)*(grdin(j + 1) - grdout(k))
          end do

        else if (grdin(j) >= grdout(k) .and. &
                 grdin(j) <= grdout(k + 1) .and. &
                 grdin(j + 1) >= grdout(k + 1)) then

          do n = 1, size(datin, 2)
            datout(k, n) = datout(k, n) + datin(j, n)*(grdout(k + 1) - grdin(j))
          end do

        else if (grdin(j) >= grdout(k) .and. &
                 grdin(j + 1) <= grdout(k + 1)) then

          do n = 1, size(datin, 2)
            datout(k, n) = datout(k, n) + datin(j, n)*(grdin(j + 1) - grdin(j))
          end do

        else if (grdin(j) <= grdout(k) .and. &
                 grdin(j + 1) >= grdout(k + 1)) then

          do n = 1, size(datin, 2)
            datout(k, n) = datout(k, n) + datin(j, n)*(grdout(k + 1) - grdout(k))

          end do
        end if

      end do

      do n = 1, size(datin, 2)
        datout(k, n) = datout(k, n)/(grdout(k + 1) - grdout(k))
      end do

    end do

  end subroutine interp_weighted_scalar_2D

!---------------------------------------------------------------------

  !> Interpolates a column by pressure-weighted averaging of the overlapping layers.
  subroutine interp_weighted_scalar_1D(grdin, grdout, datin, datout)
    real, intent(in), dimension(:) :: grdin, grdout, datin
    real, intent(out), dimension(:) :: datout

    integer :: j, k

    if (size(grdin(:)) .ne. (size(datin(:)) + 1)) &
      call mpp_error(FATAL, 'interp_weighted_scalar : input data and pressure do not have the same number of levels')
    if (size(grdout(:)) .ne. (size(datout(:)) + 1)) &
      call mpp_error(FATAL, 'interp_weighted_scalar : output data and pressure do not have the same number of levels')

    do k = 1, size(datout(:))
      datout(k) = 0.0

      do j = 1, size(datin(:))

        if (grdin(j) <= grdout(k) .and. &
            grdin(j + 1) >= grdout(k) .and. &
            grdin(j + 1) <= grdout(k + 1)) then

          datout(k) = datout(k) + datin(j)*(grdin(j + 1) - grdout(k))

        else if (grdin(j) >= grdout(k) .and. &
                 grdin(j) <= grdout(k + 1) .and. &
                 grdin(j + 1) >= grdout(k + 1)) then

          datout(k) = datout(k) + datin(j)*(grdout(k + 1) - grdin(j))

        else if (grdin(j) >= grdout(k) .and. &
                 grdin(j + 1) <= grdout(k + 1)) then

          datout(k) = datout(k) + datin(j)*(grdin(j + 1) - grdin(j))

        else if (grdin(j) <= grdout(k) .and. &
                 grdin(j + 1) >= grdout(k + 1)) then

          datout(k) = datout(k) + datin(j)*(grdout(k + 1) - grdout(k))

        end if

      end do

      datout(k) = datout(k)/(grdout(k + 1) - grdout(k))

    end do

  end subroutine interp_weighted_scalar_1D
!
!#################################################################
!
  !> Interpolates a column linearly in pressure, extrapolating from the end values.
  subroutine interp_linear(grdin, grdout, datin, datout)
    real, intent(in), dimension(:) :: grdin, grdout, datin
    real, intent(out), dimension(:) :: datout

    integer :: j, k, n
    real    :: wt

    if (size(grdin(:)) .ne. (size(datin(:)) + 1)) &
      call mpp_error(FATAL, 'interp_linear : input data and pressure do not have the same number of levels')
    if (size(grdout(:)) .ne. (size(datout(:)) + 1)) &
      call mpp_error(FATAL, 'interp_linear : output data and pressure do not have the same number of levels')

    n = size(grdin(:))

! datin(j) is taken to be at grdin(j), j = 1..n-1. Search only the brackets
! between those points, so that beyond them the end bracket extrapolates.
    if (n == 2) then
      datout(:) = datin(1)
      return
    end if

    do k = 1, size(datout(:))

      ! ascending grid values
      if (grdin(1) < grdin(n)) then
        do j = 2, n - 2
          if (grdout(k) <= grdin(j)) exit
        end do
        ! descending grid values
      else
        do j = n - 1, 3, -1
          if (grdout(k) <= grdin(j - 1)) exit
        end do
      end if

      ! linear interpolation
      wt = (grdout(k) - grdin(j - 1))/(grdin(j) - grdin(j - 1))
!print '(a,2i3,4f6.1)', 'k,j=',k,j,grdout(k),grdin(j-1),grdin(j),wt
      ! constant value extrapolation
      ! wt = min(max(wt,0.),1.)

      datout(k) = (1.-wt)*datin(j - 1) + wt*datin(j)

    end do

  end subroutine interp_linear
!
!########################################################################
!
end module mima_interpolator_mod
