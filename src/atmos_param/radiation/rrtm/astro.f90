!> Orbit, solar constant and solar zenith angle for the RRTM radiation.
!>
!> Holds the astronomical parameters used by `rrtm_radiation` and computes the cosine of the
!> solar zenith angle, instantaneous, averaged over an interval, or as a daily mean. The
!> orbit is circular; the declination follows from `obliq` and the day of the year relative
!> to the March equinox (`equinox_day`).
!>
!> Namelist: `astro_nml`
!> ([namelist reference](https://eddy-stanford.github.io/MiMA/Parameters/#astro_nml)).
!>
!> Original authors: Martin Jucker.
module rrtm_astro
!
!   Martin Jucker, 2015, https://github.com/mjucker/MiMA.
!
! Modules
  use parkind, only: im => kind_im, rb => kind_rb
  use fms_mod, only: input_nml_file, check_nml_error, &
                     error_mesg, FATAL
! Variables
  implicit none
  logical          :: astro_initialized = .false.  !! whether `astro_init` has been called
!
!---------------------------------------------------------------------------------------------------------------
!                                namelist values
!---------------------------------------------------------------------------------------------------------------
  real(kind=rb)      :: obliq = 23.439             !! [deg] obliquity
  logical            :: use_dyofyr = .false.            !! let RRTM compute the Earth-Sun distance from the day of the year
                                                        !! (assumes 365 days per year)
  real(kind=rb)      :: solr_cnst = 1370.              !! [W/m2] solar constant
  real(kind=rb)      :: solrad = 1.0                      !! Earth-Sun distance factor if `use_dyofyr = .false.`
  integer(kind=im)   :: solday = 0                        !! if > 0, perpetual run at this day of the year
  real(kind=rb)      :: equinox_day = 0.25                !! fraction of the year at which the March equinox occurs

  namelist /astro_nml/ obliq, use_dyofyr, solr_cnst, solrad, solday, equinox_day

contains
!--------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------
  !> Initializes the module: reads `astro_nml`.
  subroutine astro_init
    implicit none
    integer :: unit, ierr, io

    read (input_nml_file, nml=astro_nml, iostat=io)
    ierr = check_nml_error(io, 'astro_nml')

    astro_initialized = .true.

  end subroutine astro_init
!--------------------------------------------------------------------------------------
  !> Computes the cosine of the solar zenith angle for the RRTM shortwave radiation.
  !>
  !> If `0 < dt < 86400` the value is averaged over the interval from `Time` to `Time + dt`
  !> (zero at night); if `dt >= 86400` it is the daily mean; otherwise it is the
  !> instantaneous value. Parts of this are taken from GFDL's `astronomy.f90`.
  subroutine compute_zenith(Time, equinox_day, dt, lat, lon, cosz, dyofyr)
!
! Modules
    use time_manager_mod, only: time_type, get_time, length_of_year
    use constants_mod, only: PI
    use fms_mod, only: error_mesg, FATAL
! Local variables
    implicit none
! Inputs
    type(time_type), intent(in) :: Time        !! time of year, according to calendar
    real(kind=rb), intent(in) :: equinox_day !! fraction of the year at which the March equinox occurs
    integer(kind=im), intent(in) :: dt          !! averaging interval (if > 0) [s]
    real(kind=rb), dimension(:, :), intent(in) :: lat, lon     !! latitudes and longitudes [rad]
    real(kind=rb), dimension(:, :), intent(out):: cosz        !! cosine of the zenith angle
    integer(kind=im), intent(out):: dyofyr      !! day of the year, counted from the March equinox, at which
                                                !! `cosz` is computed
! Locals
    real(kind=rb), dimension(size(lat, 1), size(lat, 2)) :: h, cos_h, &
                                                            lat_h

    real     :: dec_sin, dec_tan, dec, dec_cos, twopi, dt_pi

    integer  :: seconds, sec2, days, daysperyear
    real, dimension(size(lon, 1), size(lon, 2)) :: time_pi, aa, bb, tt, st, stt, sh, fracday
    real     :: radsec, radday

    integer  :: i, j
! Constants
    real     :: radpersec, radperday
    real     :: eps = 1.0e-05, deg2rad
!--------------------------------------------------------------------------------------
    deg2rad = PI/180.
    twopi = 2*PI

    if (.not. astro_initialized) then
      call error_mesg('astro', 'astro_mod not initialized', FATAL)
    end if

    call get_time(length_of_year(), sec2, daysperyear)
    if (daysperyear .ne. 365 .and. use_dyofyr) then
      print *, ' number of days per year: ', daysperyear
      call error_mesg('astro', &
                      ' use_dyofyr is TRUE but the calendar year does not have 365 days. STOPPING', &
                      FATAL)
    end if

    radpersec = 2*PI/86400.
    radperday = 2*PI/daysperyear

    !get the time for origin
    call get_time(Time, seconds, days)
    !convert into radians
    radsec = seconds*radpersec
    dt_pi = dt*radpersec

    !set local time throughout the globe
    do i = 1, size(lon, 1)
      !move it into interval [-PI,PI]
      time_pi(i, :) = modulo(radsec + lon(i, :), 2*PI) - PI
    end do
    where (time_pi >= PI) time_pi = time_pi - twopi
    where (time_pi < -PI) time_pi = time_pi + twopi
    !time_pi now contains local time at each grid point
    !get day of the year relative to March equinox. We set equinox at (equinox_day,equinox_day+0.5)*daysperyear
    !note that GFDL computes relative to September equinox
    days = days - int(equinox_day*daysperyear)
    dyofyr = modulo(days, daysperyear)
    !convert into radians
    radday = dyofyr*radperday
    !get declination
    dec_sin = sin(obliq*deg2rad)*sin(radday) !("-" in GFDL's code due to differences in origin (March vs. September equinox))
    dec = asin(dec_sin)
    dec_cos = cos(dec)

    !now compute the half day, to determine if it's day or night
    dec_tan = tan(dec)
    lat_h = lat
    where (lat_h == 0.5*PI) lat_h = lat - eps
    where (lat_h == -0.5*PI) lat_h = lat + eps
    cos_h = -tan(lat_h)*dec_tan
    where (cos_h <= -1.0) h = PI
    where (cos_h >= 1.0) h = 0.0
    where (cos_h > -1.0 .and. cos_h < 1.0) &
      h = acos(cos_h)
!---------------------------------------------------------------------
!    define terms needed in the cosine zenith angle equation.
!--------------------------------------------------------------------
    aa = sin(lat)*dec_sin
    bb = cos(lat)*dec_cos

    !finally, compute the zenith angle
    if (dt > 0 .and. dt < 86400.) then !average over some given time interal dt
      tt = time_pi + dt_pi
      st = sin(time_pi)
      stt = sin(tt)
      sh = sin(h)
      cosz = 0.0
!-------------------------------------------------------------------
!    case 1: entire averaging period is before sunrise.
!-------------------------------------------------------------------
      where (time_pi < -h .and. tt < -h) cosz = 0.0

!-------------------------------------------------------------------
!    case 2: averaging period begins before sunrise, ends after sunrise
!    but before sunset
!-------------------------------------------------------------------
      where ((tt + h) /= 0.0 .and. time_pi < -h .and. abs(tt) <= h) &
        cosz = aa + bb*(stt + sh)/(tt + h)
!-------------------------------------------------------------------
!    case 3: averaging period begins before sunrise, ends after sunset,
!    but before the next sunrise. modify if averaging period extends
!    past the next day's sunrise, but if averaging period is less than
!    a half- day (pi) that circumstance will never occur.
!-------------------------------------------------------------------
      where (time_pi < -h .and. h /= 0.0 .and. h < tt) &
        cosz = aa + bb*(sh + sh)/(h + h)
!-------------------------------------------------------------------
!    case 4: averaging period begins after sunrise, ends before sunset.
!-------------------------------------------------------------------
      where (abs(time_pi) <= h .and. abs(tt) <= h) &
        cosz = aa + bb*(stt - st)/(tt - time_pi)
!-------------------------------------------------------------------
!    case 5: averaging period begins after sunrise, ends after sunset.
!    modify when averaging period extends past the next day's sunrise.
!-------------------------------------------------------------------
      where ((h - time_pi) /= 0.0 .and. abs(time_pi) <= h .and. h < tt) &
        cosz = aa + bb*(sh - st)/(h - time_pi)
!-------------------------------------------------------------------
!    case 6: averaging period begins after sunrise , ends after the
!    next day's sunrise. note that this includes the case when the
!    day length is one day (h = pi).
!-------------------------------------------------------------------
      where (twopi - h < tt .and. (tt + h - twopi) /= 0.0 .and. time_pi <= h) &
        cosz = (cosz*(h - time_pi) + (aa*(tt + h - twopi) + &
                                      bb*(stt + sh)))/((h - time_pi) + (tt + h - twopi))

!-------------------------------------------------------------------
!    case 7: averaging period begins after sunset and ends before the
!    next day's sunrise
!-------------------------------------------------------------------
      where (h < time_pi .and. twopi - h >= tt) cosz = 0.0

!-------------------------------------------------------------------
!    case 8: averaging period begins after sunset and ends after the
!    next day's sunrise but before the next day's sunset. if the
!    averaging period is less than a half-day (pi) the latter
!    circumstance will never occur.
!-----------------------------------------------------------------
      where (h < time_pi .and. twopi - h < tt) &
        cosz = aa + bb*(stt + sh)/(tt + h - twopi)

!-------------------------------------------------------------------
!    day fraction is the fraction of the averaging period contained
!    within the (-h,h) period.
!-------------------------------------------------------------------
      where (time_pi < -h .and. tt < -h) fracday = 0.0
      where (time_pi < -h .and. abs(tt) <= h) fracday = (tt + h)/dt
      where (time_pi < -h .and. h < tt) fracday = (h + h)/dt
      where (abs(time_pi) <= h .and. abs(tt) <= h) fracday = (tt - time_pi)/dt
      where (abs(time_pi) <= h .and. h < tt) fracday = (h - time_pi)/dt
      where (h < time_pi) fracday = 0.0
      where (twopi - h < tt) fracday = fracday + (tt + h - twopi)/dt
!-------------------------------------------------------------------
!    now we need to correct cosz by the fraction of day when
!     averaging
!-------------------------------------------------------------------
      cosz = cosz*fracday/radpersec
!-------------------------------------------------------------------
!    daily mean
!-------------------------------------------------------------------
    else if (dt .ge. 86400.) then
      cosz = (aa*h + bb*sin(h))/PI

!----------------------------------------------------------------------
!    if instantaneous values are desired, define cosz at time t.
!----------------------------------------------------------------------
    else !version w/o time averaging
      where (abs(time_pi) <= h)
        cosz = aa + bb*cos(time_pi)
      elsewhere
        cosz = 0.0
      end where
    end if
!----------------------------------------------------------------------
!    be sure that cosz is not negative.
!----------------------------------------------------------------------
    cosz = max(0.0, cosz)

  end subroutine compute_zenith

end module rrtm_astro
