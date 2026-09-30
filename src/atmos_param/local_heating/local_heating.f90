!> Prescribed Gaussian heating, for up to `ngauss` heat sources.
!>
!> Each source is Gaussian in longitude, latitude, log-pressure and time, and can move
!> in longitude, latitude and pressure. Sources with `pcenter` > 0 heat the atmosphere
!> (`local_heating`, called by the physics driver); sources with `pcenter` < 0 heat the
!> surface (`horizontal_heating`, called by `simple_surface`). Used with
!> `do_local_heating = .true.` in `physics_driver_nml`.
!>
!> Namelist: `local_heating_nml`
!> ([namelist reference](https://eddy-stanford.github.io/MiMA/Parameters/#local_heating_nml)).
module local_heating_mod

  use fms_mod, only: error_mesg, FATAL, &
                     input_nml_file, &
                     check_nml_error, &
                     mpp_pe, mpp_root_pe, &
                     write_version_number, stdlog, &
                     uppercase, &  !pjk
                     mpp_clock_id, mpp_clock_begin, mpp_clock_end, CLOCK_COMPONENT!,mpp_chksum

  use diag_manager_mod, only: register_diag_field, send_data
  use time_manager_mod, only: time_type

  use time_manager_mod, only: time_type, get_time, length_of_year

  use diag_manager_mod, only: register_diag_field, send_data

  use field_manager_mod, only: MODEL_ATMOS, parse

  use constants_mod, only: RADIAN, PI

  implicit none

  !-----------------------------------------------------------------------
  !---------- interfaces ------------
  public :: local_heating, local_heating_init, horizontal_heating

  !---------------------------------------------------------------------------------------------------------------
  !
  !-------------------- diagnostics fields -------------------------------
  integer :: id_tdt_lheat  !! diagnostic id of `tdt_lheat`
  character(len=14) :: mod_name = 'local_heating'  !! module name for the diagnostics
  real :: missing_value = -999.  !! missing value of the diagnostics
  !-----------------------------------------------------------------------
  !-------------------- namelist -----------------------------------------
  !-----------------------------------------------------------------------
  integer, parameter :: ngauss = 10  !! maximum number of heat sources
  real, dimension(ngauss)   :: hamp = 0.        !! [K/day] amplitude of the heating
  real, dimension(ngauss)   :: lonwidth = -1.       !! [deg] zonal width, if `loncenter` >= 0
  real, dimension(ngauss)   :: loncenter = -1.       !! [deg] longitude of the centre; zonally symmetric if < 0
  real, dimension(ngauss)   :: lonmove = 0.        !! [deg/day] zonal speed of the source
  real, dimension(ngauss)   :: latwidth = 15.       !! [deg] meridional width
  real, dimension(ngauss)   :: latcenter = 0.        !! [deg] latitude of the centre
  real, dimension(ngauss)   :: latmove = 0.        !! [deg/day] meridional speed of the source
  real, dimension(ngauss)   :: pwidth = 1.        !! [log10(hPa)] vertical width; constant in the vertical if < 0
  real, dimension(ngauss)   :: pcenter = -1.        !! [hPa] pressure of the centre; surface heating if < 0
  real, dimension(ngauss)   :: pmove = 0.        !! [hPa/day] vertical speed of the source
  logical, dimension(ngauss):: is_periodic = .false.
  !! reset the position periodically (with `tphase` and `tperiod`): periodic in longitude and
  !! pressure, back and forth in latitude
  real, dimension(ngauss)   :: twidth = -1.        !! [days] temporal width; constant in time if < 0
  real, dimension(ngauss)   :: tphase = 0.        !! [days] temporal phase
  real, dimension(ngauss)   :: tperiod = -1.
  !! temporal period: [fraction of a year] if < 0, [days] if > 0

  namelist /local_heating_nml/ hamp &
    , lonwidth, loncenter, lonmove &
    , latwidth, latcenter, latmove &
    , pwidth, pcenter, pmove &
    , is_periodic &
    , twidth, tphase, tperiod

  ! local variables
  real, dimension(ngauss) :: logpc  !! not used
  integer                :: daysperyear  !! number of days in a year
  logical                :: do_3d_heating  !! true if any source has `pcenter` > 0

contains

  !> Initializes the module: reads `local_heating_nml`, converts the namelist values to
  !> SI units and radians, and registers the diagnostic `tdt_lheat` if any source heats
  !> the atmosphere.
  subroutine local_heating_init(axes, Time)
    implicit none
    integer, intent(in), dimension(4) :: axes  !! diagnostic axes (lon, lat, pfull, phalf)
    type(time_type), intent(in)       :: Time  !! current time
    !-----------------------------------------------------------------------
    integer :: seconds
    integer :: unit, io, ierr, n

    !     ----- read namelist -----

    read (input_nml_file, nml=local_heating_nml, iostat=io)
    ierr = check_nml_error(io, 'local_heating_nml')

    ! ---- convert input units to code units  -----
    call get_time(length_of_year(), seconds, daysperyear)
    do_3d_heating = .false.
    do n = 1, ngauss
      pcenter(n) = pcenter(n)*100      ! convert hPa to Pa
      if (pcenter(n) .gt. 0.0) do_3d_heating = .true.
      hamp(n) = hamp(n)/86400.      ! convert K/d to K/s
      loncenter(n) = loncenter(n)/RADIAN ! convert degrees to radians
      lonwidth(n) = lonwidth(n)/RADIAN  ! convert degrees to radians
      lonmove(n) = lonmove(n)/RADIAN/86400. ! convert degrees/day to radians/s
      latcenter(n) = latcenter(n)/RADIAN ! convert degrees to radians
      latwidth(n) = latwidth(n)/RADIAN  ! convert degrees to radians
      latmove(n) = latmove(n)/RADIAN/86400. ! convert degrees/day to radians/s
      pmove(n) = pmove(n)*100./86400.! convert hPa/day to Pa/s
      if (tperiod(n) .lt. 0.0) tperiod(n) = -tperiod(n)*daysperyear ! convert year fraction to day of year
      twidth(n) = twidth(n)*86400     ! convert to seconds
      tphase(n) = tphase(n)*86400     ! convert to seconds
      tperiod(n) = tperiod(n)*86400    ! convert to seconds
    end do
    !----
    !------------ initialize diagnostic fields ---------------
    ! only needed if we actually do 3D heating
    if (do_3d_heating) then
      id_tdt_lheat = &
        register_diag_field(mod_name, 'tdt_lheat', axes(1:3), Time, &
                            'Temperature tendency due to local heating', &
                            'K/s', missing_value=missing_value)
    else
      id_tdt_lheat = 0
    end if

  end subroutine local_heating_init

  !-----------------------------------------------------------------------
  !-------------------- computing localized heating ----------------------
  !-----------------------------------------------------------------------

  !> Adds the heating of the atmosphere by the sources with `pcenter` > 0 to `tdt_tot`.
  subroutine local_heating(is, js, Time, lon, lat, p_full, tdt_tot)
    implicit none
    integer, intent(in)                  :: is, js  !! starting subdomain i, j indices of the physics window
    type(time_type), intent(in)           :: Time  !! current time
    real, dimension(:, :), intent(in)    :: lon, lat  !! longitudes and latitudes [rad]
    real, dimension(:, :, :), intent(in)    :: p_full  !! pressure at full levels [Pa]
    real, dimension(:, :, :), intent(inout) :: tdt_tot  !! temperature tendency, to which the heating is added [K/s]
    ! local variables
    integer :: i, j, k, n
    real, dimension(size(lon, 1), size(lon, 2)) :: horiz_tdt
    real, dimension(size(tdt_tot, 1), size(tdt_tot, 2), size(tdt_tot, 3)) :: tdt
    real    :: logp, p_factor, tcenter(3, ngauss)
    logical :: used

    tdt = 0.
    ! if local heating is 2D only it is done via horizontal_heating
    !  within simple_surface.f90
    if (do_3d_heating) then
      ! horizontal heating first
      call horizontal_heating(Time, lon, lat, horiz_tdt, tcenter)
      ! then vertical heating
      do n = 1, ngauss
        if (hamp(n) .ne. 0. .and. pcenter(n) .gt. 0.) then
          ! add vertical component
          do k = 1, size(p_full, 3)
            do j = 1, size(lon, 2)
              do i = 1, size(lon, 1)
                ! vertical component
                if (pwidth(n) .lt. 0.0) then
                  p_factor = 1.0
                else
                  logp = log10(p_full(i, j, k))
                  p_factor = exp(-(logp - tcenter(3, n))**2/(2*(pwidth(n))**2))
                end if
                ! everything together
                tdt(i, j, k) = tdt(i, j, k) + horiz_tdt(i, j)*p_factor
              end do
            end do
          end do
        end if
      end do
    end if

    tdt_tot = tdt_tot + tdt

    !------- diagnostics ------------
    if (id_tdt_lheat > 0) then
      used = send_data(id_tdt_lheat, tdt, Time, is, js, 1)
    end if

  end subroutine local_heating

  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------

  !> Returns the horizontal and temporal part of the heating, summed over the sources.
  !>
  !> Without `tcenter` (the call from `simple_surface`) it sums the surface sources
  !> (`pcenter` < 0); with `tcenter` it sums the atmospheric sources (`pcenter` > 0) and
  !> returns their current centres.
  subroutine horizontal_heating(Time, lon, lat, horiz_tdt, tcenter)
    implicit none
    type(time_type), intent(in)       :: Time  !! current time
    real, dimension(:, :), intent(in)  :: lon, lat  !! longitudes and latitudes [rad]
    real, dimension(:, :), intent(out) :: horiz_tdt  !! heating rate [K/s]
    real, dimension(3, ngauss), intent(out), optional :: tcenter
    !! current centre of each source: longitude [rad], latitude [rad] and log10 of the
    !! pressure [log10(Pa)] (zero for sources that are not used)
    ! local variables
    integer :: i, j, d, n, deltasecs
    integer :: seconds, days, fullseconds
    real    :: tcent(3), t_factor, targ, halfper
    real, dimension(size(lon, 1), size(lon, 2)) :: lon_factor, lat_factor
    logical :: do_horiz

    call get_time(Time, seconds, days)
    fullseconds = days*86400 + seconds

    horiz_tdt = 0.0
    if (present(tcenter)) then
      tcenter = 0.0
    end if
    do n = 1, ngauss
      ! check if heating should be added
      do_horiz = .false.
      !  first, is the amplitude non-zero?
      if (hamp(n) .ne. 0.) then
        !  then, 3D vs. 2D heating
        !  If pcenter < 0: 2D heating -> .not. present(tcenter)
        !  If pcenter > 0: 3D heating -> present(tcenter)
        if (present(tcenter)) then ! calling from 3D atmosphere
          if (pcenter(n) .gt. 0.0) then
            do_horiz = .true.
          end if
        elseif (pcenter(n) .lt. 0.0) then ! calling from simple_surface
          do_horiz = .true.
        end if
      end if
      ! compute horizontal heating if it should be added
      if (do_horiz) then
        ! local heating position is determined at peak heating time
        halfper = 0.5*tperiod(n)
        if (is_periodic(n)) then
          deltasecs = modulo(fullseconds - tphase(n) + halfper, tperiod(n)) - halfper
        else
          deltasecs = fullseconds
        end if
        ! compute the center position
        tcent(1) = mod(loncenter(n) + lonmove(n)*deltasecs, 2*PI)
        tcent(2) = latcenter(n) + latmove(n)*abs(deltasecs)
        if (pcenter(n) .gt. 0.0) then
          tcent(3) = log10(pcenter(n) + pmove(n)*abs(deltasecs))
        else
          tcent(3) = pcenter(n)
        end if
        ! temporal component
        if (twidth(n) .lt. 0.0) then
          t_factor = 1.0
        else
          targ = mod(fullseconds - tphase(n) + halfper, tperiod(n)) - halfper
          t_factor = exp(-(targ)**2/(2*(twidth(n))**2))
        end if
        ! meridional and zonal components
        do j = 1, size(lon, 2)
          do i = 1, size(lon, 1)
            if (loncenter(n) .ge. 0.0) then
              lon_factor(i, j) = exp(-(lon(i, j) - tcent(1))**2/(2*(lonwidth(n))**2))
              ! there is a problem when the heating is close to 360/0
              lon_factor(i, j) = max(lon_factor(i, j), &
                                     exp(-(lon(i, j) + 2*PI - tcent(1))**2/(2*(lonwidth(n))**2)))
              lon_factor(i, j) = max(lon_factor(i, j), &
                                     exp(-(lon(i, j) - 2*PI - tcent(1))**2/(2*(lonwidth(n))**2)))
            else
              lon_factor(i, j) = 1.0
            end if
            lat_factor(i, j) = exp(-(lat(i, j) - tcent(2))**2/(2*(latwidth(n))**2))
          end do
        end do
        horiz_tdt = horiz_tdt + hamp(n)*t_factor*lon_factor*lat_factor
        if (present(tcenter)) then
          do d = 1, 3
            tcenter(d, n) = tcent(d)
          end do
        end if
      end if
    end do

  end subroutine horizontal_heating

end module local_heating_mod
