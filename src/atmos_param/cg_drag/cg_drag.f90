!> Non-orographic gravity-wave drag, after Alexander and Dunkerton (1999).
!>
!> Computes the forcing of the zonal and meridional winds by a spectrum of
!> gravity waves launched at a source level in the troposphere, with the
!> changes of Cohen et al. (2013). The forcing is recalculated every
!> `cg_drag_freq` seconds and held fixed in between; it is saved in the restart
!> file `RESTART/cg_drag.res.nc`.
!>
!> In each column, and for each wind component separately (`gwfc`):
!>
!> * The source spectrum (AD99 eq. 17) at the phase speeds c = -`cmax`,
!>   -`cmax` + `dc`, ..., `cmax` is a Gaussian of amplitude `Bw` and half-width
!>   `cw`, centred on the wind u0 at the source level (or on c = 0), with the
!>   sign of c - u0.
!> * The intermittency (eq. 19) scales the spectrum so that the sum of |B0| over
!>   all phase speeds carries the source stress `Bt`.
!> * Going up from the source, a wave is removed where it is reflected (eq. 21)
!>   and deposits its momentum flux where it breaks (eq. 24) or meets a critical
!>   level. The waves that are already unstable at the source level are removed
!>   without depositing (eq. 22): with the default spectrum about 20% of `Bt`.
!> * The momentum flux that reaches the model top is deposited between the top
!>   and `damp_level_pressure`, so none leaves the model.
!> * The flux deposited in the layer between two full levels is shared between
!>   the two model levels by the pressure thickness of the layer on each side of
!>   the half level, and converted to an acceleration with the mass of each
!>   model level, so the column receives exactly the deposited momentum.
!>
!> The zonal and the meridional forcing are computed independently, each with
!> the whole of `Bt` (the spectrum centred on u0 and on v0): a column launches
!> `Bt` in the x and `Bt` in the y direction, 2 `Bt` in total over the four
!> directions.
!>
!> Namelist: `cg_drag_nml`
!> ([namelist reference](https://eddy-stanford.github.io/MiMA/Parameters/#cg_drag_nml)).
!>
!> References:
!>
!> * Alexander, M. J., and T. J. Dunkerton, 1999: A spectral parameterization
!>   of mean-flow forcing due to breaking gravity waves. J. Atmos. Sci., 56,
!>   4167-4180.
!> * Cohen, N. Y., E. P. Gerber, and O. Buhler, 2013: Compensation between
!>   resolved and unresolved wave driving in the stratosphere: implications for
!>   downward control. J. Atmos. Sci., 70, 3780-3798,
!>   https://doi.org/10.1175/JAS-D-12-0346.1.
module cg_drag_mod

  use fms_mod, only: fms_init, mpp_pe, mpp_root_pe, &
                     check_nml_error, &
                     error_mesg, FATAL, NOTE, &
                     input_nml_file, &
                     stdlog
  use time_manager_mod, only: time_manager_init, time_type, get_time, &
                              operator(-)
  use mpp_domains_mod, only: domain2d
  use restart_file_mod, only: restart_file_type, open_restart_read, &
                              open_restart_write, close_restart, &
                              read_restart_field, write_restart_field
  use diag_manager_mod, only: diag_manager_init, &
                              register_diag_field, send_data
  use constants_mod, only: constants_init, PI, RDGAS, GRAV, CP_AIR
  use namelist_check_mod, only: nml_sets

  implicit none
  private

  public cg_drag_init, cg_drag_calc, cg_drag_end, &
    cg_drag_time_vary, cg_drag_endts

!--------------------------------------------------------------------
!---- namelist -----

  integer     :: cg_drag_freq = 21600     !! [s] interval between drag calculations; the drag is held
                                          !! fixed in between (and saved in `cg_drag.res.nc`)
  integer     :: cg_drag_offset = 0   !! [s] on a cold start, time to the first calculation
                                      !! (0: `cg_drag_freq`)

  real        :: source_level_pressure = 315.e+02
  !! [Pa] the source level is the highest level with a reference pressure greater than this at
  !! the equator; towards the poles the number of levels below it shrinks with cos(latitude),
  !! down to level kmax - 1
  real       ::  damp_level_pressure = 0.85e+02
  !! [Pa] momentum flux reaching the model top is deposited from the top down to the lowest level
  !! with a reference pressure below this (the top level if there is none)
  real        :: Bt = 0.0043   !! [Pa] source momentum flux (stress), launched in each of the
                               !! x and the y direction (see the module description)
  real       :: Bw = 0.4  !! [m2/s2] amplitude of the source spectrum. The intermittency scales
                          !! the deposited flux to `Bt`, so `Bw` acts only through the breaking
                          !! criterion: a larger `Bw` makes the waves break lower.
  real       :: cw = 35.0 !! [m/s] half-width at half-maximum of the source spectrum
  real        :: cmax = 99.6          !! [m/s] maximum phase speed; 2 `cmax`/`dc` must be an integer
  real        :: dc = 1.2             !! [m/s] phase-speed resolution
  real        :: wavelength = 300.e3  !! [m] horizontal wavelength of the waves
  logical     :: centre_on_source_wind = .true.
  !! `.true.`: the source spectrum is centred on the wind at the source level (c - u0 = 0);
  !! `.false.`: on c = 0

  namelist /cg_drag_nml/ &
    cg_drag_freq, cg_drag_offset, &
    source_level_pressure, damp_level_pressure, &
    Bt, Bw, cw, cmax, dc, wavelength, centre_on_source_wind

!--------------------------------------------------------------------
!------ private data ------

  integer, allocatable :: source_level(:)  ! k index of the source level at each latitude
  integer :: damp_level  ! k index of the lowest level of the damping layer at the top

  ! gravity-wave forcing of u and v [m/s2], kept between calculation steps (and in the
  ! restart file RESTART/cg_drag.res.nc, with cgdrag_alarm)
  real, allocatable :: gwd_u(:, :, :), gwd_v(:, :, :)

  integer :: nc  ! number of phase speeds in the spectrum (symmetric around c = 0)
  real, allocatable :: c0(:)  ! phase speeds [m/s]
  real :: kwv, k2  ! horizontal wavenumber [1/m] and its square [1/m2]

  real, parameter :: bflim = 2.5e-5  ! minimum squared buoyancy frequency [1/s2]
  real, parameter :: lapse_dry = GRAV/CP_AIR  ! dry adiabatic lapse rate [K/m]

  integer          :: cgdrag_alarm       ! time remaining until the next calculation [s]
  type(time_type)  :: Time_last_call     ! model time of the previous call
  type(domain2d)   :: domain             ! grid domain, for the restart file

  integer          :: id_kedx_cgwd, id_kedy_cgwd, id_bf_cgwd, &
                      id_gwfx_cgwd, id_gwfy_cgwd
  real             :: missing_value = -999.
  character(len=7) :: mod_name = 'cg_drag'

  logical          :: module_is_initialized = .false.

contains

  !> Initializes the module: reads `cg_drag_nml`, finds the source and damping
  !> levels, sets up the phase-speed spectrum, registers the diagnostics and reads
  !> `INPUT/cg_drag.res.nc` if it exists.
  subroutine cg_drag_init(lonb, latb, domain_in, pref, Time, axes)

    real, dimension(:), intent(in)      :: lonb, latb, pref
    !! `lonb`, `latb`: longitudes and latitudes of the cell corners [rad]; `pref`: reference
    !! pressures at full levels, plus the surface value at nlev+1, for a surface pressure of
    !! 1013.25 hPa [Pa]
    type(domain2d), intent(in)      :: domain_in  !! domain decomposition of the model grid
    integer, dimension(4), intent(in)      :: axes  !! diagnostic axes (lon, lat, pfull, phalf)
    type(time_type), intent(in)      :: Time  !! current time

    integer                 :: ierr, io, logunit
    integer                 :: n, j, k
    integer                 :: idf, jdf, kmax, klevel_of_source
    real                    :: alarm, lat, nc_real
    type(restart_file_type) :: rst

    if (module_is_initialized) return

    call fms_init
    call time_manager_init
    call diag_manager_init
    call constants_init

    call check_removed_keys()
    read (input_nml_file, nml=cg_drag_nml, iostat=io)
    ierr = check_nml_error(io, 'cg_drag_nml')
    logunit = stdlog()
    if (mpp_pe() == mpp_root_pe()) write (logunit, nml=cg_drag_nml)

    kmax = size(pref(:)) - 1
    jdf = size(latb) - 1
    idf = size(lonb) - 1

!--------------------------------------------------------------------
!    the source level at the equator is the highest level with a
!    reference pressure greater than source_level_pressure. the damping
!    level is the lowest level above damp_level_pressure, or the top
!    level if there is none.
!--------------------------------------------------------------------
    damp_level = 1
    klevel_of_source = 0
    do k = 1, kmax
      if (pref(k) < damp_level_pressure) then
        damp_level = k
      end if
      if (pref(k) > source_level_pressure) then
        klevel_of_source = k
        exit
      end if
    end do
    if (klevel_of_source == 0) then
      call error_mesg('cg_drag_init', 'no model level has a reference pressure greater than '// &
                      'source_level_pressure', FATAL)
    end if

!    towards the poles the number of levels below the source level
!    shrinks with cos(latitude).
    allocate (source_level(jdf))
    do j = 1, jdf
      lat = 0.5*(latb(j + 1) + latb(j))
      source_level(j) = (kmax + 1) - nint((kmax + 1 - klevel_of_source)*cos(lat))
    end do
    source_level = min(source_level, kmax - 1)
    damp_level = min(damp_level, kmax)
!    the damping layer, where the flux reaching the top is deposited, must lie above the
!    source level of every latitude.
    if (damp_level > minval(source_level)) then
      call error_mesg('cg_drag_init', 'damp_level_pressure must not be below the source '// &
                      'level (source_level_pressure)', FATAL)
    end if

!---------------------------------------------------------------------
!    the phase speeds of the spectrum, symmetric around c = 0.
!---------------------------------------------------------------------
    nc_real = 2.0*cmax/dc
    if (abs(nc_real - nint(nc_real)) > 1.e-6*nc_real) then
      call error_mesg('cg_drag_init', '2*cmax/dc must be an integer', FATAL)
    end if
    nc = nint(nc_real) + 1
    allocate (c0(nc))
    do n = 1, nc
      c0(n) = (n - 1)*dc - cmax
    end do

    kwv = 2.*PI/wavelength
    k2 = kwv*kwv

    id_bf_cgwd = &
      register_diag_field(mod_name, 'bf_cgwd', axes(1:3), Time, &
                          'buoyancy frequency from cg_drag', '1/s', &
                          missing_value=missing_value)
    id_gwfx_cgwd = &
      register_diag_field(mod_name, 'gwfu_cgwd', axes(1:3), Time, &
                          'gravity wave forcing on mean zonal flow', &
                          'm/s2', missing_value=missing_value)
    id_gwfy_cgwd = &
      register_diag_field(mod_name, 'gwfv_cgwd', axes(1:3), Time, &
                          'gravity wave forcing on mean meridional flow', &
                          'm/s2', missing_value=missing_value)
    id_kedx_cgwd = &
      register_diag_field(mod_name, 'kedx_cgwd', axes(1:3), Time, &
                          'effective eddy viscosity from cg_drag (zonal)', 'm2/s', &
                          missing_value=missing_value)
    id_kedy_cgwd = &
      register_diag_field(mod_name, 'kedy_cgwd', axes(1:3), Time, &
                          'effective eddy viscosity from cg_drag (meridional)', 'm2/s', &
                          missing_value=missing_value)

    allocate (gwd_u(idf, jdf, kmax))
    allocate (gwd_v(idf, jdf, kmax))

!--------------------------------------------------------------------
!    if present, read the restart data file. otherwise initialize the
!    gwd fields to zero and define the time remaining until the next
!    cg_drag calculation from the namelist inputs.
!---------------------------------------------------------------------
    domain = domain_in
    if (open_restart_read(rst, 'INPUT/cg_drag.res.nc', domain)) then
      if (mpp_pe() == mpp_root_pe()) call error_mesg('cg_drag_mod', &
                                                     'Reading NetCDF formatted restart file: INPUT/cg_drag.res.nc', NOTE)
      call read_restart_field(rst, 'gwd_u', gwd_u)
      call read_restart_field(rst, 'gwd_v', gwd_v)
      call read_restart_field(rst, 'cgdrag_alarm', alarm)
      cgdrag_alarm = nint(alarm)
      call close_restart(rst)
    else
      gwd_u(:, :, :) = 0.0
      gwd_v(:, :, :) = 0.0
      if (cg_drag_offset > 0) then
        cgdrag_alarm = cg_drag_offset
      else
        cgdrag_alarm = cg_drag_freq
      end if
    end if
!---------------------------------------------------------------------
!    cg_drag_calc is passed the time at the end of each step, so the
!    first call is one model step after Time (cold start or restart).
!---------------------------------------------------------------------
    Time_last_call = Time

    module_is_initialized = .true.

  end subroutine cg_drag_init

  ! Stops with a migration message if input.nml still sets a cg_drag_nml key that was
  ! removed in v2.0 (the namelist read would otherwise only report an unknown variable).
  ! A mention of the key in a comment is ignored.
  subroutine check_removed_keys()

    if (nml_sets('bt_0') .or. nml_sets('bt_eq') .or. nml_sets('bt_nh') .or. nml_sets('bt_sh') .or. &
        nml_sets('phi0n') .or. nml_sets('phi0s') .or. nml_sets('dphin') .or. nml_sets('dphis') .or. &
        nml_sets('cwtropics') .or. nml_sets('kelvin_kludge') .or. nml_sets('bn') .or. &
        nml_sets('cn') .or. nml_sets('nk') .or. nml_sets('flag')) &
      call error_mesg('cg_drag_init', 'cg_drag_nml: Bt_0 is now Bt, nk is replaced by wavelength '// &
                      '(nk = 1: 300.e3) and flag by centre_on_source_wind (flag = 0: .true.); '// &
                      'Bt_eq, Bt_nh, Bt_sh, phi0n, phi0s, dphin, dphis, cwtropics, Bn, cn and '// &
                      'kelvin_kludge were removed (the source spectrum is the same at all '// &
                      'latitudes); see Migration_v2.md', FATAL)

  end subroutine check_removed_keys

  !> Counts down the time to the next drag calculation by the model time elapsed
  !> since the previous call.
  subroutine cg_drag_time_vary(Time)

    type(time_type), intent(in)      :: Time  !! current time

    integer :: sec, day, dt_step

!---------------------------------------------------------------------
!    decrement the time remaining until the next cg_drag calculation by
!    the model time elapsed since the previous call (or since the start
!    of the run). the physics time step is not used, as it is twice the
!    model step on leapfrog steps.
!---------------------------------------------------------------------
    call get_time(Time - Time_last_call, sec, day)
    dt_step = sec + day*86400
    Time_last_call = Time
    cgdrag_alarm = cgdrag_alarm - dt_step

  end subroutine cg_drag_time_vary

  !> Resets the countdown to `cg_drag_freq` after a calculation step.
  subroutine cg_drag_endts

    if (cgdrag_alarm <= 0) then
      cgdrag_alarm = cgdrag_alarm + cg_drag_freq
    end if

  end subroutine cg_drag_endts

  !> Returns the gravity-wave forcing of u and v.
  !>
  !> On a calculation step (every `cg_drag_freq` seconds) it computes the forcing
  !> from the buoyancy frequency, density and winds, stores it for the steps in
  !> between and sends the diagnostics; on the other steps it returns the stored
  !> forcing.
  subroutine cg_drag_calc(is, js, pfull, phalf, zfull, temp, uuu, vvv, &
                          Time, gwfcng_x, gwfcng_y)

    integer, intent(in)      :: is, js  !! indices of the first point of the physics window in the processor domain
    real, dimension(:, :, :), intent(in)      :: pfull, phalf, zfull, temp, uuu, vvv
    !! `pfull`, `phalf`: pressure at full and half levels [Pa]; `zfull`: height at full levels
    !! [m]; `temp`: temperature [K]; `uuu`, `vvv`: zonal and meridional wind [m/s]
    type(time_type), intent(in)      :: Time  !! model time at the end of the step
    real, dimension(:, :, :), intent(out)     :: gwfcng_x, gwfcng_y
    !! tendencies of u and v due to the gravity-wave forcing [m/s2]

    ! column profiles, and the wind-independent quantities of the column (column_setup)
    real, dimension(size(uuu, 3)) :: z, rho, bf, fac, omc
    real, dimension(size(uuu, 1), size(uuu, 2), size(uuu, 3)) :: ked_x, ked_y, bf_out

    real    :: dtdz
    integer :: iz0, imax, jmax, ie, je, i, j, k
    logical :: used

    imax = size(uuu, 1)
    jmax = size(uuu, 2)
    ie = is + imax - 1
    je = js + jmax - 1

    if (cgdrag_alarm <= 0) then

      bf_out = 0.0
      do j = 1, jmax
        ! the source level of this latitude (j is offset in case a physics window is used)
        iz0 = source_level(j + js - 1)
        do i = 1, imax

!----------------------------------------------------------------------
!    buoyancy frequency, from the lapse rate (one-sided at the top,
!    centred below; it is not needed at the lowest level, as the source
!    level is above it), with the squared buoyancy frequency no smaller
!    than bflim; density and height up to one level below the source.
!----------------------------------------------------------------------
          do k = 1, iz0
            if (k == 1) then
              dtdz = (temp(i, j, 1) - temp(i, j, 2))/(zfull(i, j, 1) - zfull(i, j, 2))
            else
              dtdz = (temp(i, j, k - 1) - temp(i, j, k + 1))/(zfull(i, j, k - 1) - zfull(i, j, k + 1))
            end if
            bf(k) = sqrt(max(GRAV*(dtdz + lapse_dry)/temp(i, j, k), bflim))
          end do
          do k = 1, iz0 + 1
            rho(k) = pfull(i, j, k)/(temp(i, j, k)*RDGAS)
            z(k) = zfull(i, j, k)
          end do

          call column_setup(iz0, rho, bf, z, fac, omc)
          call gwfc(iz0, pfull(i, j, :), phalf(i, j, :), bf, fac, omc, uuu(i, j, :), &
                    gwfcng_x(i, j, :), ked_x(i, j, :))
          call gwfc(iz0, pfull(i, j, :), phalf(i, j, :), bf, fac, omc, vvv(i, j, :), &
                    gwfcng_y(i, j, :), ked_y(i, j, :))
          bf_out(i, j, 1:iz0) = bf(1:iz0)
        end do
      end do

!--------------------------------------------------------------------
!    store the forcing for the steps until the next calculation, and
!    send the diagnostics.
!-------------------------------------------------------------------
      gwd_u(is:ie, js:je, :) = gwfcng_x(:, :, :)
      gwd_v(is:ie, js:je, :) = gwfcng_y(:, :, :)

      if (id_kedx_cgwd > 0) then
        used = send_data(id_kedx_cgwd, ked_x, Time, is, js, 1)
      end if
      if (id_kedy_cgwd > 0) then
        used = send_data(id_kedy_cgwd, ked_y, Time, is, js, 1)
      end if
      if (id_bf_cgwd > 0) then
        used = send_data(id_bf_cgwd, bf_out, Time, is, js)
      end if
      if (id_gwfx_cgwd > 0) then
        used = send_data(id_gwfx_cgwd, gwfcng_x, Time, is, js, 1)
      end if
      if (id_gwfy_cgwd > 0) then
        used = send_data(id_gwfy_cgwd, gwfcng_y, Time, is, js, 1)
      end if

    else
      gwfcng_x(:, :, :) = gwd_u(is:ie, js:je, :)
      gwfcng_y(:, :, :) = gwd_v(is:ie, js:je, :)
    end if

!--------------------------------------------------------------------
!    update the alarm clock, which controls how often cg_drag recalculates
!    the forcing.
!--------------------------------------------------------------------
    call cg_drag_endts
    call cg_drag_time_vary(Time)

  end subroutine cg_drag_calc

  !> Writes the restart file `RESTART/cg_drag.res.nc`.
  subroutine cg_drag_end

    type(restart_file_type) :: rst

    if (.not. module_is_initialized) return
    if (mpp_pe() == mpp_root_pe()) call error_mesg('cg_drag_mod', &
                                                   'Writing NetCDF formatted restart file: RESTART/cg_drag.res.nc', NOTE)
    call open_restart_write(rst, 'RESTART/cg_drag.res.nc', domain)
    call write_restart_field(rst, 'gwd_u', gwd_u)
    call write_restart_field(rst, 'gwd_v', gwd_v)
    call write_restart_field(rst, 'cgdrag_alarm', real(cgdrag_alarm))
    call close_restart(rst)

    module_is_initialized = .false.

  end subroutine cg_drag_end

  !> The quantities of a column that do not depend on the wind: the breaking
  !> threshold and the reflection frequency.
  subroutine column_setup(iz0, rho, bf, z, fac, omc)

    integer, intent(in) :: iz0  !! k index of the source level
    real, dimension(:), intent(in) :: rho, bf, z
    !! density [kg/m3] and height [m] at levels 1..iz0+1, buoyancy frequency [1/s] at
    !! levels 1..iz0
    real, dimension(:), intent(out) :: fac, omc
    !! `fac`: breaking threshold [s/m]; `omc`: frequency above which a wave is reflected
    !! (eq. 21) [1/s], at levels 1..iz0

    real :: alp2
    integer :: k

    do k = 1, iz0
      fac(k) = rho(k)*kwv/(2.0*rho(iz0)*bf(k))
      ! alp2 = 1/(2 Hb)**2, with the density scale height Hb = dz/ln(rho(k+1)/rho(k)) of the
      ! layer below level k
      alp2 = (log(rho(k)/rho(k + 1))/(2.0*(z(k) - z(k + 1))))**2
      omc(k) = sqrt(bf(k)**2*k2/(k2 + alp2))
    end do

  end subroutine column_setup

  !> Computes the gravity-wave forcing of one wind component in one column, and the
  !> effective eddy diffusivity, above the source level.
  !>
  !> Based on the version implemented in SKYHI (27 Oct 1998) by M. J. Alexander
  !> and L. Bruhwiler. The equation numbers are those of Alexander and Dunkerton
  !> (1999).
  subroutine gwfc(iz0, p, ph, bf, fac, omc, u, gwf, ked)

    integer, intent(in) :: iz0  !! k index of the source level
    real, dimension(:), intent(in) :: p, ph
    !! pressure at full levels 1..kmax and at half levels 1..kmax+1 [Pa]
    real, dimension(:), intent(in) :: bf, fac, omc
    !! buoyancy frequency [1/s], breaking threshold and reflection frequency (`column_setup`)
    !! at levels 1..iz0
    real, dimension(:), intent(in) :: u  !! wind component at full levels [m/s]
    real, dimension(:), intent(out) :: gwf  !! gravity-wave forcing of `u` [m/s2]
    real, dimension(:), intent(out) :: ked  !! effective eddy diffusivity [m2/s]

    ! c0mu0: phase speed relative to the wind at the source level; B0: source spectrum
    real, dimension(nc) :: c0mu0, B0
    logical, dimension(nc) :: propagating
    ! stress(k): momentum flux [Pa] deposited in layer k, between full levels k and k+1
    ! (layer 0: the flux that reaches the model top); stress_ked: the same, weighted by
    ! c - u [Pa m/s], for the eddy diffusivity
    real, dimension(0:iz0) :: stress, stress_ked
    ! mom, mom_ked: the same, collected on the model levels
    real, dimension(iz0) :: mom, mom_ked
    real :: c, c0mu, Bsum, fm, fe, w_above, w_top
    integer :: k, n

!--------------------------------------------------------------------
!    the source spectrum B0 (eq. 17) and the sum of |B0| over all phase
!    speeds, Bsum. a wave with no phase speed relative to the source wind
!    carries no flux.
!-------------------------------------------------------------------
    Bsum = 0.
    do n = 1, nc
      c0mu0(n) = c0(n) - u(iz0)
      if (c0mu0(n) == 0.0) then
        B0(n) = 0.0
      else
        if (centre_on_source_wind) then
          c = c0mu0(n)
        else
          c = c0(n)
        end if
        B0(n) = sign(Bw*exp(-log(2.0)*c**2/cw**2), c0mu0(n))
        Bsum = Bsum + abs(B0(n))
      end if
    end do
    if (Bsum == 0.0) then
      call error_mesg('cg_drag_mod', ' zero flux input at source level', FATAL)
    end if

!----------------------------------------------------------------------
!    integrate upwards from the source level. a wave that breaks at level
!    k deposits its flux in the layer below it. with the intermittency of
!    eq. 19, eps = Bt/(rho(iz0)*Bsum), the flux sum fm [m2/s2] of the
!    waves that break is the stress rho(iz0)*eps*fm = Bt*fm/Bsum [Pa].
!----------------------------------------------------------------------
    stress = 0.0
    stress_ked = 0.0
    propagating = .true.
    do k = iz0, 1, -1
      fm = 0.
      fe = 0.
      do n = 1, nc
        if (.not. propagating(n)) cycle
        c0mu = c0(n) - u(k)
        if (c0mu == 0.) then
          ! no phase speed relative to the wind: removed (the breaking test would divide
          ! by zero); this is extremely rare
          propagating(n) = .false.
        else if (abs(c0mu)*kwv >= omc(k)) then
          ! total internal reflection (eq. 21)
          propagating(n) = .false.
        else if ((B0(n)/c0mu**3 >= fac(k)) .or. (c0mu0(n)*c0mu <= 0.0)) then
          ! breaking (eq. 24), or a critical level. a wave unstable already at the source
          ! level is removed without depositing (eq. 22).
          propagating(n) = .false.
          if (k < iz0) then
            fm = fm + B0(n)
            fe = fe + c0mu*B0(n)
          end if
        end if
      end do
      stress(k) = Bt*fm/Bsum
      stress_ked(k) = Bt*fe/Bsum
    end do

!----------------------------------------------------------------------
!    the waves still propagating above the top level reach the model top.
!----------------------------------------------------------------------
    fm = 0.
    fe = 0.
    do n = 1, nc
      if (propagating(n)) then
        fm = fm + B0(n)
        fe = fe + (c0(n) - u(1))*B0(n)
      end if
    end do
    stress(0) = Bt*fm/Bsum
    stress_ked(0) = Bt*fe/Bsum

!----------------------------------------------------------------------
!    collect the deposits on the model levels. the layer between full
!    levels k and k+1 is split at the half level k+1: the fraction w_above
!    of its pressure thickness lies in model level k, the rest in model
!    level k+1. the flux that reaches the model top is spread over levels
!    1..damp_level in proportion to their mass, a uniform acceleration.
!----------------------------------------------------------------------
    mom = 0.0
    mom_ked = 0.0
    do k = 1, iz0 - 1
      w_above = (ph(k + 1) - p(k))/(p(k + 1) - p(k))
      mom(k) = mom(k) + w_above*stress(k)
      mom(k + 1) = mom(k + 1) + (1.0 - w_above)*stress(k)
      mom_ked(k) = mom_ked(k) + w_above*stress_ked(k)
      mom_ked(k + 1) = mom_ked(k + 1) + (1.0 - w_above)*stress_ked(k)
    end do
    do k = 1, damp_level
      w_top = (ph(k + 1) - ph(k))/(ph(damp_level + 1) - ph(1))
      mom(k) = mom(k) + w_top*stress(0)
      mom_ked(k) = mom_ked(k) + w_top*stress_ked(0)
    end do

!----------------------------------------------------------------------
!    the acceleration of each model level is its momentum divided by its
!    mass (ph(k+1) - ph(k))/GRAV, and the eddy diffusivity (eq. 26) is
!    that divided by bf**2.
!----------------------------------------------------------------------
    gwf = 0.0
    ked = 0.0
    do k = 1, iz0
      gwf(k) = GRAV*mom(k)/(ph(k + 1) - ph(k))
      ked(k) = GRAV*mom_ked(k)/((ph(k + 1) - ph(k))*bf(k)**2)
    end do

  end subroutine gwfc

end module cg_drag_mod
