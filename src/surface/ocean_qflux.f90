!> Prescribed ocean heat fluxes (Q-fluxes) for the slab ocean.
!>
!> The Q-flux is the sum of the sources selected in `qflux_nml`: with
!> `qflux_source = 'analytic'`, the zonally symmetric meridional Q-flux of Merlis et al.
!> (2013, Part II); with `do_warmpool`, a zonally asymmetric pattern selected by
!> `warmpool_pattern`: a tropical warm pool alone, or the warm pool and the regional
!> patterns of Garfinkel et al. (2020) (Gulf Stream, Kuroshio, tropical Atlantic, ...).
!> The analytic and warm-pool fluxes do not change in time; they are computed once by
!> `ocean_qflux_init`, which is called by `ocean_model_mod`. With `qflux_source = 'file'`,
!> the meridional Q-flux is replaced by the field `qflux_file_var` of `INPUT/qflux_file.nc`,
!> regridded to the model grid and interpolated in time at every step by
!> `ocean_qflux_update`: a file with 1 record is used at all times, one with 2 to 12 records
!> is a climatology repeated every year, and one with more records is a dated time series
!> that must span the run. The `ocean/qflux_implied` output of a prescribed-SST run (a time
!> mean, or 12 monthly means) can be used as it is.
!>
!> By default the sum is applied over the ocean only (`qflux_ocean_only`), and a constant
!> is removed from the ocean points so that its global area integral is zero
!> (`qflux_remove_mean`): the Q-flux neither heats the land nor adds net heat. Setting both
!> switches to `.false.` reproduces v1, where the Q-flux was added at every point, land
!> included, without removing its mean. With a file, the warm-pool fluxes (if `do_warmpool`)
!> are added to the file field and both switches are applied again at every step.
!>
!> ## The zonally asymmetric patterns
!>
!> `warmpool_pattern = 'warm_pool'` is one wavelength of a cosine in longitude,
!> `warmpool_amp*(1 - x**4)*cos(warmpool_k*(lon - warmpool_phase))`, with
!> `x = (lat - warmpool_centr)/warmpool_width`, for `|x| <= 1` and
!> `warmpool_phase - 90/warmpool_k <= lon <= warmpool_phase + 270/warmpool_k` (taken round
!> the circle), and 0 elsewhere. The cosine is 0 at both ends of the window, so the pattern
!> is continuous for any `warmpool_k >= 1` (for `warmpool_k < 1` the wavelength is longer
!> than the circle: the window is cut 360 degrees east of its western end, where the pattern
!> jumps). In v1
!> (`warmpool_localization_choice = 1`) the cosine went round the whole latitude band, which
!> is discontinuous at 0E unless `warmpool_k` is an integer.
!>
!> `warmpool_pattern = 'garfinkel2020'` (the default; v1's `warmpool_localization_choice = 3`)
!> is the sum of the terms listed, one per row, in `garfinkel_waves` and
!> `garfinkel_gaussians` below. v1's choice 2 (a 2018 draft of it) was removed. A term is
!> of one of two forms, with `lon` [degE] and `lat` [degN]:
!>
!> * a wave: `A*(1 - x**p)*cos(k*(lon - lon0))` (or `sin`), `x = (lat - lat0)/w`;
!> * a product of two Gaussians,
!>   `A*exp(-(a1*lon + b1*lat - c1)**2/(2*s1**2))*exp(-(a2*lon + b2*lat - c2)**2/(2*s2**2))`:
!>   `a = 0` gives a Gaussian in latitude, `b = 0` one in longitude, and both nonzero a
!>   tilted one, centred on the line `a*lon + b*lat = c`; `s` is its width (sigma) [deg],
!>   in the units of `a*lon + b*lat`. A row with no second factor is a single Gaussian.
!>
!> A term is added only in its longitude window (which may wrap through 0E; some windows
!> move with `warmpool_phase`, although the Gaussian centres do not), in its latitude window
!> (inclusive; for some waves `|x| <= 1` instead), and if its gates hold, e.g.
!> `kuroshio_amp > 0.001`. Its amplitude is linear in the namelist knobs,
!> `A = c0 + c1*knob1 + c2*knob2 + c3*knob3`, and for the Australian terms it is multiplied
!> by `qflux_amp + warmpool_amp`. The comment above each row names the feature, its sign and
!> amplitude [W/m2] at the default knobs, its centre (lon, lat) and widths (sigma in lon and
!> lat, or along the tilted coordinates) [deg], and the knobs that scale it.
!>
!> Note that this pattern also uses `qflux_amp` (Caribbean and South America, the West
!> African coast, the Barents Sea) and `qflux_amp + warmpool_amp` (Australia and the
!> Indian Ocean), even with `qflux_source = 'none'`.
!>
!> Knob budgets of the `'garfinkel2020'` pattern, at the default knobs, over the whole sphere
!> (before `qflux_ocean_only` and `qflux_remove_mean` are applied): the global mean of
!> d(pattern)/d(knob) (its net heating per unit knob) and of |d(pattern)/d(knob)| [W/m2 per
!> W/m2]. A knob with a net much smaller than its |.| mostly moves heat around.
!>
!> | knob | net (mean d/dknob) | mean of abs(d/dknob) |
!> |---|---|---|
!> | `Hawaiiextra` | -0.00063 | 0.121 |
!> | `Pac_SPCZextra` | 0.00010 | 0.078 |
!> | `Pac_ITCZextra` | 0.00047 | 0.109 |
!> | `Africaextra` | 0.00000 | 0.018 |
!> | `Sampeextra` | 0.00019 | 0.079 |
!> | `gulf_amp` | -0.00007 | 0.005 |
!> | `kuroshio_amp` | 0.00003 | 0.012 |
!> | `trop_atlantic_amp` | 0.0047 | 0.054 |
!> | `warmpool_amp` | -0.0059 | 0.195 |
!> | `qflux_amp` | 0.0011 | 0.066 |
!>
!> (On a 0.25 deg grid; for `Pac_ITCZextra` and `Sampeextra`, whose terms are off at 0, just
!> above 0.) The global mean of the default pattern is 0.009 W/m2: every knob mostly moves
!> heat around.
!>
!> Namelist: `qflux_nml`
!> ([namelist reference](https://eddy-stanford.github.io/MiMA/Parameters/#qflux_nml)).
!>
!> References:
!>
!> * Merlis, T. M., T. Schneider, S. Bordoni, and I. Eisenman, 2013: Hadley circulation
!>   response to orbital precession. Part II: Subtropical continent. J. Climate, 26,
!>   754-771.
!> * Garfinkel, C. I., I. White, E. P. Gerber, M. Jucker, and M. Erez, 2020:
!>   The building blocks of Northern Hemisphere wintertime stationary waves.
!>   J. Climate, 33, 5611-5633, https://doi.org/10.1175/JCLI-D-19-0181.1.
module ocean_qflux_mod

  use constants_mod, only: pi
  use fms_mod, only: input_nml_file, check_nml_error, &
                     error_mesg, FATAL, mpp_pe, mpp_root_pe, stdlog
  use mpp_domains_mod, only: domain2d, mpp_global_sum, BITWISE_EXACT_SUM
  use time_manager_mod, only: time_type
  use mima_interpolator_mod, only: interpolate_type, interpolator_init, CONSTANT, interpolator
  use surface_types_mod, only: surface_grid_type, nml_sets

  implicit none
  private

  public :: ocean_qflux_init, ocean_qflux_update

  real ::    qflux_amp = 26., & !! [W/m2] amplitude of the meridional Q-flux; with
          !! `warmpool_pattern = 'garfinkel2020'` it also scales several regional terms (Caribbean
          !! and South America, the West African coast, the Barents Sea and, with `warmpool_amp`,
          !! Australia and the Indian Ocean), even with `qflux_source = 'none'`
          qflux_width = 16., & !! [deg] half-width of the meridional Q-flux
          warmpool_amp = 18., & !! [W/m2] amplitude of the warm pool; with `'garfinkel2020'`, also
          !! of its wave-8 eastern edge, and the West Pacific, Pacific, Hawaii and Australian terms
          !! are only applied if `warmpool_amp` > 0.001 (the Australian ones scale with
          !! `qflux_amp + warmpool_amp`)
          warmpool_width = 35., & !! [deg] latitudinal half-width of the warm pool (and, with
          !! `'garfinkel2020'`, of its eastern edge and of the Africa wave)
          warmpool_centr = 0., & !! [deg] central latitude of the warm pool (and, with
          !! `'garfinkel2020'`, of its eastern edge and of the Africa wave)
          warmpool_phase = 140.     !! [deg] longitude of the warm-pool maximum; with
  !! `'garfinkel2020'` it also moves the longitude windows (not the centres) of the warm-pool
  !! edge, Africa, Kuroshio and Gulf Stream terms

  integer :: gulf_k = 4       !! with `'garfinkel2020'`, the zonal wave number plus 1 of the
  !! heating near 67N between 318E and 18E (the default gives wave number 3)

  real :: warmpool_k = 1.66666, & !! zonal wave number of the warm pool; with `'warm_pool'` it sets
          !! the window (one wavelength); `'garfinkel2020'` uses the fixed window
          !! [`warmpool_phase` - 54, `warmpool_phase` + 162], which assumes `warmpool_k` = 5/3
          gulf_phase = 310., & !! [deg] with `'garfinkel2020'`, the phase of the heating near 67N:
          !! its maximum is at `gulf_phase` + 38
          gulf_amp = 70., & !! [W/m2] with `'garfinkel2020'`, the Gulf Stream amplitude; the Atlantic
          !! terms (tropical Atlantic, Caribbean and South America, near 67N, West African coast,
          !! Barents Sea) are only applied if `gulf_amp` > 0.001 and `trop_atlantic_amp` > 0.001
          kuroshio_amp = 40., &  !! [W/m2] with `'garfinkel2020'`, the Kuroshio amplitude; with
          !! `kuroshio_amp` <= 0.001, a reduced set of West Pacific terms is used
          trop_atlantic_amp = 50., &  !! [W/m2] with `'garfinkel2020'`, the tropical Atlantic
          !! amplitude, also of the heating near 67N; half of it scales the Africa wave and switches
          !! the Agulhas terms on (if > 0.001); the Arctic terms are only applied if it is > 0.001
          Pac_ITCZextra = 0., & !! [W/m2] extra cooling in the subtropical South-East Pacific and
          !! zonal heating near 55S (`'garfinkel2020'` only)
          Pac_SPCZextra = 0., & !! [W/m2] extra flux in the SPCZ and South Pacific, which modulates
          !! the SPCZ (`'garfinkel2020'` only)
          Africaextra = 0., &  !! [W/m2] extra flux near the Agulhas current (`'garfinkel2020'` only)
          Sampeextra = 0., & !! [W/m2] zonal dipole near 40S and 48S (`'garfinkel2020'` only)
          Hawaiiextra = 30.0 !! [W/m2] heating in the western North Pacific near 16N, which also
  !! changes many Pacific and Indian Ocean terms (`'garfinkel2020'` only; see the table)

  character(len=32) :: qflux_source = 'analytic'
  !! `'analytic'`: the meridional Q-flux of Merlis et al. (2013) (`qflux_amp`, `qflux_width`);
  !! `'none'`: no meridional Q-flux; `'file'`: read from `qflux_file` and interpolated in time
  !! at every step (1 record: constant; 2 to 12 records: a climatology repeated every year;
  !! more: a dated time series)
  logical :: do_warmpool = .true.
  !! add the zonally asymmetric Q-flux `warmpool_pattern`, whatever `qflux_source`
  character(len=32) :: warmpool_pattern = 'garfinkel2020'
  !! the zonally asymmetric Q-flux added with `do_warmpool`: `'garfinkel2020'`: the warm pool and
  !! the regional patterns of Garfinkel et al. (2020) (v1's `warmpool_localization_choice = 3`);
  !! `'warm_pool'`: the warm pool alone, one wavelength in longitude (v1's choice 1, which went
  !! round the whole latitude band). v1's choice 2 was removed. See `ocean_qflux.f90`.
  character(len=256) :: qflux_file = ''  !! with `qflux_source = 'file'`, Q-flux file name, without `.nc`, in `INPUT/`
  character(len=256) :: qflux_file_var = 'qflux'
  !! name of the Q-flux variable [W/m2] in `qflux_file` (empty: `qflux_file`); for the
  !! output of a prescribed-SST run, `'qflux_implied'`
  logical :: qflux_ocean_only = .true.
  !! apply the Q-flux over the ocean only (the default); `.false.`: also over land, as in v1
  !! (set both `qflux_ocean_only` and `qflux_remove_mean` to `.false.` to reproduce v1)
  logical :: qflux_remove_mean = .true.
  !! subtract a constant from the Q-flux at the ocean points so that its global area
  !! integral is zero (the default): the constant is the area integral of the applied Q-flux
  !! over all points (land included with `qflux_ocean_only = .false.`) divided by the ocean
  !! area; with `qflux_ocean_only` it is the ocean mean (the global mean without land).
  !! `.false.`: no mean removed, as in v1 (set both `qflux_ocean_only` and `qflux_remove_mean`
  !! to `.false.` to reproduce v1)

  logical :: qflux_initialized = .false.

  ! The Q-flux from a file (qflux_source = 'file'), which is recomputed at every step
  type(interpolate_type), save :: qflux_interp
  character(len=256) :: qflux_var  ! name of the field in the file
  type(domain2d), save :: qdomain  ! domain of the surface grid, for the global sums
  real, allocatable :: qarea(:, :)  ! relative grid-box area
  logical, allocatable :: qmask(:, :)  ! ocean points
  real, allocatable :: q_static(:, :)  ! the warm-pool fluxes (0 without do_warmpool) [W/m2]
  real :: qarea_ocean  ! global ocean area (relative)

  namelist /qflux_nml/ qflux_source, do_warmpool, qflux_file, qflux_file_var, &
    qflux_ocean_only, qflux_remove_mean, &
    qflux_amp, qflux_width, &
    warmpool_pattern, warmpool_amp, warmpool_width, warmpool_centr, &
    warmpool_k, warmpool_phase, &
    gulf_k, gulf_phase, gulf_amp, kuroshio_amp, trop_atlantic_amp, &
    Pac_ITCZextra, Sampeextra, &
    Pac_SPCZextra, Africaextra, Hawaiiextra

! ---------------------------------------------------------------------------------------
! The zonally asymmetric patterns as tables of terms (see the module documentation)
! ---------------------------------------------------------------------------------------
! The knobs, as indices into the knob vector built by warmpool (index 0 is a zero)
  integer, parameter :: K_QFLUX = 1, K_WARMPOOL = 2, K_GULF = 3, K_KUROSHIO = 4, &
                        K_TROPATL = 5, &
                        K_AFRICA = 6, &  ! africaamp = trop_atlantic_amp/2
                        K_HAWAII = 7, K_SPCZ = 8, K_ITCZ = 9, K_AFRX = 10, K_SAMPE = 11, &
                        K_QW = 12, &  ! qflux_amp + warmpool_amp
                        K_WPCENTR = 13, K_WPWIDTH = 14, K_WPK = 15, K_WPPHASE = 16, &
                        K_GULFK = 17, K_GULFPHASE = 18, NKNOB = 18
! The gates: conditions on the knobs for a term to be applied
  integer, parameter :: ALWAYS = 0, &
                        IF_KUROSHIO = 1, &  ! kuroshio_amp > 0.001
                        IF_NO_KUROSHIO = 2, &  ! not (kuroshio_amp > 0.001)
                        IF_WARMPOOL = 3, &  ! warmpool_amp > 0.001
                        IF_GULF = 4, &  ! gulf_amp > 0.001
                        IF_TROP_ATL = 5, &  ! trop_atlantic_amp > 0.001
                        IF_AFRICA = 6, &  ! trop_atlantic_amp/2 > 0.001
                        IF_SAMPE = 7, &  ! abs(Sampeextra) > 0.001
                        IF_ITCZ = 8, &  ! abs(Pac_ITCZextra) > 0.001
                        NGATE = 8
  integer, parameter :: ON_ALWAYS(2) = [ALWAYS, ALWAYS], ON_KURO(2) = [IF_KUROSHIO, ALWAYS], &
                        ON_WP(2) = [IF_WARMPOOL, ALWAYS], ON_WP_KURO(2) = [IF_WARMPOOL, IF_KUROSHIO], &
                        ON_WP_NOKURO(2) = [IF_WARMPOOL, IF_NO_KUROSHIO], ON_GULF(2) = [IF_GULF, ALWAYS], &
                        ON_TA(2) = [IF_TROP_ATL, ALWAYS], ON_TA_GULF(2) = [IF_TROP_ATL, IF_GULF], &
                        ON_AFRICA(2) = [IF_AFRICA, ALWAYS], ON_SAMPE(2) = [IF_SAMPE, ALWAYS], &
                        ON_ITCZ(2) = [IF_ITCZ, ALWAYS]
! The kinds of longitude window (lon >= lo always; "below" is lon < hi, or lon <= hi with hi_incl)
  integer, parameter :: LON_ALL = 0, &  ! everywhere
                        LON_IN = 1, &  ! lon >= lo .and. below
                        LON_OUT = 2, &  ! lon >= lo .or. below (a window through 0E)
                        LON_ABOVE = 3, &  ! lon >= lo
                        LON_BELOW = 4, &  ! below
                        LON_PERIODIC = 5  ! lo <= lon <= hi, with lon taken in [lo, lo + 360)
  integer, parameter :: COS_WAVE = 1, SIN_WAVE = 2

  !> A number linear in the knobs: `c0 + c1*knob(k1) + c2*knob(k2) + c3*knob(k3)`.
  type :: lin_t
    real :: c0 = 0.
    integer :: k1 = 0
    real :: c1 = 0.
    integer :: k2 = 0
    real :: c2 = 0.
    integer :: k3 = 0
    real :: c3 = 0.
  end type lin_t

  !> A longitude window [degE]: its bounds are `lo + knob(ref)` and `hi + knob(ref)`.
  type :: lonwin_t
    integer :: kind = LON_ALL
    real :: lo = 0., hi = 0.
    logical :: hi_incl = .false.
    integer :: ref = 0
  end type lonwin_t

  !> A wave term `amp*(1 - x**p)*cos(k*(lon - phase))` (or `sin`), `x = (lat - lat0)/width`,
  !> where `lat0`, `width`, `k` and `phase` [deg] are each a number plus `knob(..._ref)`.
  type :: wave_t
    integer :: gate(2) = ALWAYS
    type(lin_t) :: amp
    type(lonwin_t) :: lon
    real :: lat_lo = -90., lat_hi = 90.  ! latitude window [degN], inclusive
    real :: lat0 = 0.
    integer :: lat0_ref = 0
    real :: width = 0.
    integer :: width_ref = 0
    logical :: band = .false.  ! also require |x| <= 1
    integer :: p = 4
    integer :: fn = COS_WAVE
    real :: k = 0.
    integer :: k_ref = 0
    real :: phase = 0.
    integer :: phase_ref = 0
  end type wave_t

  !> A Gaussian factor `exp(-(a*lon + b*lat - c)**2/(2*s**2))` (lon, lat in degrees), with
  !> `c` plus `c_coef*knob(c_ref)`; `s = 0` means no factor (1).
  type :: gfac_t
    real :: a = 0., b = 0., c = 0., s = 0.
    integer :: c_ref = 0
    real :: c_coef = 0.
  end type gfac_t

  !> A Gaussian term `amp*f1*f2`, `amp` multiplied by `knob(scale)` if `scale` > 0.
  type :: gauss_t
    integer :: gate(2) = ALWAYS
    type(lin_t) :: amp
    integer :: scale = 0
    type(lonwin_t) :: lon
    real :: lat_lo = -90., lat_hi = 90.  ! latitude window [degN], inclusive
    type(gfac_t) :: f1, f2
  end type gauss_t

! Longitude windows shared by several terms
  type(lonwin_t), parameter :: &
    W_KUROSHIO = lonwin_t(LON_IN, -30., 130., .false., K_WPPHASE), &
    W_WPAC = lonwin_t(LON_IN, 70., 240.), &
    W_PAC = lonwin_t(LON_IN, 129., 290.), &
    W_INDPAC = lonwin_t(LON_IN, 50., 220.), &
    W_GULF1 = lonwin_t(LON_IN, 135., 195., .false., K_WPPHASE), &
    W_GULF2 = lonwin_t(LON_IN, 158., 218., .false., K_WPPHASE), &
    W_ATL = lonwin_t(LON_ABOVE, 275.), &
    W_GREEN = lonwin_t(LON_OUT, 318., 18.), &
    W_CARIB = lonwin_t(LON_IN, 250., 344.), &
    W_AGULHAS = lonwin_t(LON_IN, 2., 100.), &
    W_SATL = lonwin_t(LON_ABOVE, 290.), &
    W_NORW_W = lonwin_t(LON_ABOVE, 310.), &
    W_NORW_E = lonwin_t(LON_BELOW, 0., 30., .true.), &
    W_HUDSON = lonwin_t(LON_IN, 260., 310., .true.), &
    W_NATL = lonwin_t(LON_IN, 275., 335., .true.)
! The absent Gaussian factor (1)
  type(gfac_t), parameter :: NOF = gfac_t()

! The wave terms of 'garfinkel2020'. Columns: gate, amp, lon window, lat_lo, lat_hi,
! lat0, lat0_ref, width, width_ref, band, p, fn, k, k_ref, phase, phase_ref.
! Signs and amplitudes [W/m2] are at the default knobs (warmpool_phase = 140, gulf_phase = 310).
  type(wave_t), parameter :: garfinkel_waves(7) = [ &
! Indo-Pacific warm pool: heating 18 W/m2, peak at warmpool_phase (140E) on warmpool_centr,
! |lat - centr| <= warmpool_width, window [phase - 54, phase + 162] = [86E, 302E]; warmpool_amp
    wave_t(ON_ALWAYS, lin_t(0., K_WARMPOOL, 1.), lonwin_t(LON_IN, -54., 162., .true., K_WPPHASE), &
           -90., 90., 0., K_WPCENTR, 0., K_WPWIDTH, .true., 4, COS_WAVE, 0., K_WPK, 0., K_WPPHASE), &
! (cig, Mar 2019) eastern edge of the warm pool, wave-8 sine: 18 W/m2, cooling 257E-279.5E,
! heating 279.5E-302E, same latitudes; window [phase + 117, phase + 162]; warmpool_amp
    wave_t(ON_ALWAYS, lin_t(0., K_WARMPOOL, 1.), lonwin_t(LON_IN, 117., 162., .true., K_WPPHASE), &
           -90., 90., 0., K_WPCENTR, 0., K_WPWIDTH, .true., 4, SIN_WAVE, 8., 0, 139.5, K_WPPHASE), &
! (cig, May 2019) Africa and the western Indian Ocean, wave 5, (1 - x**2): 25 W/m2, heating
! 10E-46E (peak 28E), cooling 46E-82E, same latitudes; window [phase - 130, phase - 58];
! trop_atlantic_amp/2
    wave_t(ON_ALWAYS, lin_t(0., K_AFRICA, 1.), lonwin_t(LON_IN, -130., -58., .true., K_WPPHASE), &
           -90., 90., 0., K_WPCENTR, 0., K_WPWIDTH, .true., 2, COS_WAVE, 5., 0, -112., K_WPPHASE), &
! Second part of the Gulf Stream at the expense of South America (it replaced a Greenland
! perturbation): Irminger and Iceland seas, heating 64.3 W/m2, peak 348E (gulf_phase + 38),
! wave number gulf_k - 1 = 3, 318E-18E, |lat - 67| <= 10; trop_atlantic_amp*36/28
    wave_t(ON_TA_GULF, lin_t(0., K_TROPATL, 36./28.), W_GREEN, &
           -90., 90., 67., 0, 10., 0, .true., 4, COS_WAVE, -1., K_GULFK, 38., K_GULFPHASE), &
! Norwegian Sea at the expense of Africa (cig, May 2019): Barents Sea, wave 2, heating 68 W/m2,
! peak 30E 76N, 345E-75E, 71N-83N (x = (lat - 76)/6.5); qflux_amp*68/26
    wave_t(ON_TA_GULF, lin_t(0., K_QFLUX, 68./26.), lonwin_t(LON_OUT, 345., 75., .true.), &
           71., 83., 76., 0, 6.5, 0, .false., 4, COS_WAVE, 2., 0, 30., 0), &
! Wave 1 over the Arctic and Hudson Bay (cig, May 2019): 25 W/m2, heating peak 10E,
! cooling peak 190E, 69N-83N (x = (lat - 76)/7); fixed
    wave_t(ON_TA, lin_t(25.), lonwin_t(), &
           69., 83., 76., 0, 7., 0, .false., 4, COS_WAVE, 1., 0, 10., 0), &
! Norwegian Sea: heating 68.2 W/m2, wave 6 (one lobe), peak 2E 68N, 347E-17E, 60N-76N
! (x = (lat - 68)/8); fixed
    wave_t(ON_TA, lin_t(68.2), lonwin_t(LON_OUT, 347., 17., .true.), &
           60., 76., 68., 0, 8., 0, .false., 4, COS_WAVE, 6., 0, 2., 0)]

! The Gaussian terms of 'garfinkel2020'. Columns: gate, amp, scale, lon window, lat_lo, lat_hi,
! f1, f2, with a factor gfac_t(a, b, c, s[, c_ref, c_coef]) = exp(-(a*lon + b*lat - c)**2/(2*s**2)).
! Signs and amplitudes [W/m2] are at the default knobs; widths are sigma [deg].
  type(gauss_t), parameter :: garfinkel_gaussians(72) = [ &
! --- Kuroshio (cig, May 2019; Pacific sector is exp), window [phase - 30, phase + 130), 5N-47N
! subtropical NE Pacific: cooling -23.8, centre (241.5E, 26.5N), tilted, sigma 7 in lon + lat,
! 25 in lon - lat; kuroshio_amp*0.594
    gauss_t(ON_KURO, lin_t(0., K_KUROSHIO, -59.4/100.), 0, W_KUROSHIO, 5., 47., &
            gfac_t(1., 1., 268., 7.), gfac_t(1., -1., 215., 25.)), &
! Kuroshio: heating 40, centre (138.75E, 31.25N), tilted, sigma 10 in lon - 3 lat, 20 in
! lon + lat; kuroshio_amp
    gauss_t(ON_KURO, lin_t(0., K_KUROSHIO, 1.), 0, W_KUROSHIO, 5., 47., &
            gfac_t(1., -3., 45., 10.), gfac_t(1., 1., 170., 20.)), &
! --- More Kuroshio at the expense of South-East Asia (cig, May 2019), 70E-240E, 10S-60N
! South-East Asia and the NW Pacific: cooling -27.6, centre (140E, 19.7N), sigma (39, 7); fixed
    gauss_t(ON_WP_KURO, lin_t(-27.60), 0, W_WPAC, -10., 60., &
            gfac_t(1., 0., 140., 39.), gfac_t(0., 1., 19.7, 7.)), &
! Philippine Sea: cooling -5.2, centre (140E, 20N), sigma (8, 4); fixed
    gauss_t(ON_WP_KURO, lin_t(-5.2), 0, W_WPAC, -10., 60., &
            gfac_t(1., 0., 140., 8.), gfac_t(0., 1., 20., 4.)), &
! Kuroshio Extension: heating 35.4, centre (160E, 35N), sigma (20, 6); fixed
    gauss_t(ON_WP_KURO, lin_t(35.4), 0, W_WPAC, -10., 60., &
            gfac_t(1., 0., 160., 20.), gfac_t(0., 1., 35., 6.)), &
! Kuroshio, East China Sea side: heating 48.8, centre (132.75E, 27.25N), tilted, sigma 10 in
! lon - 3 lat (centred on 45 + Hawaiiextra/5), 20 in lon + lat; 49.5 - 0.0228*Hawaiiextra
    gauss_t(ON_WP_KURO, lin_t(49.5, K_HAWAII, -.0228), 0, W_WPAC, -10., 60., &
            gfac_t(1., -3., 45., 10., K_HAWAII, 1./5.), gfac_t(1., 1., 160., 20.)), &
! eastern equatorial Indian Ocean: heating 22.9, centre (90E, 0), sigma (12, 5); fixed
    gauss_t(ON_WP_KURO, lin_t(22.9), 0, W_WPAC, -10., 60., &
            gfac_t(1., 0., 90., 12.), gfac_t(0., 1., 0., 5.)), &
! --- the same without the Kuroshio (kuroshio_amp <= 0.001): off at the default knobs
! South-East Asia and the NW Pacific: cooling -5.3, centre (140E, 19.7N), sigma (39, 7); fixed
    gauss_t(ON_WP_NOKURO, lin_t(-5.3), 0, W_WPAC, -10., 60., &
            gfac_t(1., 0., 140., 39.), gfac_t(0., 1., 19.7, 7.)), &
! eastern equatorial Indian Ocean: heating 22.9, centre (90E, 0), sigma (12, 5); fixed
    gauss_t(ON_WP_NOKURO, lin_t(22.9), 0, W_WPAC, -10., 60., &
            gfac_t(1., 0., 90., 12.), gfac_t(0., 1., 0., 5.)), &
! --- South Pacific and SPCZ at the expense of the cold tongue (cig, May 2019), also includes
! the Australia bit from the Kuroshio above; 129E-290E, 78S-24N
! equatorial cold tongue 270E: cooling -32, centre (270E, 0), sigma (9, 3);
! -(50 - 0.28*Pac_SPCZextra - 0.6*Hawaiiextra)
    gauss_t(ON_WP, lin_t(-50., K_SPCZ, .28, K_HAWAII, .6), 0, W_PAC, -78., 24., &
            gfac_t(1., 0., 270., 9.), gfac_t(0., 1., 0., 3.)), &
! equatorial cold tongue 250E: cooling -32, centre (250E, 1S), sigma (9, 3); as above
    gauss_t(ON_WP, lin_t(-50., K_SPCZ, .28, K_HAWAII, .6), 0, W_PAC, -78., 24., &
            gfac_t(1., 0., 250., 9.), gfac_t(0., 1., -1., 3.)), &
! equatorial cold tongue 230E: cooling -32, centre (230E, 2S), sigma (9, 3); as above
    gauss_t(ON_WP, lin_t(-50., K_SPCZ, .28, K_HAWAII, .6), 0, W_PAC, -78., 24., &
            gfac_t(1., 0., 230., 9.), gfac_t(0., 1., -2., 3.)), &
! equatorial cold tongue 210E: cooling -21, centre (210E, 2S), sigma (9, 3);
! -(39 - 0.6*Hawaiiextra)
    gauss_t(ON_WP, lin_t(-39., K_HAWAII, .6), 0, W_PAC, -78., 24., &
            gfac_t(1., 0., 210., 9.), gfac_t(0., 1., -2., 3.)), &
! central equatorial Pacific 190E: cooling -15.7, centre (190E, 0), sigma (9, 3);
! -(36 - 0.676*Hawaiiextra)
    gauss_t(ON_WP, lin_t(-36., K_HAWAII, .676), 0, W_PAC, -78., 24., &
            gfac_t(1., 0., 190., 9.), gfac_t(0., 1., 0., 3.)), &
! central equatorial Pacific 170E: heating 4.3 (cooling for Hawaiiextra < 23.7), centre
! (170E, 0), sigma (9, 3); -(16 - 0.676*Hawaiiextra)
    gauss_t(ON_WP, lin_t(-16., K_HAWAII, .676), 0, W_PAC, -78., 24., &
            gfac_t(1., 0., 170., 9.), gfac_t(0., 1., 0., 3.)), &
! Chile coast (Humboldt Current): cooling -40, centre (287E, 25S), sigma (2, 9); fixed
    gauss_t(ON_WP, lin_t(-40.), 0, W_PAC, -78., 24., &
            gfac_t(1., 0., 287., 2.), gfac_t(0., 1., -25., 9.)), &
! Peru upwelling: cooling -15, centre (282E, 15S), sigma (5, 9); fixed
    gauss_t(ON_WP, lin_t(-15.), 0, W_PAC, -78., 24., &
            gfac_t(1., 0., 282., 5.), gfac_t(0., 1., -15., 9.)), &
! subtropical SE Pacific: cooling -10, centre (240E, 21S), sigma (40, 11);
! -(25 + Pac_ITCZextra + Pac_SPCZextra - 0.5*Hawaiiextra)
    gauss_t(ON_WP, lin_t(-25., K_ITCZ, -1., K_SPCZ, -1., K_HAWAII, .5), 0, W_PAC, -78., 24., &
            gfac_t(1., 0., 240., 40.), gfac_t(0., 1., -21., 11.)), &
! subtropical North Pacific 195E: cooling -17, centre (195E, 16N), sigma (13, 7);
! -(38 - 0.7*Hawaiiextra)
    gauss_t(ON_WP, lin_t(-38., K_HAWAII, .7), 0, W_PAC, -78., 24., &
            gfac_t(1., 0., 195., 13.), gfac_t(0., 1., 16., 7.)), &
! subtropical North Pacific 225E: cooling -30.4, centre (225E, 16N), sigma (13, 7);
! -(51.4 - 0.7*Hawaiiextra)
    gauss_t(ON_WP, lin_t(-51.4, K_HAWAII, .7), 0, W_PAC, -78., 24., &
            gfac_t(1., 0., 225., 13.), gfac_t(0., 1., 16., 7.)), &
! Southern Ocean, Pacific sector: heating 6.6, centre (220E, 57S), sigma (40, 15);
! 28.2 - (25/30)*0.8623*Hawaiiextra
    gauss_t(ON_WP, lin_t(28.2, K_HAWAII, -25./30.*.8623), 0, W_PAC, -78., 24., &
            gfac_t(1., 0., 220., 40.), gfac_t(0., 1., -57., 15.)), &
! SPCZ 165E: heating 14, centre (165E, 20S), sigma (20, 5); 14 + 0.9*Pac_SPCZextra
    gauss_t(ON_WP, lin_t(14., K_SPCZ, .9), 0, W_PAC, -78., 24., &
            gfac_t(1., 0., 165., 20.), gfac_t(0., 1., -20., 5.)), &
! SPCZ, south-eastern end: heating 16, centre (195E, 33S), sigma (20, 7);
! 16 + 0.9*Pac_SPCZextra
    gauss_t(ON_WP, lin_t(16., K_SPCZ, .9), 0, W_PAC, -78., 24., &
            gfac_t(1., 0., 195., 20.), gfac_t(0., 1., -33., 7.)), &
! East Australian Current: heating 50, centre (155E, 30S), sigma (3, 7);
! 50 + 0.9*Pac_SPCZextra
    gauss_t(ON_WP, lin_t(50., K_SPCZ, .9), 0, W_PAC, -78., 24., &
            gfac_t(1., 0., 155., 3.), gfac_t(0., 1., -30., 7.)), &
! east of New Zealand: heating 40, centre (180E, 40S), sigma (5, 5); 40 + 0.9*Pac_SPCZextra
    gauss_t(ON_WP, lin_t(40., K_SPCZ, .9), 0, W_PAC, -78., 24., &
            gfac_t(1., 0., 180., 5.), gfac_t(0., 1., -40., 5.)), &
! Southern Ocean, Amundsen-Bellingshausen: heating 16, centre (240E, 62S), sigma (30, 8);
! 41 - (25/30)*Hawaiiextra
    gauss_t(ON_WP, lin_t(41., K_HAWAII, -25./30.), 0, W_PAC, -78., 24., &
            gfac_t(1., 0., 240., 30.), gfac_t(0., 1., -62., 8.)), &
! NH ITCZ 7N, 180E: heating 30, centre (180E, 6.97N), sigma (13, 2); 60 - Hawaiiextra
    gauss_t(ON_WP, lin_t(60., K_HAWAII, -1.), 0, W_PAC, -78., 24., &
            gfac_t(1., 0., 180., 13.), gfac_t(0., 1., 6.97, 2.)), &
! NH ITCZ 7N, 210E: heating 17, centre (210E, 6.97N), sigma (13, 2); 47 - Hawaiiextra
    gauss_t(ON_WP, lin_t(47., K_HAWAII, -1.), 0, W_PAC, -78., 24., &
            gfac_t(1., 0., 210., 13.), gfac_t(0., 1., 6.97, 2.)), &
! NH ITCZ 7N, 240E: heating 15, centre (240E, 6.97N), sigma (13, 2); 45 - Hawaiiextra
    gauss_t(ON_WP, lin_t(45., K_HAWAII, -1.), 0, W_PAC, -78., 24., &
            gfac_t(1., 0., 240., 13.), gfac_t(0., 1., 6.97, 2.)), &
! West Pacific warm pool 3N: heating 4.0, centre (145E, 3N), sigma (14, 4);
! 19.5 + 0.4*Pac_SPCZextra - 0.5158*Hawaiiextra
    gauss_t(ON_WP, lin_t(19.5, K_SPCZ, .4, K_HAWAII, -.5158), 0, W_PAC, -78., 24., &
            gfac_t(1., 0., 145., 14.), gfac_t(0., 1., 3., 4.)), &
! West Pacific warm pool 7N: heating 24.5, centre (150E, 7N), sigma (13, 3);
! 40 + 0.35*Pac_SPCZextra - 0.5158*Hawaiiextra
    gauss_t(ON_WP, lin_t(40., K_SPCZ, .35, K_HAWAII, -.5158), 0, W_PAC, -78., 24., &
            gfac_t(1., 0., 150., 13.), gfac_t(0., 1., 7., 3.)), &
! --- Hawaii (cig, May 2019), 50E-220E, 16S-30N
! western North Pacific 16N: heating 27, centre (145E, 16N), sigma (20, 5); 0.9*Hawaiiextra
    gauss_t(ON_WP, lin_t(0., K_HAWAII, .9), 0, W_INDPAC, -16., 30., &
            gfac_t(1., 0., 145., 20.), gfac_t(0., 1., 16., 5.)), &
! --- South Pacific and South Indian Ocean at the expense of Australia (cig, May 2019),
! 50E-220E, 36S-10N; these scale with qflux_amp + warmpool_amp (44 at the default)
! Australia: cooling -44.9, centre (135E, 20S), sigma (15, 6); -1.02*(qflux_amp + warmpool_amp)
    gauss_t(ON_WP, lin_t(-1.02), K_QW, W_INDPAC, -36., 10., &
            gfac_t(1., 0., 135., 15.), gfac_t(0., 1., -20., 6.)), &
! eastern Australia: cooling -10.2, centre (147E, 27S), sigma (8, 7); -1.02*10/44*(q + wp)
    gauss_t(ON_WP, lin_t(-1.02*10./44.), K_QW, W_INDPAC, -36., 10., &
            gfac_t(1., 0., 147., 8.), gfac_t(0., 1., -27., 7.)), &
! NW Australia and the eastern Indian Ocean: heating 16.9, centre (120E, 20S), sigma (30, 6);
! 1.02*16.6/44*(q + wp)
    gauss_t(ON_WP, lin_t(1.02*16.6/44.), K_QW, W_INDPAC, -36., 10., &
            gfac_t(1., 0., 120., 30.), gfac_t(0., 1., -20., 6.)), &
! SE Indian Ocean, south of Java: cooling -14.8 (heating for Hawaiiextra < 19.8), centre
! (100E, 10S), sigma (10, 4); 1.02*(27.89 - 1.412*Hawaiiextra)/44*(q + wp)
    gauss_t(ON_WP, lin_t(1.02*27.89/44., K_HAWAII, -1.02*1.412/44.), K_QW, W_INDPAC, -36., 10., &
            gfac_t(1., 0., 100., 10.), gfac_t(0., 1., -10., 4.)), &
! Maritime Continent, equator: heating 5.0, centre (135E, 0), sigma (15, 4); 1.02*4.9/44*(q + wp)
    gauss_t(ON_WP, lin_t(1.02*4.9/44.), K_QW, W_INDPAC, -36., 10., &
            gfac_t(1., 0., 135., 15.), gfac_t(0., 1., 0., 4.)), &
! Timor and Arafura seas: 0 (heating for Pac_SPCZextra > 0), centre (130E, 8S), sigma (14, 9);
! 0.7137*Pac_SPCZextra
    gauss_t(ON_WP, lin_t(0., K_SPCZ, .7137), 0, W_INDPAC, -36., 10., &
            gfac_t(1., 0., 130., 14.), gfac_t(0., 1., -8., 9.)), &
! --- Gulf Stream (cig, May 2019; Atlantic sector is exp), 10N-52N
! Gulf Stream: heating 70, centre (296.7E, 38.3N), tilted, sigma 3 in lon - 2 lat, 25 in
! lon + lat; window [phase + 135, phase + 195); gulf_amp
    gauss_t(ON_GULF, lin_t(0., K_GULF, 1.), 0, W_GULF1, 10., 52., &
            gfac_t(1., -2., 220., 3.), gfac_t(1., 1., 335., 25.)), &
! subtropical NE Atlantic (Canary Current): cooling -63.9, centre (337.5E, 25N), tilted,
! sigma 3 in lon - lat/2, 7 in lat; window [phase + 158, phase + 218); gulf_amp*63.9/70
    gauss_t(ON_GULF, lin_t(0., K_GULF, -63.9/70.), 0, W_GULF2, 10., 52., &
            gfac_t(1., -.5, 325., 3.), gfac_t(0., 1., 25., 7.)), &
! --- More Gulf Stream at the expense of the tropical Atlantic (cig, May 2019), east of 275E,
! 35S-77N
! Gulf of Guinea 342E: cooling -50, centre (342E, 5S), sigma (9, 5); trop_atlantic_amp
    gauss_t(ON_TA_GULF, lin_t(0., K_TROPATL, -1.), 0, W_ATL, -35., 77., &
            gfac_t(1., 0., 342., 9.), gfac_t(0., 1., -5., 5.)), &
! Gulf of Guinea 360E: cooling -50, centre (360E, 5S), sigma (8, 5); trop_atlantic_amp
    gauss_t(ON_TA_GULF, lin_t(0., K_TROPATL, -1.), 0, W_ATL, -35., 77., &
            gfac_t(1., 0., 360., 8.), gfac_t(0., 1., -5., 5.)), &
! tropical South Atlantic: cooling -12.6, centre (345E, 16S), sigma (16, 8); fixed
    gauss_t(ON_TA_GULF, lin_t(-12.6), 0, W_ATL, -35., 77., &
            gfac_t(1., 0., 345., 16.), gfac_t(0., 1., -16., 8.)), &
! North Atlantic Current: heating 54.7, centre (323.3E, 51.7N), tilted, sigma 10 in
! lon - 2 lat, 30 in lon + lat; trop_atlantic_amp*30.65/28
    gauss_t(ON_TA_GULF, lin_t(0., K_TROPATL, 30.65/28.), 0, W_ATL, -35., 77., &
            gfac_t(1., -2., 220., 10.), gfac_t(1., 1., 375., 30.)), &
! --- Second part of the Gulf Stream at the expense of South America, 318E-18E, 35S-15N
! (see also the wave near 67N)
! Gulf of Guinea 0E: cooling -50, centre (0E, 5S), sigma (8, 5); trop_atlantic_amp
    gauss_t(ON_TA_GULF, lin_t(0., K_TROPATL, -1.), 0, W_GREEN, -35., 15., &
            gfac_t(1., 0., 0., 8.), gfac_t(0., 1., -5., 5.)), &
! tropical South Atlantic: cooling -12.6, centre (-15E, 16S), sigma (16, 8); fixed
    gauss_t(ON_TA_GULF, lin_t(-12.6), 0, W_GREEN, -35., 15., &
            gfac_t(1., 0., -15., 16.), gfac_t(0., 1., -16., 8.)), &
! --- Caribbean and South America at the expense of the tropical North and South Atlantic
! (cig, May 2019), 250E-344E, 35S-40N
! subtropical South America: cooling -23.9, centre (290E, 20S), sigma (20, 7); qflux_amp*0.92
    gauss_t(ON_TA_GULF, lin_t(0., K_QFLUX, -.92), 0, W_CARIB, -35., 40., &
            gfac_t(1., 0., 290., 20.), gfac_t(0., 1., -20., 7.)), &
! tropical North Atlantic: cooling -16.8, centre (325E, 19.5N), sigma (22, 8); fixed
    gauss_t(ON_TA_GULF, lin_t(-16.8), 0, W_CARIB, -35., 40., &
            gfac_t(1., 0., 325., 22.), gfac_t(0., 1., 19.5, 8.)), &
! Gulf of Mexico: heating 31.2, centre (270E, 22N), sigma (7, 5); qflux_amp*1.2
    gauss_t(ON_TA_GULF, lin_t(0., K_QFLUX, 1.2), 0, W_CARIB, -35., 40., &
            gfac_t(1., 0., 270., 7.), gfac_t(0., 1., 22., 5.)), &
! north-western South America (Ecuador): heating 41.1, centre (283E, 0), sigma (5, 6);
! qflux_amp*1.58
    gauss_t(ON_TA_GULF, lin_t(0., K_QFLUX, 1.58), 0, W_CARIB, -35., 40., &
            gfac_t(1., 0., 283., 5.), gfac_t(0., 1., 0., 6.)), &
! Amazon: heating 27.7, centre (304E, 2S), sigma (6, 7); qflux_amp*1.06415
    gauss_t(ON_TA_GULF, lin_t(0., K_QFLUX, 1.06415), 0, W_CARIB, -35., 40., &
            gfac_t(1., 0., 304., 6.), gfac_t(0., 1., -2., 7.)), &
! Peru and the Andes: heating 22.1, centre (284E, 10S), sigma (5, 6); qflux_amp*0.85
    gauss_t(ON_TA_GULF, lin_t(0., K_QFLUX, .85), 0, W_CARIB, -35., 40., &
            gfac_t(1., 0., 284., 5.), gfac_t(0., 1., -10., 6.)), &
! north-eastern Brazil: heating 16.4, centre (317E, 6S), sigma (5, 4); qflux_amp*0.63
    gauss_t(ON_TA_GULF, lin_t(0., K_QFLUX, .63), 0, W_CARIB, -35., 40., &
            gfac_t(1., 0., 317., 5.), gfac_t(0., 1., -6., 4.)), &
! Atlantic ITCZ: heating 42.5, centre (325E, 4.2N), sigma (11, 2); fixed
    gauss_t(ON_TA_GULF, lin_t(42.54), 0, W_CARIB, -35., 40., &
            gfac_t(1., 0., 325., 11.), gfac_t(0., 1., 4.2, 2.)), &
! --- Agulhas Current at the expense of Africa (cig, May 2019), 2E-100E, 60S-35N; applied if
! trop_atlantic_amp/2 > 0.001
! Sahel and North Africa: cooling -30, centre (28E, 18N), sigma (10, 7.07); fixed
    gauss_t(ON_AFRICA, lin_t(-30.), 0, W_AGULHAS, -60., 35., &
            gfac_t(1., 0., 28., 10.), gfac_t(0., 1., 18., sqrt(50.))), &
! southern Africa: cooling -30, centre (28E, 18S), sigma (10, 7.75); fixed
    gauss_t(ON_AFRICA, lin_t(-30.), 0, W_AGULHAS, -60., 35., &
            gfac_t(1., 0., 28., 10.), gfac_t(0., 1., -18., sqrt(60.))), &
! Angola coast (Benguela): cooling -38.5, centre (11E, 15S), sigma (2, 10);
! -(38.5 + 0.7709*Africaextra)
    gauss_t(ON_AFRICA, lin_t(-38.5, K_AFRX, -.7709), 0, W_AGULHAS, -60., 35., &
            gfac_t(1., 0., 11., 2.), gfac_t(0., 1., -15., 10.)), &
! Agulhas Return Current: heating 83, centre (50E, 40S), sigma (25, 4); 83 + Africaextra
    gauss_t(ON_AFRICA, lin_t(83., K_AFRX, 1.), 0, W_AGULHAS, -60., 35., &
            gfac_t(1., 0., 50., 25.), gfac_t(0., 1., -40., 4.)), &
! Southern Ocean south of Africa: cooling -64.2, centre (50E, 48S), sigma (20, 4);
! -(64.22 + 1.3*Africaextra)
    gauss_t(ON_AFRICA, lin_t(-64.22, K_AFRX, -1.3), 0, W_AGULHAS, -60., 35., &
            gfac_t(1., 0., 50., 20.), gfac_t(0., 1., -48., 4.)), &
! Agulhas Current: heating 38, centre (38.2E, 28.2S), tilted, sigma 4 in lon - 2 lat/3,
! 15 in lon + lat; 38 + Africaextra/3
    gauss_t(ON_AFRICA, lin_t(38., K_AFRX, 1./3.), 0, W_AGULHAS, -60., 35., &
            gfac_t(1., -2./3., 57., 4.), gfac_t(1., 1., 10., 15.)), &
! Gulf of Guinea and the Congo: heating 20, centre (14E, 0), sigma (5.48, 7.07); fixed
    gauss_t(ON_AFRICA, lin_t(20.), 0, W_AGULHAS, -60., 35., &
            gfac_t(1., 0., 14., sqrt(30.)), gfac_t(0., 1., 0., sqrt(50.))), &
! East Africa: heating 11, centre (36E, 0), sigma (5.48, 7.07); fixed
    gauss_t(ON_AFRICA, lin_t(11.), 0, W_AGULHAS, -60., 35., &
            gfac_t(1., 0., 36., sqrt(30.)), gfac_t(0., 1., 0., sqrt(50.))), &
! --- Zonal dipole (cig, Feb 2020; "Agulhas current at the expense of Africa"), 61S-30S;
! applied if abs(Sampeextra) > 0.001 (off at the default)
! zonal 40S: heating for Sampeextra > 0, centre 40S, sigma 4; 0.8822*Sampeextra
    gauss_t(ON_SAMPE, lin_t(0., K_SAMPE, .8822), 0, lonwin_t(), -61., -30., &
            gfac_t(0., 1., -40., 4.), NOF), &
! zonal 48S: cooling for Sampeextra > 0, centre 48S, sigma 4; -Sampeextra
    gauss_t(ON_SAMPE, lin_t(0., K_SAMPE, -1.), 0, lonwin_t(), -61., -30., &
            gfac_t(0., 1., -48., 4.), NOF), &
! --- Zonally symmetric heating (cig, Jun 2020), 65S-45S; applied if abs(Pac_ITCZextra) > 0.001
! Southern Ocean 55S: heating for Pac_ITCZextra > 0, centre 55S, sigma 7; 0.74537*Pac_ITCZextra
    gauss_t(ON_ITCZ, lin_t(0., K_ITCZ, .74537), 0, lonwin_t(), -65., -45., &
            gfac_t(0., 1., -55., 7.), NOF), &
! --- Dipole in the South Atlantic (cig, Aug 2019), east of 290E, 61S-30S
! Brazil-Malvinas confluence: heating 37.4, centre (323E, 36S), sigma (11, 4); fixed
    gauss_t(ON_ALWAYS, lin_t(37.4), 0, W_SATL, -61., -30., &
            gfac_t(1., 0., 323., 11.), gfac_t(0., 1., -36., 4.)), &
! Malvinas Current: cooling -40, centre (311E, 45S), sigma (11, 4); fixed
    gauss_t(ON_ALWAYS, lin_t(-40.), 0, W_SATL, -61., -30., &
            gfac_t(1., 0., 311., 11.), gfac_t(0., 1., -45., 4.)), &
! --- Norwegian Sea at the expense of Africa (cig, May 2019), 10N-35N (see also the Barents
! Sea wave)
! West African coast, west of 0E: cooling -14.5, centre (357E, 20N), sigma (20, 7), east of
! 310E; qflux_amp*14.5/26
    gauss_t(ON_TA_GULF, lin_t(0., K_QFLUX, -14.5/26.), 0, W_NORW_W, 10., 35., &
            gfac_t(1., 0., 357., 20.), gfac_t(0., 1., 20., 7.)), &
! West African coast, east of 0E: the same Gaussian, centred on -3E, west of 30E
    gauss_t(ON_TA_GULF, lin_t(0., K_QFLUX, -14.5/26.), 0, W_NORW_E, 10., 35., &
            gfac_t(1., 0., -3., 20.), gfac_t(0., 1., 20., 7.)), &
! --- Arctic and Hudson Bay (cig, May 2019; see also the waves 1 and 6)
! Hudson Bay: cooling -38, centre (278.7E, 63.3N), tilted, sigma 10 in lon - 2 lat, 20 in
! lon + lat; 260E-310E, 55N-85N; fixed
    gauss_t(ON_TA, lin_t(-38.), 0, W_HUDSON, 55., 85., &
            gfac_t(1., -2., 152., 10.), gfac_t(1., 1., 342., 20.)), &
! southern Hudson Bay: cooling -100, centre (275E, 58N), sigma (5, 4); 260E-310E, 55N-85N; fixed
    gauss_t(ON_TA, lin_t(-100.), 0, W_HUDSON, 55., 85., &
            gfac_t(1., 0., 275., 5.), gfac_t(0., 1., 58., 4.)), &
! Gulf Stream: heating 10.8, centre (296.7E, 38.3N), tilted, sigma 10 in lon - 2 lat, 25 in
! lon + lat; 275E-335E, 10N-52N; fixed
    gauss_t(ON_TA, lin_t(10.8), 0, W_NATL, 10., 52., &
            gfac_t(1., -2., 220., 10.), gfac_t(1., 1., 335., 25.))]

contains

!########################################################
  !> Reads `qflux_nml` and returns the Q-flux on this PE's grid.
  !>
  !> The raw Q-flux is the sum of the meridional Q-flux (`qflux_source = 'analytic'`) or
  !> the file field at `Time` (`qflux_source = 'file'`) and the zonally asymmetric fluxes
  !> (`do_warmpool`); it is zero if none is selected. With `qflux_ocean_only` it is then set
  !> to 0 off the ocean, and with `qflux_remove_mean` its global area integral (over all
  !> points, land included where it is applied there), divided by the ocean area, is
  !> subtracted at the ocean points, so that the global net Q-flux is zero. The global sums
  !> are bitwise independent of the domain decomposition. Statistics of the raw and the
  !> applied Q-flux are written to the log file. With a file, `ocean_qflux_update` repeats
  !> this at every step.
  subroutine ocean_qflux_init(grid, mask, Time, flux)
    type(surface_grid_type), intent(in) :: grid  !! the surface grid
    logical, dimension(:, :), intent(in) :: mask  !! ocean points
    type(time_type), intent(in) :: Time  !! current time (for the file field)
    real, dimension(:, :), intent(out) :: flux
    !! Q-flux applied at each point (ocean heat flux convergence; 0 over land with
    !! `qflux_ocean_only`) [W/m2]
    integer :: ierr, io
    real, dimension(size(flux, 1), size(flux, 2)) :: ocean, land_q
    real :: area_total, area_ocean, area_land_q, mean_ocean, mean_land_q, removed

    call check_removed_keys()
    read (input_nml_file, nml=qflux_nml, iostat=io)
    ierr = check_nml_error(io, 'qflux_nml')
    if (mpp_pe() == mpp_root_pe()) write (stdlog(), nml=qflux_nml)

    select case (trim(qflux_source))
    case ('analytic', 'none')
    case ('file')
      if (len_trim(qflux_file) == 0) &
        call error_mesg('ocean_qflux_init', "qflux_source = 'file' needs qflux_file", FATAL)
    case default
      call error_mesg('ocean_qflux_init', 'unknown qflux_source = "'//trim(qflux_source)// &
                      '"; use analytic, none or file', FATAL)
    end select
    select case (trim(warmpool_pattern))
    case ('garfinkel2020')
    case ('warm_pool')
      if (do_warmpool .and. .not. warmpool_k > 0.) &
        call error_mesg('ocean_qflux_init', "warmpool_pattern = 'warm_pool' needs warmpool_k > 0", FATAL)
    case default
      call error_mesg('ocean_qflux_init', 'unknown warmpool_pattern = "'//trim(warmpool_pattern)// &
                      '"; use garfinkel2020 (v1 warmpool_localization_choice = 3) or warm_pool '// &
                      '(v1 choice 1, now one wavelength); v1 choice 2 was removed', FATAL)
    end select

    qflux_initialized = .true.

    ocean = merge(1., 0., mask)
    area_total = mpp_global_sum(grid%domain, grid%area, flags=BITWISE_EXACT_SUM)
    area_ocean = mpp_global_sum(grid%domain, grid%area*ocean, flags=BITWISE_EXACT_SUM)

    flux = 0.
    if (trim(qflux_source) == 'file') then
! Q-flux from a file, plus the (static) warm-pool fluxes; kept for ocean_qflux_update
      qflux_var = qflux_file_var
      if (len_trim(qflux_var) == 0) qflux_var = qflux_file
      call interpolator_init(qflux_interp, trim(qflux_file)//'.nc', grid%lon_bnd, grid%lat_bnd, &
                             data_names=(/qflux_var/), data_out_of_bounds=(/CONSTANT/))
      qdomain = grid%domain
      allocate (qarea(size(flux, 1), size(flux, 2)), qmask(size(flux, 1), size(flux, 2)), &
                q_static(size(flux, 1), size(flux, 2)))
      qarea = grid%area
      qmask = mask
      qarea_ocean = area_ocean
      q_static = 0.
      if (do_warmpool) call warmpool(grid%lon_bnd, grid%lat_bnd, q_static)
      call interpolator(qflux_interp, Time, flux, trim(qflux_var))
      flux = flux + q_static
    else
!mj q-flux as in Merlis et al (2013) [Part II]
      if (trim(qflux_source) == 'analytic') call merlis_qflux(grid%lat_bnd, flux)
!mj q-flux to create a tropical temperature perturbation
      if (do_warmpool) call warmpool(grid%lon_bnd, grid%lat_bnd, flux)
    end if

! Statistics of the raw Q-flux (they do not change flux): its mean over the ocean, and the
! land points where it is nonzero, with their area and mean
    land_q = merge(1., 0., (.not. mask) .and. flux /= 0.)
    area_land_q = mpp_global_sum(grid%domain, grid%area*land_q, flags=BITWISE_EXACT_SUM)
    mean_ocean = 0.
    if (area_ocean > 0.) &
      mean_ocean = mpp_global_sum(grid%domain, grid%area*ocean*flux, flags=BITWISE_EXACT_SUM)/area_ocean
    mean_land_q = 0.
    if (area_land_q > 0.) &
      mean_land_q = mpp_global_sum(grid%domain, grid%area*land_q*flux, flags=BITWISE_EXACT_SUM)/area_land_q

    call apply_switches(grid%domain, grid%area, mask, area_ocean, flux, removed)

    if (mpp_pe() == mpp_root_pe()) then
      write (stdlog(), '(a)') 'ocean_qflux_init: Q-flux statistics (area-weighted)'
      if (trim(qflux_source) == 'file') &
        write (stdlog(), '(a)') '  (Q-flux from '//trim(qflux_file)//'.nc at the initial time; '// &
        'it is recomputed at every step)'
      write (stdlog(), '(a,es14.6)') '  raw Q-flux mean over the ocean [W/m2]:          ', mean_ocean
      write (stdlog(), '(a,es14.6)') '  removed from the ocean points (net Q / ocean area) [W/m2]: ', removed
      write (stdlog(), '(a,es14.6)') '  land area with a nonzero raw Q-flux [fraction of the globe]: ', &
        area_land_q/area_total
      if (area_total > area_ocean) &
        write (stdlog(), '(a,es14.6)') '  land area with a nonzero raw Q-flux [fraction of the land]:  ', &
        area_land_q/(area_total - area_ocean)
      write (stdlog(), '(a,es14.6)') '  raw Q-flux mean over that land area [W/m2]:     ', mean_land_q
      if (qflux_ocean_only) then
        write (stdlog(), '(a)') '  qflux_ocean_only: the Q-flux over land is set to 0'
      else
        write (stdlog(), '(a)') '  qflux_ocean_only = .false.: the Q-flux is also applied over land (v1)'
      end if
    end if

  end subroutine ocean_qflux_init
!########################################################

  !> Stops with a migration message if `input.nml` still sets a `qflux_nml` key that was
  !> removed in v2.0 (the namelist read would otherwise only report an unknown variable).
  !> A mention of the key in a comment is ignored.
  subroutine check_removed_keys()

    if (nml_sets('warmpool_localization_choice')) &
      call error_mesg('ocean_qflux_init', 'qflux_nml: warmpool_localization_choice was replaced '// &
                      "by warmpool_pattern ('garfinkel2020' for choice 3, 'warm_pool' for choice 1; "// &
                      'choice 2 was removed); see Migration_v2.md', FATAL)
    if (nml_sets('north_sea_heat')) &
      call error_mesg('ocean_qflux_init', 'qflux_nml: north_sea_heat was removed with '// &
                      'warmpool_localization_choice = 2; delete it (see Migration_v2.md)', FATAL)

  end subroutine check_removed_keys
!########################################################

  !> Recomputes the Q-flux at `Time` when it is read from a file (`qflux_source = 'file'`);
  !> otherwise does nothing, and returns `.false.`.
  !>
  !> The file field is interpolated to `Time`, the warm-pool fluxes are added (if
  !> `do_warmpool`), and `qflux_ocean_only` and `qflux_remove_mean` are applied as in
  !> `ocean_qflux_init`, with bitwise reproducible global sums.
  logical function ocean_qflux_update(Time, flux) result(updated)
    type(time_type), intent(in) :: Time  !! time of the step
    real, dimension(:, :), intent(inout) :: flux
    !! Q-flux applied at each point; unchanged unless it is read from a file [W/m2]
    real :: removed

    updated = trim(qflux_source) == 'file'
    if (.not. updated) return
    call interpolator(qflux_interp, Time, flux, trim(qflux_var))
    flux = flux + q_static
    call apply_switches(qdomain, qarea, qmask, qarea_ocean, flux, removed)

  end function ocean_qflux_update
!########################################################

  !> Applies `qflux_ocean_only` and `qflux_remove_mean` to the raw Q-flux `flux`.
  !>
  !> With `qflux_ocean_only`, `flux` is set to 0 off the ocean. With `qflux_remove_mean`,
  !> its global area integral (over all points), divided by the ocean area, is subtracted
  !> at the ocean points; the global sum is bitwise independent of the domain
  !> decomposition.
  subroutine apply_switches(domain, area, mask, area_ocean, flux, removed)
    type(domain2d), intent(in) :: domain  !! domain of the surface grid
    real, dimension(:, :), intent(in) :: area  !! relative grid-box area
    logical, dimension(:, :), intent(in) :: mask  !! ocean points
    real, intent(in) :: area_ocean  !! global ocean area (relative, as `area`)
    real, dimension(:, :), intent(inout) :: flux  !! raw Q-flux on input, applied Q-flux on output [W/m2]
    real, intent(out) :: removed  !! the constant subtracted at the ocean points [W/m2]

! Ocean only: no Q-flux over land
    if (qflux_ocean_only) then
      where (.not. mask) flux = 0.
    end if
! Make the global net Q-flux zero: subtract, at the ocean points, the area integral of the
! applied Q-flux over all points (land included where it is applied there) divided by the
! ocean area. With qflux_ocean_only this is the ocean mean (the whole sphere without land).
    removed = 0.
    if (qflux_remove_mean .and. area_ocean > 0.) then
      removed = mpp_global_sum(domain, area*flux, flags=BITWISE_EXACT_SUM)/area_ocean
      where (mask) flux = flux - removed
    end if

  end subroutine apply_switches
!########################################################

  !> Subtracts the meridional Q-flux of Merlis et al. (2013, Part II) from `flux`.
  subroutine merlis_qflux(latb, flux)
    real, dimension(:), intent(in)    :: latb   !! latitudes of the cell boundaries [rad]
    real, dimension(:, :), intent(inout) :: flux   !! total ocean heat flux [W/m2]

    integer :: j
    real :: lat, coslat

    if (.not. qflux_initialized) then
      call error_mesg('merlis_qflux', 'ocean_qflux_mod is not initialized', FATAL)
    end if

    do j = 1, size(latb) - 1
      lat = 0.5*(latb(j + 1) + latb(j))
      coslat = cos(lat)
      lat = lat*180./pi
      flux(:, j) = flux(:, j) - qflux_amp*(1 - 2.*lat**2/qflux_width**2)* &
                   exp(-((lat)**2/(qflux_width)**2))/coslat
    end do

  end subroutine merlis_qflux

!########################################################

  !> Adds the zonally asymmetric Q-flux `warmpool_pattern` to `flux` (see the module
  !> documentation).
  subroutine warmpool(lonb, latb, flux)
    real, dimension(:), intent(in)   :: lonb, latb  !! longitudes and latitudes of the cell boundaries [rad]
    real, dimension(:, :), intent(inout):: flux       !! total ocean heat flux [W/m2]

    real :: kn(0:NKNOB)
    logical :: on(0:NGATE)
    type(wave_t) :: pool(1)

! The knob vector (kn(0) = 0 is the reference of terms that use no knob)
    kn = 0.
    kn(K_QFLUX) = qflux_amp
    kn(K_WARMPOOL) = warmpool_amp
    kn(K_GULF) = gulf_amp
    kn(K_KUROSHIO) = kuroshio_amp
    kn(K_TROPATL) = trop_atlantic_amp
    kn(K_AFRICA) = trop_atlantic_amp*2/4.
    kn(K_HAWAII) = Hawaiiextra
    kn(K_SPCZ) = Pac_SPCZextra
    kn(K_ITCZ) = Pac_ITCZextra
    kn(K_AFRX) = Africaextra
    kn(K_SAMPE) = Sampeextra
    kn(K_QW) = qflux_amp + warmpool_amp
    kn(K_WPCENTR) = warmpool_centr
    kn(K_WPWIDTH) = warmpool_width
    kn(K_WPK) = warmpool_k
    kn(K_WPPHASE) = warmpool_phase
    kn(K_GULFK) = real(gulf_k)
    kn(K_GULFPHASE) = gulf_phase
! The gates
    on(ALWAYS) = .true.
    on(IF_KUROSHIO) = kuroshio_amp > 0.001
    on(IF_NO_KUROSHIO) = .not. on(IF_KUROSHIO)
    on(IF_WARMPOOL) = warmpool_amp > 0.001
    on(IF_GULF) = gulf_amp > 0.001
    on(IF_TROP_ATL) = trop_atlantic_amp > 0.001
    on(IF_AFRICA) = kn(K_AFRICA) > 0.001
    on(IF_SAMPE) = abs(Sampeextra) > 0.001
    on(IF_ITCZ) = abs(Pac_ITCZextra) > 0.001

    select case (trim(warmpool_pattern))
    case ('garfinkel2020')
      call add_waves(garfinkel_waves, kn, on, lonb, latb, flux)
      call add_gaussians(garfinkel_gaussians, kn, on, lonb, latb, flux)
    case ('warm_pool')
! One wavelength, [phase - 90/k, phase + 270/k], round the circle
      pool(1) = wave_t(ON_ALWAYS, lin_t(0., K_WARMPOOL, 1.), &
                       lonwin_t(LON_PERIODIC, -90./warmpool_k, 270./warmpool_k, .true., K_WPPHASE), &
                       -90., 90., 0., K_WPCENTR, 0., K_WPWIDTH, .true., 4, COS_WAVE, 0., K_WPK, 0., K_WPPHASE)
      call add_waves(pool, kn, on, lonb, latb, flux)
    end select
  end subroutine warmpool

  !> The value of `x` for the knobs `kn`.
  pure real function lin_value(x, kn)
    type(lin_t), intent(in) :: x
    real, intent(in) :: kn(0:)
    lin_value = x%c0 + x%c1*kn(x%k1) + x%c2*kn(x%k2) + x%c3*kn(x%k3)
  end function lin_value

  !> Whether the longitude `lon` [rad] is in the window `w` with bounds `lo`, `hi` [rad];
  !> `lon_e` is `lon`, or for a periodic window `lon` shifted by a multiple of 2 pi into
  !> [`lo`, `lo` + 2 pi).
  pure subroutine lon_window(w, lo, hi, lon, inside, lon_e)
    type(lonwin_t), intent(in) :: w
    real, intent(in) :: lo, hi, lon
    logical, intent(out) :: inside
    real, intent(out) :: lon_e
    logical :: below

    lon_e = lon
    if (w%hi_incl) then
      below = lon <= hi
    else
      below = lon < hi
    end if
    select case (w%kind)
    case (LON_IN)
      inside = lon >= lo .and. below
    case (LON_OUT)
      inside = lon >= lo .or. below
    case (LON_ABOVE)
      inside = lon >= lo
    case (LON_BELOW)
      inside = below
    case (LON_PERIODIC)
      if (lon < lo .or. lon >= lo + 2.*pi) lon_e = lo + modulo(lon - lo, 2.*pi)
      inside = lon_e <= hi
    case default
      inside = .true.
    end select
  end subroutine lon_window

  !> Adds the wave terms `rows` to `flux`, for the knobs `kn` and the gates `on`.
  subroutine add_waves(rows, kn, on, lonb, latb, flux)
    type(wave_t), intent(in) :: rows(:)
    real, intent(in) :: kn(0:)
    logical, intent(in) :: on(0:)
    real, dimension(:), intent(in) :: lonb, latb  ! cell boundaries [rad]
    real, dimension(:, :), intent(inout) :: flux
    integer :: i, j, n
    real :: lat, lon, lon_e, x, arg
    logical :: inside
    logical, dimension(size(rows)) :: use
    real, dimension(size(rows)) :: amp, lo, hi, lat0, width, k, ref, ph0

    do n = 1, size(rows)
      use(n) = on(rows(n)%gate(1)) .and. on(rows(n)%gate(2))
      amp(n) = lin_value(rows(n)%amp, kn)
      lo(n) = (rows(n)%lon%lo + kn(rows(n)%lon%ref))*pi/180.
      hi(n) = (rows(n)%lon%hi + kn(rows(n)%lon%ref))*pi/180.
      lat0(n) = rows(n)%lat0 + kn(rows(n)%lat0_ref)
      width(n) = rows(n)%width + kn(rows(n)%width_ref)
      k(n) = rows(n)%k + kn(rows(n)%k_ref)
      ref(n) = kn(rows(n)%phase_ref)*pi/180.
      ph0(n) = rows(n)%phase*pi/180.
    end do
    do j = 1, size(latb) - 1
      lat = 0.5*(latb(j + 1) + latb(j))*180./pi
      do i = 1, size(lonb) - 1
        lon = 0.5*(lonb(i + 1) + lonb(i))
        do n = 1, size(rows)
          if (.not. use(n)) cycle
          if (lat < rows(n)%lat_lo .or. lat > rows(n)%lat_hi) cycle
          x = (lat - lat0(n))/width(n)
          if (rows(n)%band .and. abs(x) > 1.) cycle
          call lon_window(rows(n)%lon, lo(n), hi(n), lon, inside, lon_e)
          if (.not. inside) cycle
          arg = k(n)*(lon_e - ref(n) - ph0(n))
          if (rows(n)%fn == SIN_WAVE) then
            flux(i, j) = flux(i, j) + (1.-x**rows(n)%p)*amp(n)*sin(arg)
          else
            flux(i, j) = flux(i, j) + (1.-x**rows(n)%p)*amp(n)*cos(arg)
          end if
        end do
      end do
    end do
  end subroutine add_waves

  !> Adds the Gaussian terms `rows` to `flux`, for the knobs `kn` and the gates `on`.
  subroutine add_gaussians(rows, kn, on, lonb, latb, flux)
    type(gauss_t), intent(in) :: rows(:)
    real, intent(in) :: kn(0:)
    logical, intent(in) :: on(0:)
    real, dimension(:), intent(in) :: lonb, latb  ! cell boundaries [rad]
    real, dimension(:, :), intent(inout) :: flux
    integer :: i, j, n
    real :: lat, lon, lon_e, lon_d
    logical :: inside
    logical, dimension(size(rows)) :: use
    real, dimension(size(rows)) :: amp, lo, hi, c1, c2

    do n = 1, size(rows)
      use(n) = on(rows(n)%gate(1)) .and. on(rows(n)%gate(2))
      amp(n) = lin_value(rows(n)%amp, kn)
      if (rows(n)%scale > 0) amp(n) = amp(n)*kn(rows(n)%scale)
      lo(n) = (rows(n)%lon%lo + kn(rows(n)%lon%ref))*pi/180.
      hi(n) = (rows(n)%lon%hi + kn(rows(n)%lon%ref))*pi/180.
      c1(n) = rows(n)%f1%c + rows(n)%f1%c_coef*kn(rows(n)%f1%c_ref)
      c2(n) = rows(n)%f2%c + rows(n)%f2%c_coef*kn(rows(n)%f2%c_ref)
    end do
    do j = 1, size(latb) - 1
      lat = 0.5*(latb(j + 1) + latb(j))*180./pi
      do i = 1, size(lonb) - 1
        lon = 0.5*(lonb(i + 1) + lonb(i))
        lon_d = lon*180./pi
        do n = 1, size(rows)
          if (.not. use(n)) cycle
          if (lat < rows(n)%lat_lo .or. lat > rows(n)%lat_hi) cycle
          call lon_window(rows(n)%lon, lo(n), hi(n), lon, inside, lon_e)
          if (.not. inside) cycle
          flux(i, j) = flux(i, j) + amp(n)*gfac(rows(n)%f1, c1(n), lon_d, lat) &
                       *gfac(rows(n)%f2, c2(n), lon_d, lat)
        end do
      end do
    end do
  end subroutine add_gaussians

  !> The Gaussian factor `f` with centre `c`, at `lon`, `lat` [deg] (1 if `f%s` = 0).
  pure real function gfac(f, c, lon, lat)
    type(gfac_t), intent(in) :: f
    real, intent(in) :: c, lon, lat
    gfac = 1.
    if (f%s > 0.) gfac = exp(-(f%a*lon + f%b*lat - c)**2/(2.*f%s**2))
  end function gfac

!########################################################

end module ocean_qflux_mod
