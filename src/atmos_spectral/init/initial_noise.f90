!> Random perturbations of the initial state of the spectral dynamical core.
!>
!> Adds a random perturbation to the spectral temperature, and optionally to the vorticity,
!> when the model starts, to break the symmetry of an idealized initial state or to make
!> the members of an ensemble. The perturbation is band-limited to the total wavenumbers
!> `l_min` to `l_max` (never the zonal mean with the default `m_min = 1`, never the
!> global mean, never the last total wavenumber of the truncation) and rescaled so that
!> its global grid-point rms on each level is exactly the amplitude, times a vertical
!> taper:
!>
!> * `'white'`: every spherical harmonic in the band has the same variance and a random
!>   phase, which makes an isotropic field;
!> * `'spectral'`: an isotropic field whose variance per total wavenumber L follows
!>   `L**(-spectral_slope)` or a Gaussian around `l_peak`;
!> * `'smooth_grid'`: smooth gradient (Perlin) noise with several octaves, evaluated in
!>   three dimensions on the unit sphere (so it has no seam and no polar artefacts),
!>   optionally localized in a Gaussian blob, then transformed to spectral space and
!>   band-limited.
!>
!> The levels are independent, correlated in log pressure, or share one horizontal
!> pattern (`vertical_structure`). The random numbers are a hash of the seed, the ensemble
!> member and the global indices of each coefficient (or lattice point), not the intrinsic
!> random number generator, so a namelist gives the same perturbation for any number of
!> processors and with any compiler. A scalar `initial_noise_applied` in
!> `spectral_dynamics.res.nc` records that a run was perturbed, so that the perturbation
!> is not added again at a restart unless `apply = 'every_start'`.
!>
!> The options are described in
!> [Adding noise to the initial conditions](https://eddy-stanford.github.io/MiMA/Configurations/#adding-noise-to-the-initial-conditions).
!>
!> Namelist: `initial_noise_nml`
!> ([namelist reference](https://eddy-stanford.github.io/MiMA/Parameters/#initial_noise_nml)).
module initial_noise_mod

  use, intrinsic :: iso_fortran_env, only: int64

  use fms_mod, only: mpp_pe, mpp_root_pe, error_mesg, FATAL, NOTE, stdlog, &
                     input_nml_file, check_nml_error, lowercase

  use mpp_mod, only: mpp_sum

  use constants_mod, only: pi

  use transforms_mod, only: get_grid_domain, get_spec_domain, get_num_fourier, get_num_spherical, &
                            get_fourier_inc, get_spherical_wave, get_lon_max, get_deg_lon, &
                            get_deg_lat, trans_grid_to_spherical

  implicit none
  private

  public :: initial_noise_init, add_initial_noise

!===============================================================================================
! namelist variables

  character(len=16) :: noise_type = 'none'
  !! `'none'`: no perturbation; `'white'`: every spherical harmonic in the band has the same
  !! variance; `'spectral'`: variance per total wavenumber set by `spectral_shape`;
  !! `'smooth_grid'`: smooth gradient (Perlin) noise on the sphere, band-limited
  real :: amplitude_t = 0.5
  !! [K] global grid-point rms of the temperature perturbation on each level (times the
  !! taper of `p_top`, `p_bottom` and `taper_depth`)
  logical :: perturb_vorticity = .false. !! also perturb the vorticity (the winds stay non-divergent)
  real :: amplitude_vor = 1.e-6 !! [1/s] global grid-point rms of the vorticity perturbation on each level
  integer :: seed = 1 !! seed of the random numbers (>= 0); with `ensemble_member` it fixes the realisation
  integer :: ensemble_member = 0 !! ensemble member (>= 0): each member gets an independent perturbation
  integer :: l_min = 1 !! lowest total wavenumber perturbed (>= 1: the global mean is never perturbed)
  integer :: l_max = 21
  !! highest total wavenumber perturbed; the last total wavenumber of the truncation is
  !! never perturbed
  integer :: m_min = 1 !! lowest zonal wavenumber perturbed; 1 keeps the zonal mean, 0 also perturbs it
  character(len=16) :: spectral_shape = 'power'
  !! `noise_type = 'spectral'`: the variance per total wavenumber L, `'power'`
  !! (`L**(-spectral_slope)`) or `'gaussian'` (centred on `l_peak`, width `l_width`)
  real :: spectral_slope = 3.
  !! `spectral_shape = 'power'`: exponent p of the variance per total wavenumber `L**(-p)`;
  !! 0 is a flat spectrum in L, 3 is like the enstrophy-cascade range
  real :: l_peak = 10. !! `spectral_shape = 'gaussian'`: total wavenumber of the peak of the variance
  real :: l_width = 4. !! `spectral_shape = 'gaussian'`: width (standard deviation) of the peak, in total wavenumber
  real :: base_wavenumber = 4.
  !! `'smooth_grid'`: approximate total wavenumber of the first (largest) octave of the
  !! gradient noise; each further octave doubles it
  integer :: octaves = 4 !! `'smooth_grid'`: number of octaves of the gradient noise
  real :: gain = 0.5 !! `'smooth_grid'`: amplitude ratio of successive octaves (between 0 and 1)
  real :: lon0 = 0. !! [deg] `'smooth_grid'`: longitude of the centre of the localization
  real :: lat0 = 0. !! [deg] `'smooth_grid'`: latitude of the centre of the localization
  real :: radius = 0.
  !! [deg] `'smooth_grid'`: e-folding great-circle distance of the Gaussian localization
  !! `exp(-d**2/(2 radius**2))` around `lon0`, `lat0`; 0: no localization
  character(len=16) :: vertical_structure = 'correlated'
  !! `'independent'`: an independent pattern on each level; `'correlated'`: correlation
  !! `exp(-|ln p_k - ln p_j|/vertical_scale)` between levels; `'column'`: one pattern on
  !! every level
  real :: vertical_scale = 1. !! [scale heights] `'correlated'`: correlation length in ln p
  real :: p_top = 0. !! [Pa] top of the perturbed layer (reference pressure of the full levels)
  real :: p_bottom = 1.e5 !! [Pa] bottom of the perturbed layer (reference pressure of the full levels)
  real :: taper_depth = 0.
  !! [Pa] depth of the `sin**2` ramp of the amplitude at the top and the bottom of the
  !! perturbed layer; 0: the amplitude is 1 inside the layer and 0 outside
  character(len=16) :: apply = 'cold_start'
  !! when the perturbation is added: `'cold_start'`: only on a cold start; `'once'`: also on
  !! a restart whose restart file does not record an earlier perturbation (to branch an
  !! ensemble from a spun-up state); `'every_start'`: at every start, restarts included

  namelist /initial_noise_nml/ noise_type, amplitude_t, perturb_vorticity, amplitude_vor, &
    seed, ensemble_member, l_min, l_max, m_min, spectral_shape, spectral_slope, l_peak, &
    l_width, base_wavenumber, octaves, gain, lon0, lat0, radius, vertical_structure, &
    vertical_scale, p_top, p_bottom, taper_depth, apply

!===============================================================================================

  ! field number in the hash keys
  integer, parameter :: FIELD_T = 1, FIELD_VOR = 2
  ! kind of random number in the hash keys
  integer, parameter :: TAG_HARMONIC = 1, TAG_LATTICE = 2, TAG_OFFSET = 3
  ! keeps the lattice coordinates of the gradient noise positive in the hash keys
  integer, parameter :: lattice_offset = 65536
  integer(int64), parameter :: mask32 = 4294967295_int64

  logical :: module_is_initialized = .false.
  integer :: ms, me, ns, ne, is, ie, js, je, num_fourier, num_spherical, fourier_inc, num_levels
  logical, allocatable, dimension(:, :) :: in_band  ! (ms:me, ns:ne) the coefficient is perturbed
  real, allocatable, dimension(:, :) :: sigma       ! (ms:me, ns:ne) standard deviation of the harmonic
  real, allocatable, dimension(:) :: taper          ! amplitude factor of each level
  real, allocatable, dimension(:) :: rho            ! correlation with the level above ('correlated')

contains

!===============================================================================================

  !> Reads and checks `initial_noise_nml` and sets up the band, the spectrum and the
  !> vertical structure of the perturbation.
  !>
  !> Must be called after `transforms_init`.
  subroutine initial_noise_init(num_levels_in, triang_trunc, p_ref)
    integer, intent(in) :: num_levels_in  !! number of model levels
    logical, intent(in) :: triang_trunc  !! triangular (`.true.`) or rhomboidal truncation
    real, intent(in), dimension(:) :: p_ref
    !! pressure of the full levels, top to bottom, for the reference surface pressure [Pa]

    integer :: ierr, io, m, n, k, l, num_selected, l_top
    integer, allocatable, dimension(:, :) :: wave
    real :: x, f_max
    character(len=160) :: msg

    if (module_is_initialized) return

    read (input_nml_file, nml=initial_noise_nml, iostat=io)
    ierr = check_nml_error(io, 'initial_noise_nml')
    if (mpp_pe() == mpp_root_pe()) write (stdlog(), nml=initial_noise_nml)

    noise_type = lowercase(noise_type)
    spectral_shape = lowercase(spectral_shape)
    vertical_structure = lowercase(vertical_structure)
    apply = lowercase(apply)

    num_levels = num_levels_in
    if (size(p_ref) /= num_levels) call error_mesg('initial_noise_init', 'size(p_ref) /= num_levels', FATAL)
    call check_noise_nml()

    module_is_initialized = .true.
    if (noise_type == 'none') return

    call get_spec_domain(ms, me, ns, ne)
    call get_grid_domain(is, ie, js, je)
    call get_num_fourier(num_fourier)
    call get_num_spherical(num_spherical)
    call get_fourier_inc(fourier_inc)

    if (noise_type == 'smooth_grid') then
      if (fourier_inc /= 1) call error_mesg('initial_noise_init', "initial_noise_nml: noise_type = 'smooth_grid' "// &
                                            'needs the whole sphere (fourier_inc = 1)', FATAL)
      f_max = base_wavenumber/pi*2.**(octaves - 1)
      if (f_max > 0.1*lattice_offset) call error_mesg('initial_noise_init', 'initial_noise_nml: '// &
                                                      'base_wavenumber*2**(octaves-1) is too large', FATAL)
    end if

    ! the band: l_min <= L <= l_max, m_min <= m*fourier_inc, inside the truncation and
    ! below its last total wavenumber (or its last n and m, for rhomboidal truncation)
    allocate (wave(ms:me, ns:ne), in_band(ms:me, ns:ne), sigma(ms:me, ns:ne))
    call get_spherical_wave(wave)
    num_selected = 0
    do n = ns, ne
      do m = ms, me
        l = wave(m, n)
        in_band(m, n) = l >= l_min .and. l <= l_max .and. m*fourier_inc >= m_min
        if (triang_trunc) then
          in_band(m, n) = in_band(m, n) .and. l <= num_spherical - 2
        else
          in_band(m, n) = in_band(m, n) .and. n <= num_spherical - 2 .and. m <= num_fourier - 1
        end if
        sigma(m, n) = 0.
        if (.not. in_band(m, n)) cycle
        num_selected = num_selected + 1
        ! the variance of each real harmonic is V(L)/(2L+1), V(L) the variance per total wavenumber
        if (noise_type == 'spectral') then
          if (spectral_shape == 'power') then
            sigma(m, n) = sqrt(real(l)**(-spectral_slope)/real(2*l + 1))
          else
            sigma(m, n) = sqrt(exp(-0.5*((real(l) - l_peak)/l_width)**2)/real(2*l + 1))
          end if
        else
          sigma(m, n) = 1.
        end if
      end do
    end do
    deallocate (wave)
    call mpp_sum(num_selected)
    if (num_selected == 0) call error_mesg('initial_noise_init', 'initial_noise_nml: no spherical harmonic '// &
                                           'is in the band set by l_min, l_max and m_min (the last total '// &
                                           'wavenumber of the truncation is never perturbed)', FATAL)

    ! vertical taper and correlation, from the reference pressures
    allocate (taper(num_levels), rho(num_levels))
    do k = 1, num_levels
      if (p_ref(k) < p_top .or. p_ref(k) > p_bottom) then
        taper(k) = 0.
      else if (taper_depth == 0.) then
        taper(k) = 1.
      else
        x = min(1., (p_ref(k) - p_top)/taper_depth, (p_bottom - p_ref(k))/taper_depth)
        taper(k) = sin(0.5*pi*x)**2
      end if
    end do
    if (all(taper == 0.)) call error_mesg('initial_noise_init', 'initial_noise_nml: no level is inside '// &
                                          'the perturbed layer set by p_top, p_bottom and taper_depth', FATAL)
    rho = 0.
    if (vertical_structure == 'correlated') then
      do k = 2, num_levels
        rho(k) = exp(-abs(log(p_ref(k)) - log(p_ref(k - 1)))/vertical_scale)
      end do
    end if

    if (mpp_pe() == mpp_root_pe()) then
      l_top = num_spherical - 2
      if (.not. triang_trunc) l_top = fourier_inc*(num_fourier - 1) + num_spherical - 2
      write (msg, '(a,i0,a,i0,a,i0,a)') 'initial_noise_nml: ', num_selected, &
        ' spherical harmonics per level in the band from L = ', l_min, ' to ', min(l_max, l_top), &
        ', apply = '''//trim(apply)//''''
      call error_mesg('initial_noise_init', trim(msg), NOTE)
    end if

  end subroutine initial_noise_init

!===============================================================================================

  !> Adds the perturbation to both time levels of the spectral temperature (and, with
  !> `perturb_vorticity`, of the vorticity) if `noise_type` and `apply` ask for it.
  !>
  !> The arrays are this processor's part of the spectral fields (`ms:me`, `ns:ne`, all
  !> levels). The increments do not depend on the decomposition.
  subroutine add_initial_noise(cold_start, applied_before, ts, vors, applied, applied_vor, dts, dvors)
    logical, intent(in) :: cold_start  !! there was no restart file
    logical, intent(in) :: applied_before
    !! the restart file records that a perturbation was added in an earlier segment
    complex, intent(inout), dimension(ms:, ns:, :, :) :: ts  !! spectral temperature [K], both time levels
    complex, intent(inout), dimension(ms:, ns:, :, :) :: vors  !! spectral vorticity [1/s], both time levels
    logical, intent(out) :: applied  !! the temperature was perturbed
    logical, intent(out) :: applied_vor  !! the vorticity was perturbed
    complex, intent(out), dimension(ms:, ns:, :) :: dts  !! the increment of the temperature [K] (0 if none)
    complex, intent(out), dimension(ms:, ns:, :) :: dvors  !! the increment of the vorticity [1/s] (0 if none)

    integer :: nt

    if (.not. module_is_initialized) call error_mesg('add_initial_noise', 'initial_noise_init has not been called', FATAL)

    applied = .false.
    applied_vor = .false.
    dts = cmplx(0., 0.)
    dvors = cmplx(0., 0.)
    if (noise_type == 'none') return

    select case (apply)
    case ('cold_start')
      applied = cold_start
    case ('once')
      applied = cold_start .or. .not. applied_before
    case default
      applied = .true.
    end select
    if (.not. applied) then
      if (mpp_pe() == mpp_root_pe()) call error_mesg('add_initial_noise', 'initial_noise_nml: no perturbation '// &
                                                      'added at this restart (apply = '''//trim(apply)//''')', NOTE)
      return
    end if

    call make_perturbation(FIELD_T, amplitude_t, dts)
    do nt = 1, size(ts, 4)
      ts(:, :, :, nt) = ts(:, :, :, nt) + dts
    end do
    if (perturb_vorticity) then
      call make_perturbation(FIELD_VOR, amplitude_vor, dvors)
      do nt = 1, size(vors, 4)
        vors(:, :, :, nt) = vors(:, :, :, nt) + dvors
      end do
      applied_vor = .true.
    end if

  end subroutine add_initial_noise

!===============================================================================================

  ! Checks the namelist values.
  subroutine check_noise_nml()

    select case (noise_type)
    case ('none', 'white', 'spectral', 'smooth_grid')
    case default
      call error_mesg('initial_noise_init', 'initial_noise_nml: noise_type = '''//trim(noise_type)// &
                      ''' is not one of ''none'', ''white'', ''spectral'', ''smooth_grid''', FATAL)
    end select
    select case (spectral_shape)
    case ('power', 'gaussian')
    case default
      call error_mesg('initial_noise_init', 'initial_noise_nml: spectral_shape = '''//trim(spectral_shape)// &
                      ''' is not one of ''power'', ''gaussian''', FATAL)
    end select
    select case (vertical_structure)
    case ('independent', 'correlated', 'column')
    case default
      call error_mesg('initial_noise_init', 'initial_noise_nml: vertical_structure = '''// &
                      trim(vertical_structure)//''' is not one of ''independent'', ''correlated'', ''column''', FATAL)
    end select
    select case (apply)
    case ('cold_start', 'once', 'every_start')
    case default
      call error_mesg('initial_noise_init', 'initial_noise_nml: apply = '''//trim(apply)// &
                      ''' is not one of ''cold_start'', ''once'', ''every_start''', FATAL)
    end select
    if (amplitude_t < 0. .or. amplitude_vor < 0.) &
      call error_mesg('initial_noise_init', 'initial_noise_nml: amplitude_t and amplitude_vor must be >= 0', FATAL)
    if (seed < 0 .or. ensemble_member < 0) &
      call error_mesg('initial_noise_init', 'initial_noise_nml: seed and ensemble_member must be >= 0', FATAL)
    if (l_min < 1 .or. l_max < l_min) &
      call error_mesg('initial_noise_init', 'initial_noise_nml: need 1 <= l_min <= l_max', FATAL)
    if (m_min < 0) call error_mesg('initial_noise_init', 'initial_noise_nml: m_min must be >= 0', FATAL)
    if (l_width <= 0.) call error_mesg('initial_noise_init', 'initial_noise_nml: l_width must be > 0', FATAL)
    if (base_wavenumber <= 0. .or. octaves < 1 .or. gain <= 0. .or. gain > 1.) &
      call error_mesg('initial_noise_init', 'initial_noise_nml: need base_wavenumber > 0, octaves >= 1 '// &
                      'and 0 < gain <= 1', FATAL)
    if (radius < 0.) call error_mesg('initial_noise_init', 'initial_noise_nml: radius must be >= 0', FATAL)
    if (radius > 0. .and. noise_type /= 'smooth_grid' .and. noise_type /= 'none') &
      call error_mesg('initial_noise_init', 'initial_noise_nml: the localization (radius > 0) needs '// &
                      'noise_type = ''smooth_grid''', FATAL)
    if (vertical_scale <= 0.) call error_mesg('initial_noise_init', 'initial_noise_nml: vertical_scale must be > 0', FATAL)
    if (p_top < 0. .or. p_bottom <= p_top .or. taper_depth < 0. .or. 2.*taper_depth > p_bottom - p_top) &
      call error_mesg('initial_noise_init', 'initial_noise_nml: need 0 <= p_top < p_bottom and '// &
                      '0 <= 2*taper_depth <= p_bottom - p_top', FATAL)

  end subroutine check_noise_nml

!===============================================================================================

  ! The perturbation of one field: independent horizontal patterns, combined in the
  ! vertical, then rescaled to the rms amplitude*taper(k) on each level.
  subroutine make_perturbation(field, amplitude, c)
    integer, intent(in) :: field
    real, intent(in) :: amplitude
    complex, intent(out), dimension(ms:, ns:, :) :: c

    complex, allocatable, dimension(:, :, :) :: xi
    real, dimension((num_fourier + 1)*num_levels) :: partial
    real :: mean_square, weight
    integer :: k, m, n, num_patterns, i

    num_patterns = num_levels
    if (vertical_structure == 'column') num_patterns = 1
    allocate (xi(ms:me, ns:ne, num_patterns))
    if (noise_type == 'smooth_grid') then
      call smooth_patterns(field, xi)
    else
      call random_patterns(field, xi)
    end if

    select case (vertical_structure)
    case ('column')
      do k = 1, num_levels
        c(:, :, k) = xi(:, :, 1)
      end do
    case ('independent')
      c = xi
    case ('correlated')
      ! a first-order autoregression down the column: its covariance is exactly
      ! exp(-|ln p_k - ln p_j|/vertical_scale), and the recursion is the product with
      ! the Cholesky factor of that matrix
      c(:, :, 1) = xi(:, :, 1)
      do k = 2, num_levels
        c(:, :, k) = rho(k)*c(:, :, k - 1) + sqrt(1. - rho(k)**2)*xi(:, :, k)
      end do
    end select
    deallocate (xi)

    ! global mean square of each level: (1/2) sum c**2 for m = 0, sum |c|**2 for m > 0
    ! (with the normalization of the Legendre functions); the partial sums of each m come
    ! from one processor and are added in the same order on all
    partial = 0.
    do k = 1, num_levels
      do m = ms, me
        weight = 1.
        if (m == 0) weight = 0.5
        i = m + 1 + (k - 1)*(num_fourier + 1)
        do n = ns, ne
          partial(i) = partial(i) + weight*(real(c(m, n, k))**2 + aimag(c(m, n, k))**2)
        end do
      end do
    end do
    call mpp_sum(partial, size(partial))
    do k = 1, num_levels
      mean_square = 0.
      do m = 0, num_fourier
        mean_square = mean_square + partial(m + 1 + (k - 1)*(num_fourier + 1))
      end do
      if (taper(k) == 0.) then
        c(:, :, k) = cmplx(0., 0.)
      else if (mean_square > 0.) then
        c(:, :, k) = c(:, :, k)*(amplitude*taper(k)/sqrt(mean_square))
      else
        call error_mesg('add_initial_noise', 'initial_noise_nml: the perturbation is zero in the band '// &
                        '(is the Gaussian spectrum of l_peak, l_width far outside l_min to l_max?)', FATAL)
      end if
    end do

  end subroutine make_perturbation

!===============================================================================================

  ! Gaussian random harmonics with the standard deviation sigma(m, n) of each real
  ! harmonic: real for m = 0 (with the factor sqrt(2) of the normalization), random
  ! real and imaginary parts for m > 0. Keyed on the global m, n and the pattern number.
  subroutine random_patterns(field, xi)
    integer, intent(in) :: field
    complex, intent(out), dimension(ms:, ns:, :) :: xi

    integer :: j, m, n, pattern
    integer(int64) :: h

    xi = cmplx(0., 0.)
    do j = 1, size(xi, 3)
      pattern = j
      if (size(xi, 3) == 1) pattern = 0
      do n = ns, ne
        do m = ms, me
          if (.not. in_band(m, n)) cycle
          h = hash_prefix(field, TAG_HARMONIC, pattern)
          h = mix(mix(h, m), n)
          if (m == 0) then
            xi(m, n, j) = cmplx(sqrt(2.)*sigma(m, n)*gaussian(mix(h, 0)), 0.)
          else
            xi(m, n, j) = sigma(m, n)*cmplx(gaussian(mix(h, 0)), gaussian(mix(h, 1)))
          end if
        end do
      end do
    end do

  end subroutine random_patterns

!===============================================================================================

  ! Smooth gradient noise (fractional Brownian motion of 3-D Perlin noise) evaluated at the
  ! grid points on the unit sphere, optionally localized, transformed to spectral space
  ! and limited to the band.
  subroutine smooth_patterns(field, xi)
    integer, intent(in) :: field
    complex, intent(out), dimension(ms:, ns:, :) :: xi

    real, allocatable, dimension(:, :, :) :: grid
    real, allocatable, dimension(:) :: deg_lon, deg_lat
    real, allocatable, dimension(:, :) :: offset
    integer(int64), allocatable, dimension(:) :: lattice
    real :: lon, lat, x(3), d, cos_d, frequency, amp, total, blob
    integer :: lon_max, i, j, jj, o, pattern, axis

    call get_lon_max(lon_max)
    allocate (deg_lon(lon_max), deg_lat(je - js + 1))
    call get_deg_lon(deg_lon)
    call get_deg_lat(deg_lat)
    allocate (grid(is:ie, js:je, size(xi, 3)))
    allocate (offset(3, 0:octaves - 1), lattice(0:octaves - 1))

    do j = 1, size(xi, 3)
      pattern = j
      if (size(xi, 3) == 1) pattern = 0
      ! a random shift of the lattice of each octave and pattern, so that their zeros
      ! (at the lattice points) do not coincide
      do o = 0, octaves - 1
        lattice(o) = mix(hash_prefix(field, TAG_LATTICE, pattern), o)
        do axis = 1, 3
          offset(axis, o) = uniform(mix(mix(hash_prefix(field, TAG_OFFSET, pattern), o), axis))
        end do
      end do
      do jj = js, je
        lat = deg_lat(jj - js + 1)*pi/180.
        do i = is, ie
          lon = deg_lon(i)*pi/180.
          x = (/cos(lat)*cos(lon), cos(lat)*sin(lon), sin(lat)/)
          total = 0.
          amp = 1.
          frequency = base_wavenumber/pi
          do o = 0, octaves - 1
            total = total + amp*perlin3(lattice(o), frequency*x + offset(:, o))
            amp = amp*gain
            frequency = 2.*frequency
          end do
          if (radius > 0.) then
            cos_d = sin(lat)*sin(lat0*pi/180.) + cos(lat)*cos(lat0*pi/180.)*cos(lon - lon0*pi/180.)
            d = acos(max(-1., min(1., cos_d)))
            blob = exp(-0.5*(d/(radius*pi/180.))**2)
            total = total*blob
          end if
          grid(i, jj, j) = total
        end do
      end do
    end do

    call trans_grid_to_spherical(grid, xi)
    do j = 1, size(xi, 3)
      where (.not. in_band) xi(:, :, j) = cmplx(0., 0.)
    end do
    if (ms == 0) xi(0, :, :) = cmplx(real(xi(0, :, :)), 0.)

  end subroutine smooth_patterns

!===============================================================================================

  ! Improved Perlin noise (Perlin 2002) at the point p, with the gradients of the lattice
  ! points chosen by the hash h of the lattice; about -1 to 1, zero at the lattice points.
  real function perlin3(h, p)
    integer(int64), intent(in) :: h
    real, intent(in) :: p(3)

    integer :: i0(3), dx, dy, dz
    real :: f(3), w(3), corner(0:1, 0:1, 0:1)
    integer(int64) :: hc

    i0 = floor(p)
    f = p - i0
    w = f*f*f*(f*(f*6. - 15.) + 10.)
    do dz = 0, 1
      do dy = 0, 1
        do dx = 0, 1
          hc = mix(mix(mix(h, i0(1) + dx + lattice_offset), i0(2) + dy + lattice_offset), i0(3) + dz + lattice_offset)
          corner(dx, dy, dz) = gradient_dot(hc, f(1) - dx, f(2) - dy, f(3) - dz)
        end do
      end do
    end do
    corner(0, :, :) = corner(0, :, :) + w(1)*(corner(1, :, :) - corner(0, :, :))
    corner(0, 0, :) = corner(0, 0, :) + w(2)*(corner(0, 1, :) - corner(0, 0, :))
    perlin3 = corner(0, 0, 0) + w(3)*(corner(0, 0, 1) - corner(0, 0, 0))

  end function perlin3

  ! The dot product of (x, y, z) with one of the 12 gradients (+-1, +-1, 0), ... of Perlin
  ! (2002), chosen by the low four bits of h.
  pure real function gradient_dot(h, x, y, z)
    integer(int64), intent(in) :: h
    real, intent(in) :: x, y, z

    integer :: b
    real :: u, v

    b = int(iand(h, 15_int64))
    if (b < 8) then
      u = x
    else
      u = y
    end if
    if (b < 4) then
      v = y
    else if (b == 12 .or. b == 14) then
      v = x
    else
      v = z
    end if
    if (iand(b, 1) /= 0) u = -u
    if (iand(b, 2) /= 0) v = -v
    gradient_dot = u + v

  end function gradient_dot

!===============================================================================================
! Counter-based random numbers: a 32-bit hash of the keys, so each number depends only on
! its keys (and not on the order of the calls, the decomposition or the compiler).

  ! The hash of the seed, the ensemble member, the field, the kind of number and the pattern.
  pure integer(int64) function hash_prefix(field, tag, pattern)
    integer, intent(in) :: field, tag, pattern

    hash_prefix = mix(mix(mix(mix(mix(2654435769_int64, seed), ensemble_member), field), tag), pattern)

  end function hash_prefix

  ! Adds the key k (>= 0) to the hash h.
  pure integer(int64) function mix(h, k)
    integer(int64), intent(in) :: h
    integer, intent(in) :: k

    mix = lowbias32(ieor(h, iand(int(k, int64), mask32)))

  end function mix

  ! The "lowbias32" integer hash of C. Wellons: a bijection of the 32-bit integers with
  ! a low bias, here on integers 0 to 2**32-1 held in int64.
  pure integer(int64) function lowbias32(x_in)
    integer(int64), intent(in) :: x_in

    integer(int64) :: x

    x = iand(x_in, mask32)
    x = ieor(x, shiftr(x, 16))
    x = mul32(x, 2146121005_int64)  ! 0x7feb352d
    x = ieor(x, shiftr(x, 15))
    x = mul32(x, 2221713035_int64)  ! 0x846ca68b
    x = ieor(x, shiftr(x, 16))
    lowbias32 = x

  end function lowbias32

  ! a*b modulo 2**32 for 0 <= a, b < 2**32, by 16-bit halves of b so that no intermediate
  ! exceeds 2**49 (signed integer overflow is undefined in Fortran)
  pure integer(int64) function mul32(a, b)
    integer(int64), intent(in) :: a, b

    mul32 = iand(a*iand(b, 65535_int64) + shiftl(iand(a*shiftr(b, 16), 65535_int64), 16), mask32)

  end function mul32

  ! A uniform number in [0, 1) with 53 random bits, from the hash h and a further hash of it.
  pure real function uniform(h)
    integer(int64), intent(in) :: h

    integer(int64) :: h2

    h2 = lowbias32(ieor(h, 1779033703_int64))
    uniform = (real(shiftr(h, 5))*67108864. + real(shiftr(h2, 6)))/9007199254740992.

  end function uniform

  ! A standard normal number from the hash h (Box-Muller, with two uniform numbers).
  pure real function gaussian(h)
    integer(int64), intent(in) :: h

    real :: u1, u2

    u1 = uniform(mix(h, 0))
    u2 = uniform(mix(h, 1))
    gaussian = sqrt(-2.*log(1. - u1))*cos(2.*pi*u2)

  end function gaussian

end module initial_noise_mod
