!> Domain decompositions of the grid and of the spectral coefficients, for the spectral
!> transforms (transpose method).
!>
!> The grid domain is decomposed in latitude (layout `1, npes`; the number of latitudes must
!> be a multiple of the number of processors) and the spectral domain in the zonal index m
!> (layout `npes, 1`), into contiguous blocks of m chosen by `spectral_decomposition`.
!> `io_layout` sets the I/O domains of both.
!>
!> Namelist: `spec_mpp_nml`
!> ([namelist reference](https://eddy-stanford.github.io/MiMA/Parameters/#spec_mpp_nml)).
!>
!> Original authors: V. Balaji.
module spec_mpp_mod
  use fms_mod, only: mpp_pe, mpp_root_pe, mpp_npes, mpp_error, FATAL, &
                     input_nml_file, check_nml_error, stdlog

  use mpp_domains_mod, only: mpp_domains_init, domain1D, domain2D, &
                             mpp_define_domains, mpp_get_compute_domain, mpp_get_compute_domains, &
                             mpp_get_domain_components, mpp_get_pelist, mpp_define_io_domain

  implicit none
  private

  type(domain2D), save, public :: grid_domain, spectral_domain
  !! `grid_domain`: decomposition of the Gaussian grid (layout `1, npes`); `spectral_domain`:
  !! decomposition of the spectral coefficients (m, n) (layout `npes, 1`)
  logical, private :: module_is_initialized = .false.
  integer, private :: pe, npes

  integer, private :: io_layout(2) = (/1, 1/)
  !! I/O layout of the grid domain: each I/O domain (group of processes) writes one file of the
  !! diagnostics and restarts. `1,1` gives single files. Each entry must divide the processor
  !! layout, which is `1, npes`; with more than one I/O domain the files are split
  !! (`atmos_daily.nc.0000`, ...) and must be joined with `mppnccombine`. Split restart files can
  !! only be read with the same `io_layout`.
  ! The spectral domain, which is decomposed along m rather than latitude,
  ! uses the transposed layout, (io_layout(2), io_layout(1)).

  character(len=16), private :: spectral_decomposition = 'balanced'
  !! how the zonal indices m are divided among the processes, in contiguous blocks:
  !! `'equal'`: blocks of equal numbers of m; `'balanced'`: blocks of equal estimated work.
  !! With triangular truncation the Legendre transforms of zonal index m cost about
  !! `num_spherical - m` (the number of its retained coefficients), so equal blocks leave the
  !! first process with the most work. `'balanced'` gives each m the cost
  !! `num_spherical - m + spectral_cost_per_m` and minimizes the largest block cost. It applies
  !! only to triangular truncation with `fourier_inc = 1` (otherwise `'equal'` is used). The
  !! results and restart files are the same as with `'equal'`, bit for bit, and a run may restart
  !! with either.
  real, private :: spectral_cost_per_m = 20.
  !! for `spectral_decomposition = 'balanced'`: the cost of each zonal index m that does not
  !! depend on its number of retained coefficients (the time step and damping work on all
  !! `num_spherical + 1` coefficients, the transposes on all latitudes), in units of the cost of
  !! one retained coefficient. Measured at T42 on 4 processes: about 16 with 1 OpenMP thread and
  !! 28 with 2, as the threads speed up the Legendre transforms more than the rest.

  namelist /spec_mpp_nml/ io_layout, spectral_decomposition, spectral_cost_per_m

  public :: spec_mpp_init, get_grid_domain, get_spec_domain, spec_mpp_end

contains

!=======================================================================================================================

  !> Reads `spec_mpp_nml` and defines the grid and spectral domains and their I/O domains.
  subroutine spec_mpp_init(num_fourier, num_spherical, num_lon, lat_max, fourier_inc, triang_trunc)
    integer, intent(in) ::  num_fourier, num_spherical, num_lon, lat_max
    !! `num_fourier`, `num_spherical`: largest spectral indices m and n; `num_lon`, `lat_max`:
    !! numbers of longitudes and latitudes
    integer, intent(in), optional :: fourier_inc
    !! zonal wavenumber increment (default 1); `spectral_decomposition = 'balanced'` needs 1
    logical, intent(in), optional :: triang_trunc
    !! triangular truncation (default `.true.`); `spectral_decomposition = 'balanced'` needs it
    integer :: io, ierr
    integer :: layout(2)
    integer, allocatable :: xextent(:)
    logical :: balanced
    character(len=4) :: chtmp1, chtmp2

    if (module_is_initialized) return
    call mpp_domains_init()
    pe = mpp_pe()
    npes = mpp_npes()

    read (input_nml_file, nml=spec_mpp_nml, iostat=io)
    ierr = check_nml_error(io, 'spec_mpp_nml')
    if (spectral_decomposition /= 'equal' .and. spectral_decomposition /= 'balanced') &
      call mpp_error(FATAL, 'SPEC_MPP_INIT: spectral_decomposition must be ''equal'' or ''balanced'', not '''// &
                     trim(spectral_decomposition)//'''')
    if (spectral_cost_per_m < 0.) call mpp_error(FATAL, 'SPEC_MPP_INIT: spectral_cost_per_m must not be negative')

    if (pe == mpp_root_pe()) write (stdlog(), nml=spec_mpp_nml)

!grid domain: 1D decomposition along Y (transforms_mod relies on whole latitude rows on each PE)
    layout = (/1, npes/)
    call mpp_define_domains((/1, num_lon, 1, lat_max/), layout, grid_domain)
    if (pe == mpp_root_pe()) call print_decomp(npes, layout, grid_domain, 'ATMOS MODEL DOMAIN DECOMPOSITION')

!requirement of equal domains: can be generalized to retain mirror symmetry between N/S if unequal.
!the equal-domains requirement permits us to eliminate one buffer/unbuffer in the transpose_fourier routines.
    if (mod(lat_max, layout(2)) /= 0) then
      write (chtmp1, '(i4)') layout(2)
      write (chtmp2, '(i4)') lat_max
      call mpp_error(FATAL, 'SPEC_MPP_INIT:Requires num_lat_rows/num_pes=int;num_pes='&
     &//chtmp1//';num_lat_rows='//chtmp2)
    end if
    call mpp_define_io_domain(grid_domain, io_layout)

!spectral domain: 1D decomposition along M (transforms_mod relies on all n on each PE)
    layout = (/npes, 1/)
    balanced = spectral_decomposition == 'balanced' .and. layout(1) > 1 .and. layout(1) <= num_fourier + 1
    if (present(fourier_inc)) balanced = balanced .and. fourier_inc == 1
    if (present(triang_trunc)) balanced = balanced .and. triang_trunc
    if (balanced) then
      allocate (xextent(layout(1)))
      call balanced_m_blocks(num_fourier, num_spherical, xextent)
      call mpp_define_domains((/0, num_fourier, 0, num_spherical/), layout, spectral_domain, xextent=xextent)
    else
      call mpp_define_domains((/0, num_fourier, 0, num_spherical/), layout, spectral_domain)
    end if
    call mpp_define_io_domain(spectral_domain, (/io_layout(2), io_layout(1)/))
    if (pe == mpp_root_pe()) call print_decomp(npes, layout, spectral_domain, 'SPECTRAL DOMAIN DECOMPOSITION (X = M, Y = N)')

    module_is_initialized = .true.
    return
  end subroutine spec_mpp_init
!=======================================================================================================================

  !> Divides the zonal indices 0 to `num_fourier` into `size(xextent)` contiguous, nonempty
  !> blocks that minimize the largest block cost, the cost of m being
  !> `num_spherical - m + spectral_cost_per_m` (triangular truncation).
  subroutine balanced_m_blocks(num_fourier, num_spherical, xextent)
    integer, intent(in) :: num_fourier, num_spherical
    integer, intent(out) :: xextent(:)
    real :: cost(0:num_fourier), lo, hi, mid
    integer :: m, iter
    logical :: ok

    do m = 0, num_fourier
      cost(m) = real(max(num_spherical - m, 0)) + spectral_cost_per_m
    end do
! bisect for the smallest largest block cost that a partition achieves
    lo = maxval(cost)
    hi = sum(cost)
    do iter = 1, 100
      mid = 0.5*(lo + hi)
      if (mid <= lo .or. mid >= hi) exit
      call fill_blocks(mid, ok)
      if (ok) then
        hi = mid
      else
        lo = mid
      end if
    end do
    call fill_blocks(hi, ok)
    if (.not. ok) call mpp_error(FATAL, 'SPEC_MPP_INIT: balanced_m_blocks failed')

  contains

    ! Fills each block in turn with as many m as fit in max_cost, leaving at least one m
    ! for each later block.
    subroutine fill_blocks(max_cost, ok)
      real, intent(in) :: max_cost
      logical, intent(out) :: ok
      real :: block_cost
      integer :: p, m, nb

      nb = size(xextent)
      xextent = 0
      m = 0
      do p = 1, nb
        block_cost = 0.
        do while (m <= num_fourier - (nb - p))
          if (block_cost + cost(m) > max_cost) exit
          block_cost = block_cost + cost(m)
          xextent(p) = xextent(p) + 1
          m = m + 1
        end do
        if (xextent(p) == 0) exit
      end do
      ok = m == num_fourier + 1 .and. all(xextent > 0)

    end subroutine fill_blocks

  end subroutine balanced_m_blocks
!=======================================================================================================================

  subroutine print_decomp(npes, layout, Domain, title)
    integer, intent(in) :: npes, layout(2)
    type(domain2d), intent(in) :: Domain
    character(len=*), intent(in) :: title
    integer, dimension(0:npes - 1) :: xsize, ysize
    integer :: i, j, xlist(layout(1)), ylist(layout(2))
    type(domain1D) :: Xdom, Ydom

    call mpp_get_compute_domains(Domain, xsize=xsize, ysize=ysize)
    call mpp_get_domain_components(Domain, Xdom, Ydom)
    call mpp_get_pelist(Xdom, xlist)
    call mpp_get_pelist(Ydom, ylist)

    write (*, '(a)') title
    write (*, 110) (xsize(xlist(i)), i=1, layout(1))
    write (*, 120) (ysize(ylist(j)), j=1, layout(2))

110 format('  X-AXIS = ', 24i4, /, (11x, 24i4))
120 format('  Y-AXIS = ', 24i4, /, (11x, 24i4))

  end subroutine print_decomp
!=======================================================================================================================

  !> Returns the index bounds of this processor's part of the grid.
  subroutine get_grid_domain(is, ie, js, je)
    integer, intent(out) :: is, ie, js, je  !! `is`, `ie`: first and last longitude index; `js`, `je`: first and last latitude index

    if (.not. module_is_initialized) call mpp_error(FATAL, 'subroutine get_grid_domain: spec_mpp is not initialized')

    call mpp_get_compute_domain(grid_domain, is, ie, js, je)

    return
  end subroutine get_grid_domain
!=======================================================================================================================

  !> Returns the index bounds of this processor's part of the spectral coefficients.
  subroutine get_spec_domain(ms, me, ns, ne)
    integer, intent(out) :: ms, me, ns, ne  !! `ms`, `me`: first and last index m; `ns`, `ne`: first and last index n

    if (.not. module_is_initialized) call mpp_error(FATAL, 'subroutine get_spec_domain: spec_mpp is not initialized')

    call mpp_get_compute_domain(spectral_domain, ms, me, ns, ne)

    return
  end subroutine get_spec_domain
!=======================================================================================================================

  !> Marks the module as not initialized.
  subroutine spec_mpp_end

    module_is_initialized = .false.

    return
  end subroutine spec_mpp_end
!=======================================================================================================================

end module spec_mpp_mod
