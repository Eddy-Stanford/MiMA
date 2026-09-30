!> Domain decompositions of the grid and of the spectral coefficients, for the spectral
!> transforms (transpose method).
!>
!> The grid domain is decomposed in latitude (layout `1, npes`; the number of latitudes must
!> be a multiple of the number of processors) and the spectral domain in the zonal index m
!> (layout `npes, 1`). `io_layout` sets the I/O domains of both.
!>
!> Namelist: `spec_mpp_nml`
!> ([namelist reference](https://eddy-stanford.github.io/MiMA/Parameters/#spec_mpp_nml)).
!>
!> Original authors: V. Balaji.
module spec_mpp_mod
  use fms_mod, only: mpp_pe, mpp_root_pe, mpp_npes, write_version_number, mpp_error, FATAL, &
                     input_nml_file, check_nml_error, stdlog

  use mpp_domains_mod, only: mpp_domains_init, domain1D, domain2D, GLOBAL_DATA_DOMAIN, &
                             mpp_define_domains, mpp_get_compute_domain, mpp_get_compute_domains, &
                             mpp_get_domain_components, mpp_get_pelist, mpp_define_io_domain

  implicit none
  private

  character(len=128), private :: version = '$Id: spec_mpp.f90,v 10.0 2003/10/24 22:01:02 fms Exp $'
  character(len=128), private :: tagname = '$Name: lima $'
  type(domain2D), save, public :: grid_domain, spectral_domain, global_spectral_domain
  !! `grid_domain`: decomposition of the Gaussian grid; `spectral_domain`: decomposition of the
  !! spectral coefficients (m, n); `global_spectral_domain`: as `spectral_domain`, with global data
  !! boundaries
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

  namelist /spec_mpp_nml/ io_layout

  public :: spec_mpp_init, get_grid_domain, get_spec_domain, spec_mpp_end

contains

!=======================================================================================================================

  !> Reads `spec_mpp_nml` and defines the grid and spectral domains and their I/O domains.
  subroutine spec_mpp_init(num_fourier, num_spherical, num_lon, lat_max, grid_layout, spectral_layout)
    integer, intent(in) ::  num_fourier, num_spherical, num_lon, lat_max
    !! `num_fourier`, `num_spherical`: largest spectral indices m and n; `num_lon`, `lat_max`:
    !! numbers of longitudes and latitudes
    integer, intent(in), optional :: grid_layout(2), spectral_layout(2)
    !! `grid_layout`, `spectral_layout`: processor layouts of the grid and spectral domains (default:
    !! `1, npes` and `npes, 1`)
    integer :: i, io, ierr
    integer :: layout(2)
    character(len=4) :: chtmp1, chtmp2

    if (module_is_initialized) return
    call mpp_domains_init()
    pe = mpp_pe()
    npes = mpp_npes()

    read (input_nml_file, nml=spec_mpp_nml, iostat=io)
    ierr = check_nml_error(io, 'spec_mpp_nml')

    call write_version_number(version, tagname)
    if (pe == mpp_root_pe()) write (stdlog(), nml=spec_mpp_nml)

!grid domain: by default, 1D decomposition along Y
    layout = (/1, npes/)
    if (present(grid_layout)) layout = grid_layout
    call mpp_define_domains((/1, num_lon, 1, lat_max/), layout, grid_domain)
    if (pe == mpp_root_pe()) call print_decomp(npes, layout, grid_domain)

!requirement of equal domains: can be generalized to retain mirror symmetry between N/S if unequal.
!the equal-domains requirement permits us to eliminate one buffer/unbuffer in the transpose_fourier routines.
    if (mod(lat_max, layout(2)) .ne. 0) then
!       call mpp_error( FATAL, 'SPEC_MPP_INIT: currently requires equal grid domains on all PEs.' )
      write (chtmp1, '(i4)') layout(2)
      write (chtmp2, '(i4)') lat_max
      call mpp_error(FATAL, 'SPEC_MPP_INIT:Requires num_lat_rows/num_pes=int;num_pes='&
     &//chtmp1//';num_lat_rows='//chtmp2)
    end if
    call mpp_define_io_domain(grid_domain, io_layout)

!spectral domain: by default, 1D decomposition along M
    layout = (/npes, 1/)
    if (present(spectral_layout)) layout = spectral_layout
    call mpp_define_domains((/0, num_fourier, 0, num_spherical/), layout, spectral_domain)
    call mpp_define_io_domain(spectral_domain, (/io_layout(2), io_layout(1)/))

!global spectral domains (may be used for I/O) are the same as spectral domains, with global data boundaries
    call mpp_define_domains((/0, num_fourier, 0, num_spherical/), layout, global_spectral_domain, &
                            xflags=GLOBAL_DATA_DOMAIN, yflags=GLOBAL_DATA_DOMAIN)

    module_is_initialized = .true.
    return
  end subroutine spec_mpp_init
!=======================================================================================================================

  subroutine print_decomp(npes, layout, Domain)
    integer, intent(in) :: npes, layout(2)
    type(domain2d), intent(in) :: Domain
    integer, dimension(0:npes - 1) :: xsize, ysize
    integer :: i, j, xlist(layout(1)), ylist(layout(2))
    type(domain1D) :: Xdom, Ydom

    call mpp_get_compute_domains(Domain, xsize=xsize, ysize=ysize)
    call mpp_get_domain_components(Domain, Xdom, Ydom)
    call mpp_get_pelist(Xdom, xlist)
    call mpp_get_pelist(Ydom, ylist)

    write (*, 100)
    write (*, 110) (xsize(xlist(i)), i=1, layout(1))
    write (*, 120) (ysize(ylist(j)), j=1, layout(2))

100 format('ATMOS MODEL DOMAIN DECOMPOSITION')
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
