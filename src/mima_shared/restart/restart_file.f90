module restart_file_mod

!-----------------------------------------------------------------------
! Restart files in the layout written by the old fms_io. Each field is
! stored as (x, y, z, Time), scalars and 1-D arrays included, with the
! dimensions xaxis_N, yaxis_N and zaxis_N numbered in order of first use.
! Fields on the file's domain are domain decomposed (written through its
! io_domain); fields on another domain, such as spectral coefficients,
! are gathered and stored whole. Output is buffered and written by
! close_restart, so that the file is fully defined before any data.
!-----------------------------------------------------------------------

  use fms2_io_mod, only: FmsNetcdfFile_t, FmsNetcdfDomainFile_t, open_file, close_file, register_axis, &
                         register_field, write_data, read_data, unlimited, file_exists, &
                         get_variable_size, get_variable_num_dimensions, &
                         get_variable_dimension_names, is_dimension_registered, &
                         get_global_io_domain_indices, variable_exists, &
                         variable_att_exists, get_variable_attribute
  use mpp_mod, only: mpp_error, FATAL, mpp_sum, mpp_npes
  use mpp_domains_mod, only: domain2d, mpp_global_field, mpp_get_compute_domain, &
                             mpp_get_global_domain, mpp_get_io_domain_layout
  use netcdf, only: nf90_set_fill, NF90_NOFILL, NF90_NOERR, NF90_MAX_NAME

  implicit none
  private

  public :: restart_file_type, open_restart_read, open_restart_write, close_restart
  public :: read_restart_field, write_restart_field, get_restart_field_size
  public :: restart_field_exists
  public :: check_field_size

  type restart_field_type
    character(len=64) :: name
    integer :: level                ! record of the Time dimension
    logical :: decomposed           ! x and y are on the file's domain
    integer :: axis(3)              ! x, y and z axis numbers
    real, allocatable :: data(:, :, :)
  end type restart_field_type

  type restart_file_type
    private
    type(FmsNetcdfDomainFile_t) :: fileobj
    logical :: writing = .false.
    integer :: num_fields = 0
    type(restart_field_type), allocatable :: fields(:)
  end type restart_file_type

  interface read_restart_field
    module procedure read_restart_0d, read_restart_1d, read_restart_2d, read_restart_3d
  end interface

  interface write_restart_field
    module procedure write_restart_0d, write_restart_1d, write_restart_2d, write_restart_3d
  end interface

  character(len=*), parameter :: axis_prefix(3) = (/'xaxis_', 'yaxis_', 'zaxis_'/)

contains

!#######################################################################
! Open path (e.g. 'INPUT/atmosphere.res.nc') for reading. Returns
! .false. if it does not exist. It is an error if only a native-format
! restart (the same name without '.nc') exists, if the file is split
! into pieces (path.0000, ...) that the io_layout does not read, or if
! the file cannot be opened on every PE.

  logical function open_restart_read(rst, path, domain)
    type(restart_file_type), intent(inout) :: rst
    character(len=*), intent(in)    :: path
    type(domain2d), intent(in)    :: domain
    integer :: io_layout(2), num_pieces, num_open
    character(len=16) :: ch1, ch2

    if (.not. file_exists(path)) then
      num_pieces = 0
      do while (file_exists(piece_name(path, num_pieces)))
        num_pieces = num_pieces + 1
      end do
      io_layout = mpp_get_io_domain_layout(domain)
      if (num_pieces > 0 .and. (num_pieces /= product(io_layout) .or. product(io_layout) == 1)) then
        write (ch1, '(i0)') num_pieces
        write (ch2, '(i0,",",i0)') io_layout
        call mpp_error(FATAL, 'restart_file_mod: '//trim(path)//' is split into '//trim(ch1)// &
                       ' files ('//trim(piece_name(path, 0))//', ...), which io_layout = '//trim(ch2)// &
                       ' cannot read. Combine them with mppnccombine, or set spec_mpp_nml io_layout'// &
                       ' to the layout that wrote them.')
      end if
    end if

    open_restart_read = open_file(rst%fileobj, path, 'read', domain)
    num_open = merge(1, 0, open_restart_read)
    call mpp_sum(num_open)
    if (num_open /= 0 .and. num_open /= mpp_npes()) then
      write (ch1, '(i0,"/",i0)') num_open, mpp_npes()
      call mpp_error(FATAL, 'restart_file_mod: '//trim(path)//' could be opened on only '// &
                     trim(ch1)//' PEs; the restart files do not match the io_layout.')
    end if
    if (.not. open_restart_read) then
      if (file_exists(path(1:len_trim(path) - 3))) call mpp_error(FATAL, 'restart_file_mod: '// &
                                                                  path(1:len_trim(path) - 3)// &
                                                                  ' is a native-format restart file, which is no longer '// &
                                                                  'supported. Restart from netCDF files ('//trim(path)//').')
      return
    end if
    rst%writing = .false.

  end function open_restart_read

!#######################################################################

  subroutine open_restart_write(rst, path, domain)
    type(restart_file_type), intent(inout) :: rst
    character(len=*), intent(in)    :: path
    type(domain2d), intent(in)    :: domain
    integer :: old_mode

    if (.not. open_file(rst%fileobj, path, 'overwrite', domain)) &
      call mpp_error(FATAL, 'restart_file_mod: cannot create '//trim(path))
! Records not written for a field hold zeros, as fms_io wrote, not fill values.
    if (rst%fileobj%is_root) then
      if (nf90_set_fill(rst%fileobj%ncid, NF90_NOFILL, old_mode) /= NF90_NOERR) &
        call mpp_error(FATAL, 'restart_file_mod: nf90_set_fill failed for '//trim(path))
    end if
    rst%writing = .true.
    rst%num_fields = 0
    allocate (rst%fields(16))

  end subroutine open_restart_write

!#######################################################################

  subroutine close_restart(rst)
    type(restart_file_type), intent(inout) :: rst

    if (rst%writing) then
      call write_fields(rst)
      deallocate (rst%fields)
      rst%num_fields = 0
      rst%writing = .false.
    end if
    call close_file(rst%fileobj)

  end subroutine close_restart

!#######################################################################

  subroutine get_restart_field_size(rst, name, siz)
    type(restart_file_type), intent(inout) :: rst
    character(len=*), intent(in)    :: name
    integer, intent(out)   :: siz(4)

    call global_size(rst%fileobj, name, siz)

  end subroutine get_restart_field_size

!#######################################################################

  logical function restart_field_exists(rst, name)
    type(restart_file_type), intent(inout) :: rst
    character(len=*), intent(in)    :: name

    restart_field_exists = variable_exists(rst%fileobj, name)

  end function restart_field_exists

!#######################################################################
! Stop unless variable name in fileobj has the sizes expected. Any
! further dimensions must have size 1, except that with level the last
! one is the Time dimension, which must hold that record. Also used for
! input files that are not restarts.

  subroutine check_field_size(fileobj, name, expected, level)
    class(FmsNetcdfFile_t), intent(inout) :: fileobj
    character(len=*), intent(in)    :: name
    integer, intent(in)    :: expected(:)
    integer, optional, intent(in)    :: level
    integer, allocatable :: siz(:)
    integer :: n
    logical :: ok

    allocate (siz(get_variable_num_dimensions(fileobj, name)))
    call global_size(fileobj, name, siz)
    n = size(expected)
    if (present(level)) then
      ok = size(siz) == n + 1
      if (ok) ok = siz(n + 1) >= level
    else
      ok = size(siz) >= n
      if (ok) ok = all(siz(n + 1:) == 1)
    end if
    if (ok) ok = all(siz(1:n) == expected)
    if (.not. ok) call mpp_error(FATAL, 'restart_file_mod: size mismatch: file '//trim(fileobj%path)// &
                                 ' variable '//trim(name)//' has '//trim(int_list(siz))//' but the model expects '// &
                                 trim(int_list(expected))//trim(record_text(level)))

  end subroutine check_field_size

  function record_text(level)
    integer, optional, intent(in) :: level
    character(len=32) :: record_text

    record_text = ''
    if (present(level)) write (record_text, '(a,i0)') ' and record ', level

  end function record_text

!#######################################################################
! Global size of each dimension. In a file of an io_layout other than
! (1,1) the decomposed axes cover part of the grid; their
! domain_decomposition attribute holds the global extent.

  subroutine global_size(fileobj, name, siz)
    class(FmsNetcdfFile_t), intent(inout) :: fileobj
    character(len=*), intent(in)    :: name
    integer, intent(out)   :: siz(:)
    character(len=NF90_MAX_NAME) :: dim_names(size(siz))
    integer :: decomposition(4), a

    call get_variable_size(fileobj, name, siz)
    call get_variable_dimension_names(fileobj, name, dim_names)
    do a = 1, min(2, size(siz))
      if (variable_exists(fileobj, dim_names(a))) then
        if (variable_att_exists(fileobj, dim_names(a), 'domain_decomposition')) then
          call get_variable_attribute(fileobj, dim_names(a), 'domain_decomposition', decomposition)
          siz(a) = decomposition(2) - decomposition(1) + 1
        end if
      end if
    end do

  end subroutine global_size

!#######################################################################
! Readers. level is the Time record (default 1). domain is the data's
! domain when it is not the file's; the whole field is then read and
! the compute domain part kept.

  subroutine read_restart_0d(rst, name, data, level)
    type(restart_file_type), intent(inout) :: rst
    character(len=*), intent(in)    :: name
    real, intent(out)   :: data
    integer, optional, intent(in)    :: level
    real :: buf(1, 1, 1)

    call read_field(rst, name, buf, record(level), .false.)
    data = buf(1, 1, 1)

  end subroutine read_restart_0d

  subroutine read_restart_1d(rst, name, data, level)
    type(restart_file_type), intent(inout) :: rst
    character(len=*), intent(in)    :: name
    real, intent(out)   :: data(:)
    integer, optional, intent(in)    :: level
    real :: buf(size(data), 1, 1)

    call read_field(rst, name, buf, record(level), .false.)
    data = buf(:, 1, 1)

  end subroutine read_restart_1d

  subroutine read_restart_2d(rst, name, data, level, domain)
    type(restart_file_type), intent(inout) :: rst
    character(len=*), intent(in)    :: name
    real, intent(out)   :: data(:, :)
    integer, optional, intent(in)    :: level
    type(domain2d), optional, intent(in)    :: domain
    real :: buf(size(data, 1), size(data, 2), 1)

    call read_restart_3d(rst, name, buf, level, domain)
    data = buf(:, :, 1)

  end subroutine read_restart_2d

  subroutine read_restart_3d(rst, name, data, level, domain)
    type(restart_file_type), intent(inout) :: rst
    character(len=*), intent(in)    :: name
    real, intent(out)   :: data(:, :, :)
    integer, optional, intent(in)    :: level
    type(domain2d), optional, intent(in)    :: domain
    real, allocatable :: global(:, :, :)
    integer :: is, ie, js, je, gis, gjs, nx, ny

    if (present(domain)) then
      call mpp_get_global_domain(domain, xbegin=gis, ybegin=gjs, xsize=nx, ysize=ny)
      call mpp_get_compute_domain(domain, is, ie, js, je)
      allocate (global(nx, ny, size(data, 3)))
      call read_field(rst, name, global, record(level), .false.)
      data = global(is - gis + 1:ie - gis + 1, js - gjs + 1:je - gjs + 1, :)
    else
      call read_field(rst, name, data, record(level), .true.)
    end if

  end subroutine read_restart_3d

!#######################################################################
! Writers, with the same arguments as the readers.

  subroutine write_restart_0d(rst, name, data, level)
    type(restart_file_type), intent(inout) :: rst
    character(len=*), intent(in)    :: name
    real, intent(in)    :: data
    integer, optional, intent(in)    :: level

    call add_field(rst, name, reshape((/data/), (/1, 1, 1/)), record(level), .false.)

  end subroutine write_restart_0d

  subroutine write_restart_1d(rst, name, data, level)
    type(restart_file_type), intent(inout) :: rst
    character(len=*), intent(in)    :: name
    real, intent(in)    :: data(:)
    integer, optional, intent(in)    :: level

    call add_field(rst, name, reshape(data, (/size(data), 1, 1/)), record(level), .false.)

  end subroutine write_restart_1d

  subroutine write_restart_2d(rst, name, data, level, domain)
    type(restart_file_type), intent(inout) :: rst
    character(len=*), intent(in)    :: name
    real, intent(in)    :: data(:, :)
    integer, optional, intent(in)    :: level
    type(domain2d), optional, intent(in)    :: domain

    call write_restart_3d(rst, name, reshape(data, (/size(data, 1), size(data, 2), 1/)), level, domain)

  end subroutine write_restart_2d

  subroutine write_restart_3d(rst, name, data, level, domain)
    type(restart_file_type), intent(inout) :: rst
    character(len=*), intent(in)    :: name
    real, intent(in)    :: data(:, :, :)
    integer, optional, intent(in)    :: level
    type(domain2d), optional, intent(in)    :: domain
    real, allocatable :: global(:, :, :)
    integer :: nx, ny

    if (present(domain)) then
      call mpp_get_global_domain(domain, xsize=nx, ysize=ny)
      allocate (global(nx, ny, size(data, 3)))
      call mpp_global_field(domain, data, global)
      call add_field(rst, name, global, record(level), .false.)
    else
      call add_field(rst, name, data, record(level), .true.)
    end if

  end subroutine write_restart_3d

!#######################################################################

  integer function record(level)
    integer, optional, intent(in) :: level

    record = 1
    if (present(level)) record = level

  end function record

!#######################################################################

  subroutine read_field(rst, name, data, level, decomposed)
    type(restart_file_type), intent(inout) :: rst
    character(len=*), intent(in)    :: name
    real, intent(inout) :: data(:, :, :)
    integer, intent(in)    :: level
    logical, intent(in)    :: decomposed
    character(len=NF90_MAX_NAME), allocatable :: dim_names(:)
    integer :: nx, ny

    if (decomposed) then
      call mpp_get_global_domain(rst%fileobj%domain, xsize=nx, ysize=ny)
      call check_field_size(rst%fileobj, name, (/nx, ny, size(data, 3)/), level)
      ! register the field's x and y dimensions, by their names in the file
      allocate (dim_names(get_variable_num_dimensions(rst%fileobj, name)))
      call get_variable_dimension_names(rst%fileobj, name, dim_names)
      if (.not. is_dimension_registered(rst%fileobj, dim_names(1))) &
        call register_axis(rst%fileobj, dim_names(1), 'x')
      if (.not. is_dimension_registered(rst%fileobj, dim_names(2))) &
        call register_axis(rst%fileobj, dim_names(2), 'y')
    else
      call check_field_size(rst%fileobj, name, shape(data), level)
    end if
    call read_data(rst%fileobj, name, data, unlim_dim_level=level)

  end subroutine read_field

!#######################################################################

  subroutine add_field(rst, name, data, level, decomposed)
    type(restart_file_type), intent(inout) :: rst
    character(len=*), intent(in)    :: name
    real, intent(in)    :: data(:, :, :)
    integer, intent(in)    :: level
    logical, intent(in)    :: decomposed
    type(restart_field_type), allocatable :: tmp(:)

    if (.not. rst%writing) call mpp_error(FATAL, 'restart_file_mod: '//trim(name)// &
                                          ' written to a restart file not opened for writing')
    if (rst%num_fields == size(rst%fields)) then
      allocate (tmp(2*size(rst%fields)))
      tmp(1:rst%num_fields) = rst%fields
      call move_alloc(tmp, rst%fields)
    end if
    rst%num_fields = rst%num_fields + 1
    rst%fields(rst%num_fields)%name = name
    rst%fields(rst%num_fields)%level = level
    rst%fields(rst%num_fields)%decomposed = decomposed
    rst%fields(rst%num_fields)%data = data

  end subroutine add_field

!#######################################################################
! Define the axes, Time and the fields, then write them.

  subroutine write_fields(rst)
    type(restart_file_type), intent(inout) :: rst
    integer :: key(rst%num_fields, 3), num_axes(3)
    integer :: n, m, a, i, is, ie
    character(len=16) :: dims(4)

! Axis keys: the length, or 0 for a decomposed x or y.
    num_axes = 0
    do n = 1, rst%num_fields
      do a = 1, 3
        m = size(rst%fields(n)%data, a)
        if (a < 3 .and. rst%fields(n)%decomposed) m = 0
        rst%fields(n)%axis(a) = findloc(key(1:num_axes(a), a), m, dim=1)
        if (rst%fields(n)%axis(a) == 0) then
          num_axes(a) = num_axes(a) + 1
          key(num_axes(a), a) = m
          rst%fields(n)%axis(a) = num_axes(a)
        end if
      end do
    end do

    do a = 1, 3
      do i = 1, num_axes(a)
        if (key(i, a) == 0) then
          call register_axis(rst%fileobj, axis_name(a, i), merge('x', 'y', a == 1))
        else
          call register_axis(rst%fileobj, axis_name(a, i), key(i, a))
        end if
        call register_field(rst%fileobj, axis_name(a, i), 'float', (/axis_name(a, i)/))
      end do
    end do
    call register_axis(rst%fileobj, 'Time', unlimited)
    call register_field(rst%fileobj, 'Time', 'double', (/'Time'/))

    do n = 1, rst%num_fields
      if (any(rst%fields(1:n - 1)%name == rst%fields(n)%name)) cycle
      do a = 1, 3
        dims(a) = axis_name(a, rst%fields(n)%axis(a))
      end do
      dims(4) = 'Time'
      call register_field(rst%fileobj, rst%fields(n)%name, 'double', dims)
    end do

! Axis values are the (global) indices.
    do a = 1, 3
      do i = 1, num_axes(a)
        if (key(i, a) == 0) then
          call get_global_io_domain_indices(rst%fileobj, axis_name(a, i), is, ie)
        else
          is = 1; ie = key(i, a)
        end if
        call write_data(rst%fileobj, axis_name(a, i), (/(real(m, kind=4), m=is, ie)/))
      end do
    end do
    call write_data(rst%fileobj, 'Time', (/(real(m), m=1, maxval(rst%fields(1:rst%num_fields)%level))/))

    do n = 1, rst%num_fields
      call write_data(rst%fileobj, rst%fields(n)%name, rst%fields(n)%data, &
                      unlim_dim_level=rst%fields(n)%level)
    end do

  end subroutine write_fields

!#######################################################################

  function piece_name(path, i)
    character(len=*), intent(in) :: path
    integer, intent(in) :: i
    character(len=len_trim(path) + 5) :: piece_name

    write (piece_name, '(a,".",i4.4)') trim(path), i

  end function piece_name

!#######################################################################

  function int_list(n)
    integer, intent(in) :: n(:)
    character(len=13*size(n) + 2) :: int_list

    write (int_list, '("(", *(i0, :, ", "))') n
    int_list = trim(int_list)//')'

  end function int_list

!#######################################################################

  function axis_name(a, i)
    integer, intent(in) :: a, i
    character(len=16) :: axis_name

    write (axis_name, '(a,i0)') axis_prefix(a), i

  end function axis_name

!#######################################################################

end module restart_file_mod
