!   Author: ZE
!   Date: 5/20/21
module dd_drag_mod
    use forpy_mod
    use iso_fortran_env,  only: real64
    use time_manager_mod, only: time_type 

    implicit none

    public dd_drag_init, dd_drag_calc, dd_drag_end

    logical                             :: module_is_initialized =.false.
    integer                             :: ierror
    type(module_py)                     :: wavenet
    type(list)                          :: paths 

    contains
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!           PUBLIC SUBROUTINES
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine dd_drag_init ()
    if (module_is_initialized) return

    ierror = forpy_initialize ()

    !---------------------------------------------------------------------
    ! Load Python Module: wavenet.py 
    !---------------------------------------------------------------------
    ierror = get_sys_path(paths)
    ierror = paths%append('/scratch/zespinos/models/code/MiMA/src/atmos_param')
    ierror = import_py(wavenet, "wavenet")
    if (ierror/=0) then; call err_print; endif
    ierror = print_py(wavenet)
    if (ierror/=0) then; call err_print; endif
    write(*,*) "############ PRINT WAVENET PACKAGE ############"
    ierror = call_py_noret(wavenet, "init_wavenet")
    write(*,*) "############ INIT WAVENET ############"
    if (ierror/=0) then; call err_print; endif

    module_is_initialized = .true.

end subroutine dd_drag_init

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subroutine dd_drag_calc (is, js, lat, pfull, zfull, temp, uuu, vvv,  &
                         Time, delt, gwfcng_x, gwfcng_y)
    !---------------------------------------------------------------------
    ! Calculate Data-Driven Gravity Wave Drag
    !---------------------------------------------------------------------
    ! Note: May need asynchronous flag in front of elements
    
    !-------------------------------------------------------------------
    !    Arguments (Intent in):
    !
    !       is,js    starting subdomain i,j indices of data in 
    !                the physics_window being integrated
    !       lat      array of model latitudes at cell boundaries [radians]
    !       pfull    pressure at model full levels [ Pa  ]
    !       zfull    height at model full levels [ m  ]
    !       temp     temperature at model levels [ deg K ]
    !       uuu      zonal wind  [ m/s  ]
    !       vvv      meridional wind  [ m/s  ]
    !       Time     current time, needed for diagnostics [ time_type  ]
    !       delt     physics time step [ s  ]
    !
    !    Arguments (Intent out):
    !
    !       gwfcng_x time tendency for u eqn due to gravity-wave forcing
    !                [ m/s^2  ]
    !       gwfcng_y time tendency for v eqn due to gravity-wave forcing
    !                [ m/s^2  ]
    !
    !    Local Variables: 
    !       
    !       paths environment variable PYTHONPATH - add current directory "." to
    !       sys.path
    !-------------------------------------------------------------------
    ! Intent In
    integer,                              intent(in)      :: is, js
    real,               dimension(:,:),   intent(in)      :: lat
    real,               dimension(:,:,:), intent(in)      :: pfull, zfull, temp, uuu, vvv
    type(time_type),                      intent(in)      :: Time
    !real, intent(in) :: Time
    real           ,                      intent(in)      :: delt
    ! Intent Out  
    real,               dimension(:,:,:), intent(out)     :: gwfcng_x, gwfcng_y
    ! Local Variables
    type(ndarray)                                         :: lat_py, pfull_py, zfull_py, temp_py, uuu_py, vvv_py, gwfcng_x_py, gwfcng_y_py, gwfu_arr, gwfv_arr
    type(object)                                          :: gwfu, gwfv
    type(tuple)                                           :: args
    real(kind=real64),  dimension(:,:,:), pointer         :: matrix_u, matrix_v
        
    !---------------------------------------------------------------------
    ! Covert input into Forpy Numpy Arrays 
    !---------------------------------------------------------------------
    !write(*,*) "############ Convert Input to Forpy  ############"
    ierror = ndarray_create(lat_py, lat)
    ierror = ndarray_create(pfull_py, pfull)
    ierror = ndarray_create(zfull_py, zfull)
    ierror = ndarray_create(temp_py, temp)
    ierror = ndarray_create(uuu_py, uuu)
    ierror = ndarray_create(vvv_py, vvv)
    ierror = ndarray_create(gwfcng_x_py, gwfcng_x)
    ierror = ndarray_create(gwfcng_y_py, gwfcng_y)
    if (ierror/=0) then; call err_print; endif

    !---------------------------------------------------------------------
    ! Create Python Argument 
    !---------------------------------------------------------------------
    !write(*,*) "############ Create Python Argument ############"
    ierror = tuple_create(args, 11)
    ierror = args%setitem(0, is)
    ierror = args%setitem(1, js)
    ierror = args%setitem(2, lat_py)
    ierror = args%setitem(3, pfull_py)
    ierror = args%setitem(4, zfull_py)
    ierror = args%setitem(5, temp_py)
    ierror = args%setitem(6, uuu_py)
    ierror = args%setitem(7, vvv_py)
    !ierror = args%setitem(8, Time)
    ierror = args%setitem(8, delt)
    ierror = args%setitem(9, gwfcng_x_py)
    ierror = args%setitem(10, gwfcng_y_py)
    if (ierror/=0) then; call err_print; endif

    !---------------------------------------------------------------------
    ! Calculate GWFD
    !---------------------------------------------------------------------
    !write(*,*) "############ Calculate GWFD  ############"
    ierror = call_py(gwfu, wavenet, "predict_u", args)
    ierror = call_py(gwfv, wavenet, "predict_v", args)
    if (ierror/=0) then; call err_print; endif

    ierror = cast(gwfu_arr, gwfu)
    ierror = cast(gwfv_arr, gwfv)
    if (ierror/=0) then; call err_print; endif
    
    ierror = gwfu_arr%get_data(matrix_u)
    ierror = gwfv_arr%get_data(matrix_v)
    if (ierror/=0) then; call err_print; endif

    gwfcng_x = matrix_u
    gwfcng_y = matrix_v
    if (ierror/=0) then; call err_print; endif

    !---------------------------------------------------------------------
    ! Destroy Objects
    !---------------------------------------------------------------------
    call lat_py%destroy
    call pfull_py%destroy
    call zfull_py%destroy
    call temp_py%destroy
    call uuu_py%destroy
    call vvv_py%destroy
    call gwfcng_x_py%destroy
    call gwfcng_y_py%destroy
    call gwfu_arr%destroy
    call gwfv_arr%destroy
    call gwfu%destroy
    call gwfv%destroy
    call args%destroy

end subroutine dd_drag_calc


subroutine dd_drag_end()
    write(*,*) "############ Finalize Forpy ############"
    call wavenet%destroy
    call paths%destroy

    call forpy_finalize

    module_is_initialized = .false.

end subroutine dd_drag_end 

end module dd_drag_mod
