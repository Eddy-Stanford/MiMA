module dd_drag_mod 
    use forpy_mod
    use iso_fortran_env, only: real64
    !REAL use time_manager_mod, only:  time_type

    implicit none

    logical :: module_is_initialized=.false.
    integer :: ierror
    type(list)                              :: paths
    type(module_py)                         :: wavenet                                           

    public dd_drag_init, dd_drag_calc, dd_drag_end

    contains

    subroutine dd_drag_init ()
        !(lonb, latb, pref, Time, axes)
        !------------------------------------------------------------
        !   INITIALIZE Data-Driven Gravity Wave Drag
        !------------------------------------------------------------

        !------------------------------------------------------------
        !   Arguments (Intent in)
        !       lonb = longitude in radians of the grid box edges
        !       latb = latitude in radians of the grid box edges
        !       pref = array of reference pressures at full levels (plus surface value at nlev+1), 
        !           based on 1012.25hpa pstar [ Pa ]
        !       Time = current time (time_type)
        !       axes = data axes for diagnostics
        !------------------------------------------------------------
        ! real,    dimension(:),  intent(in)   :: lonb, latb, pref
        ! integer, dimension(4),  intent(in)   :: axes
        ! real,                   intent(in)   :: Time
        ! !REAL type(time_type), intent(in) :: Time
        

        if (module_is_initialized) return

        ierror = forpy_initialize ()

        !-------------------------------------------------------------------                         
        !   Load Python Module: wavenet.py
        !-------------------------------------------------------------------                         
        ierror = get_sys_path(paths)
        ierror = paths%append('/home/zespinosa/Stanford/research/Coupling_MIMA_GW/fortran_wavenet_adapter/')
        ierror = import_py(wavenet, 'wavenet_2') 
        if(ierror/=0) then;call err_print; write(*,*) ierror; endif
        ierror = print_py(wavenet)
        ierror = call_py_noret(wavenet, 'init_wavenet')

        module_is_initialized = .true.

    end subroutine dd_drag_init 


    subroutine dd_drag_calc (is, js, lat, pfull, zfull, temp, uuu, vvv, Time, delt, gwfcng_x, gwfcng_y)  
        !------------------------------------------------------------
        !   Calculate Data-Driven Gravity Wave Drag
        !------------------------------------------------------------
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
        !     Arguments (Intent out):                                                                  
        !                                                                                            
        !       gwfcng_x time tendency for u eqn due to gravity-wave forcing                         
        !                [ m/s^2  ]                                                                  
        !       gwfcng_y time tendency for v eqn due to gravity-wave forcing                         
        !                [ m/s^2  ]                                                                  
        ! 
        !     Local Variables: 
        !       paths    environment variable PYTHONPATH - add current directory "." to sys.path
        ! 
        !-------------------------------------------------------------------                         
        integer,                intent(in)      :: is, js 
        real, dimension(:,:),   intent(in)      :: lat                                               
        real, dimension(:,:,:), intent(in)      :: pfull, zfull, temp, uuu, vvv                      
        !REAL type(time_type),        intent(in)      :: Time                                              
        real,                   intent(in)      :: Time                                              
        real,                   intent(in)      :: delt                                              

        type(ndarray)                           :: lat_py, pfull_py, zfull_py, temp_py, uuu_py, vvv_py, &
                                                    gwfu_arr, gwfv_arr
        type(object)                            :: gwfuv, gwfu, gwfv
        type(tuple)                             :: predictions, args
                                                                                                        
        real, dimension(:,:,:), intent(out)     :: gwfcng_x, gwfcng_y
        real(kind=real64), dimension(:,:,:), pointer :: matrix_u, matrix_v, matrix_uv 


        !-------------------------------------------------------------------                         
        !   Convert input to forpy numpy types
        !-------------------------------------------------------------------                         
        ierror = ndarray_create(lat_py, lat)
        ierror = ndarray_create(pfull_py, pfull)
        ierror = ndarray_create(zfull_py, zfull)
        ierror = ndarray_create(temp_py, temp)
        ierror = ndarray_create(uuu_py, uuu)
        ierror = ndarray_create(vvv_py, vvv)


        !-------------------------------------------------------------------                         
        !  Create Python Argument
        !-------------------------------------------------------------------                         
        ierror = tuple_create(args, 10)
        if(ierror/=0) then;call err_print; stop;endif
        ierror = args%setitem(0, is)
        ierror = args%setitem(1, js)
        ierror = args%setitem(2, lat_py)
        ierror = args%setitem(3, pfull_py)
        ierror = args%setitem(4, zfull_py)
        ierror = args%setitem(5, temp_py)
        ierror = args%setitem(6, uuu_py)
        ierror = args%setitem(7, vvv_py)
        ierror = args%setitem(8, Time)
        ierror = args%setitem(9, delt)

        !-------------------------------------------------------------------                         
        !  Calculate GWFD 
        !------ -------------------------------------------------------------                         
        ierror = call_py(gwfuv, wavenet, "predict", args)                                                                                    
        if(ierror/=0) then;call err_print;stop;endif

        ! extra step when using combined
        ierror = cast(predictions, gwfuv)
        ierror = predictions%getitem(gwfu, 0)
        ierror = predictions%getitem(gwfv, 1)
        if(ierror/=0) then;call err_print;stop;endif

        ierror = cast(gwfu_arr, gwfu)
        ierror = cast(gwfv_arr, gwfv)
        if(ierror/=0) then;call err_print;stop;endif
        
        ierror = gwfu_arr%get_data(matrix_u)
        ierror = gwfv_arr%get_data(matrix_v)
        if(ierror/=0) then;call err_print;stop;endif

        gwfcng_x = matrix_u
        gwfcng_y = matrix_v
        if(ierror/=0) then;call err_print;stop;endif

        !-------------------------------------------------------------------                         
        !  Destroy Objects 
        !------ -------------------------------------------------------------                         
        ! Might need to do this
        write(*,*) "#### Delete Objects ####"
        call lat_py%destroy
        call pfull_py%destroy
        call zfull_py%destroy
        call temp_py%destroy
        call uuu_py%destroy
        call vvv_py%destroy
        call gwfu_arr%destroy
        call gwfv_arr%destroy
        call gwfu%destroy
        call gwfv%destroy
        call args%destroy

    end subroutine dd_drag_calc 


    subroutine dd_drag_end ()
        call forpy_finalize
        
        module_is_initialized = .false.

    end subroutine dd_drag_end


end module dd_drag_mod
