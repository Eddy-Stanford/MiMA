module dd_drag_mod

use forpy_mod

implicit none

public const_dd_drag_calc
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!           PUBLIC SUBROUTINES
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subroutine const_dd_drag_calc (is, js, lat, pfull, zfull, temp, uuu, vvv,  &
                         Time, delt, gwfcng_x, gwfcng_y)
    !--------------------------------------------------------------------  
    !    cg_drag_calc defines the arrays needed to
    !    calculate the convective
    !    gravity wave forcing, calls gwfc to calculate
    !    the forcing, returns 
    !    the desired output fields, and saves the
    !    values for later retrieval
    !    if they are not calculated on every timestep.
    !
    !---------------------------------------------------------------------

    !---------------------------------------------------------------------
    integer,                intent(in)      :: is, js
    real, dimension(:,:),   intent(in)      :: lat
    real, dimension(:,:,:), intent(in)      :: pfull, zfull, temp, uuu, vvv
    type(time_type),        intent(in)      :: Time
    real           ,        intent(in)      :: delt
    real, dimension(:,:,:), intent(out)     :: gwfcng_x, gwfcng_y

    !-------------------------------------------------------------------
    !    intent(in) variables:
    !
    !       is,js    starting subdomain i,j indices of
    !       data in 
    !                the physics_window being
    !                integrated
    !       lat      array of model latitudes at cell
    !       boundaries [radians]
    !       pfull    pressure at model full levels [ Pa  ]
    !       zfull    height at model full levels [ m  ]
    !       temp     temperature at model levels [ deg
    !       K ]
    !       uuu      zonal wind  [ m/s  ]
    !       vvv      meridional wind  [ m/s  ]
    !       Time     current time, needed for
    !       diagnostics [ time_type  ]
    !       delt     physics time step [ s  ]
    !
    !    intent(out) variables:
    !
    !       gwfcng_x time tendency for u eqn due to
    !       gravity-wave forcing
    !                [ m/s^2  ]
    !       gwfcng_y time tendency for v eqn due to
    !       gravity-wave forcing
    !                [ m/s^2  ]
    !
    !-------------------------------------------------------------------

    !-------------------------------------------------------------------
    !    local variables:
