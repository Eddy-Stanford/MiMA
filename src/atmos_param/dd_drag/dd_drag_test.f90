program test_wavenet
use dd_drag_mod, only: dd_drag_init, dd_drag_calc, dd_drag_end
implicit none

! Local Features
integer                           :: is, js, I
real,      dimension(64,128)     :: lat
real,      dimension(64,128,40)  :: pfull, zfull, temp, uuu, vvv
real                              :: Time, delta

! GWD
real,      dimension(64,128,40)  :: gwfcng_x, gwfcng_y

! Initialize Features
is = 3
js = 4
lat = 1
pfull = 1
zfull = 1
temp = 1
uuu = 1
vvv = 1

! Call dd_drag_init
call dd_drag_init()

! Call dd_drag_calc
call dd_drag_calc(is, js, lat, pfull, zfull, temp, uuu, vvv, Time, delta, gwfcng_x, gwfcng_y)
print *, "gwfu = ", gwfcng_x(2,2,2)
print *,"gwfv = ", gwfcng_y(2,2,2)


! Call dd_drag_end
call dd_drag_end()


end program
