module cg_drag_king_mod
  use time_manager_mod,  only:  time_manager_init, time_type
  use fms_mod, only: fms_init, mpp_pe, mpp_root_pe,  &
  file_exist, check_nml_error,  &
  error_mesg,  FATAL, WARNING, NOTE, &
  close_file, open_namelist_file, &
  stdlog, write_version_number, &
  read_data, write_data,   &
  open_restart_file, field_size
  use time_manager_mod,       only:  time_manager_init, time_type
  use diag_manager_mod,       only:  diag_manager_init,   &
  register_diag_field, send_data
  use constants_mod,          only:  constants_init, PI, RDGAS, GRAV, CP_AIR, &
  SECONDS_PER_DAY
  use transforms_mod,         only: grid_domain
  implicit none
  private 
  
  character(len=128)  :: version =  '$Id: cg_drag_king.F90,v 19.0 2024/08/16 $'
  character(len=128)  :: tagname =  '$Name: Robert C. King $'
  character(len=12) :: mod_name = 'cg_drag_king'
  logical          :: module_is_initialized=.false.
  !! NAME OF INPUT FILE
  character(len=32) :: data_name = 'cg_drag_king_params'
  
  !! NAME OF CW FIELD
  character(len=16) :: data_cw_field = 'cw'
  
  !! NAME OF BT FIELD
  character(len=16) :: data_bt_field = 'bt'
  
  !! NAME OF SOURCE LEVEL FIELD
  character(len=16) :: data_src_field = 'srclvl'
  
  !! missing value
  real             :: missing_value = -999.
  
  !! WAVE DRAG VARIABLES
  integer, allocatable, dimension(:,:) ::  source_level, damp_level
  real, allocatable, dimension(:,:,:) ::  gwd_u, gwd_v
  real, allocatable, dimension(:,:) :: source_amp !! TOTAL AMPLITUDE OF GRAVITY WAVE STRESS
  real, allocatable, dimension(:,:) :: half_width !! HALF WIDTH SPECTRUM
  
  
  !! GRAVITY WAVE PHYSICAL PARAMTERS
  real :: source_level_pressure=315.e+02 ! source level here at the equator.
  
  real :: damp_level_pressure=0.8e+02 ! any wave reaching the top level will be deposited down to this level
  
  real :: cmax = 99.6 ! maximum phase speed of gravity waves (m/s) 
  
  integer :: nk=1 ! number of vertical wavelengths in the spectrum
  
  real :: dc =1.2 ! Gravity wave spectral resolition (m/s)
  
  real :: Bw = 0.4 ! wave momentum flux amplitude (m/s)^2
  
  real :: kelvin_kludge=1.
  logical :: calculate_ked=.false. 
  
  
  !! other run time
  real, dimension(:), allocatable :: c0,kwv,k2
  integer :: nc !! number of wave speeds in spectrum (symmetric around c = 0)
  integer :: klevel_of_source, klevel_of_damp 
  ! k index of the gravity wave source level
  ! AND k index of the level where mesosphere drag is dumped
  
  
  integer     :: cg_drag_freq=0     ! calculation frequency [ s ]
  integer     :: cg_drag_offset=0   ! offset of calculation from 00Z [ s ]
  ! only has use if restarts are written
  ! at 00Z and calculations are not done
  ! every time step
  integer :: cgdrag_alarm ! time remaining unit next cg drag calculation 
  
  !---------------------------------------------------------------------
  !   variables for netcdf diagnostic fields.
  !---------------------------------------------------------------------
  integer          :: id_kedx_cgwd, id_kedy_cgwd, id_bf_cgwd, &
  id_gwfx_cgwd, id_gwfy_cgwd
  
  !! ---------------------------------------------------------------------
  !  Namelist for cg_drag_king_mod
  !! ---------------------------------------------------------------------
  
  namelist /cg_drag_king_nml/ & 
  source_level_pressure, &
  damp_level_pressure, cmax, & 
  nk, dc, Bw, &
  cg_drag_freq, cg_drag_offset
  
  !! ---------------------------------------------------------------------
  ! Subroutines for cg_drag_king_mod
  !! ---------------------------------------------------------------------
  public   cg_drag_king_init, cg_drag_king_calc, cg_drag_king_end, &
  cg_drag_king_time_vary, cg_drag_king_endts
  
  private gwfc
  
  contains
  
  subroutine cg_drag_king_init (lonb,latb,pref,Time,axes)
    !-------------------------------------------------------------------
    !   cg_drag_king_init is the constructor for cg_drag_king_mod.
    !-------------------------------------------------------------------
    
    !-------------------------------------------------------------------
    !   intent(in) variables:
    !
    !       lonb      1d array of model longitudes on cell corners [radians]
    !       latb      1d array of model latitudes at cell corners [radians]
    !       pref      array of reference pressures at full levels (plus
    !                 surface value at nlev+1), based on 1013.25hPa pstar
    !                 [ Pa ]
    !       Time      current time (time_type)
    !       axes      data axes for diagnostics
    !
    !------------------------------------------------------------------
    
    real,    dimension(:),   intent(in)      :: lonb, latb, pref
    integer, dimension(4),   intent(in)      :: axes
    type(time_type),         intent(in)      :: Time
    
    
    
    
    !-------------------------------------------------------------------
    !   local variables: 
    
    integer                 :: unit, ierr, io, logunit
    integer                 :: n, i, j, k
    integer                 :: idf, jdf, kmax
    real                    :: lat(size(lonb) - 1, size(latb) - 1)
    real                    :: thislatdeg
    integer, dimension(4) :: sizbt
    integer, dimension(4) :: sizcw 
    !-------------------------------------------------------------------    
    !  constants
    real                    :: pif = 3.14159265358979/180.
    real                    :: pifinv = 180./3.14159265358979
    
    !-------------------------------------------------------------------
    !-------------------------------------------------------------------
    !   local variables: 
    !   
    !       unit           unit number for nml file 
    !       ierr           error return flag 
    !       io             error return code 
    !       n              loop index
    !       k              loop index
    !       idf            number of i points on this processor
    !       jdf            number of j points on this processor
    !       kmax           number of k points on this processor
    !
    !---------------------------------------------------------------------
    if (module_is_initialized) return ! return if already initialized
    
    ! initalize all preqrequesites
    call fms_init
    call time_manager_init
    call diag_manager_init
    call constants_init
    
    ! read namelist
    if (file_exist('input.nml')) then
      unit = open_namelist_file()
      ierr=1
      do while(ierr /= 0)
        read(unit,nml=cg_drag_king_nml,iostat=io,end=10)
        ierr = check_nml_error(io,'cg_drag_king_nml')
      enddo
      10          call close_file(unit)
    endif
    
    call write_version_number (version, tagname)
    logunit = stdlog()
    if (mpp_pe() == mpp_root_pe()) then
      write (logunit, nml=cg_drag_king_nml)
    endif 
    kmax = size(pref) -1
    jdf = size(latb) - 1
    idf = size(lonb) - 1
    
    allocate( source_level(idf,jdf))
    allocate( damp_level(idf,jdf))
    allocate(source_amp(idf,jdf))
    allocate(half_width(idf,jdf))
    
    if(file_exist('INPUT/cg_drag_king_params.nc')) then
      if(mpp_pe() == mpp_root_pe()) then
        call error_mesg('cg_drag_king_mod','Reading NETCDF parameter file: INPUT/cg_drag_king_params.nc',NOTE)
      endif
      call field_size('INPUT/cg_drag_king_params.nc', data_bt_field, sizbt)
      call field_size('INPUT/cg_drag_king_params.nc', data_cw_field, sizcw)
      
      call read_data('INPUT/cg_drag_king_params.nc', data_bt_field, source_amp,domain=grid_domain)
      call read_data('INPUT/cg_drag_king_params.nc', data_cw_field, half_width,domain=grid_domain)
    else 
      call error_mesg('cg_drag_king_mod','No NETCDF parameter file found',FATAL)
    endif
    
    
    !! descend from top layer to bottom
    do k=1,kmax
      
      if (pref(k) < damp_level_pressure) then
        klevel_of_damp = k
      endif 
      !! Get the model level of the source level, it is the level that is just VERTICALY BELOW the source level
      if (pref(k) > source_level_pressure) then
        ! above source level
        klevel_of_source = k 
        exit 
      endif 
    enddo 
    
    do j=1,jdf
      do i=1,idf
        lat(i,j) = 0.5*(latb(j) + latb(j+1))
        !! Not sure about this?
        source_level(i,j) = (kmax + 1) - ( (kmax+1 - klevel_of_source)*cos(lat(i,j)) + 0.5)
        
        ! fixed dampign level
        damp_level(i,j) = klevel_of_damp
      enddo
    enddo 
    source_level = min(source_level,kmax-1)
    damp_level = min(damp_level,kmax)
    
    !! ALLOCATE phase speed bins
    nc = 2.0*cmax/dc + 1
    allocate ( c0(nc) )
    do n=1,nc
      c0(n) = (n-1)*dc - cmax
    end do
    
    !! ALLOCATE wavenumber
    allocate (kwv(nk))
    allocate (k2(nk))
    do n=1,nk
      kwv(n) = 2.*PI/(30.*(10.**n)*1.e3) ! largest wavenumber is 30 km
      k2(n) = kwv(n)**2
    enddo 
    
    id_bf_cgwd =  &
    register_diag_field (mod_name, 'bf_cgwd', axes(1:3), Time, &
    'buoyancy frequency from cg_drag_king', ' /s',   &
    missing_value=missing_value)
    id_gwfx_cgwd =  &
    register_diag_field (mod_name, 'gwfu_cgwd', axes(1:3), Time, &
    'gravity wave forcing on mean zonal flow', &
    'm/s^2',  missing_value=missing_value)
    id_gwfy_cgwd =  &
    register_diag_field (mod_name, 'gwfv_cgwd', axes(1:3), Time, &
    'gravity wave forcing on mean meridional flow', &
    'm/s^2',  missing_value=missing_value)
    id_kedx_cgwd =  &
    register_diag_field (mod_name, 'kedx_cgwd', axes(1:3), Time, &
    'effective eddy viscosity from cg_drag_king', 'm^2/s',   &
    missing_value=missing_value)
    id_kedy_cgwd =  &
    register_diag_field (mod_name, 'kedy_cgwd', axes(1:3), Time, &
    'effective eddy viscosity from cg_drag_king', 'm^2/s',   &
    missing_value=missing_value)
    
    allocate ( gwd_u(idf,jdf,kmax) )
    allocate ( gwd_v(idf,jdf,kmax) )
    
    if( cg_drag_freq /= 0 ) then
      if( modulo(86400,cg_drag_freq) /= 0 ) then
        call error_mesg('cg_drag_king','cg_drag_freq must divide 86400 (full day) or equal 0', FATAL)
      endif
    endif
    
    ! set GWD for u/v to 0.0
    
    gwd_u(:,:,:) = 0.0
    gwd_v(:,:,:) = 0.0
    
    
    if (cg_drag_offset > 0) then
      cgdrag_alarm = cg_drag_offset
    else 
      cgdrag_alarm = cg_drag_freq
    endif
    
    module_is_initialized = .true.
  end subroutine cg_drag_king_init
  
  subroutine cg_drag_king_time_vary (delt)
    
    real           ,        intent(in)      :: delt
    
    !---------------------------------------------------------------------
    !    decrement the time remaining until the next cg_drag_king calculation.
    !---------------------------------------------------------------------
    cgdrag_alarm = cgdrag_alarm - delt
    
    !---------------------------------------------------------------------
    
  end subroutine cg_drag_king_time_vary
  
  subroutine cg_drag_king_endts
    
    !--------------------------------------------------------------------
    !    if this was a calculation step, reset cgdrag_alarm to indicate 
    !    the time remaining before the next calculation of gravity wave 
    !    forcing.
    !--------------------------------------------------------------------
    if (cgdrag_alarm <= 0 ) then
      cgdrag_alarm = cgdrag_alarm + cg_drag_freq
    endif
    
  end subroutine cg_drag_king_endts
  
  subroutine cg_drag_king_calc (is, js, lat, pfull, zfull, temp, uuu, vvv,  &
    Time, delt, gwfcng_x, gwfcng_y)
    !--------------------------------------------------------------------  
    !    cg_drag_king_calc defines the arrays needed to calculate the convective
    !    gravity wave forcing, calls gwfc to calculate the forcing, returns 
    !    the desired output fields, and saves the values for later retrieval
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
    !       is,js    starting subdomain i,j indices of data in 
    !                the physics_window being integrated
    !       lat      array of model latitudes at cell boundaries [radians]
    !       pfull    pressure at model full levels [ Pa ]
    !       zfull    height at model full levels [ m ]
    !       temp     temperature at model levels [ deg K ]
    !       uuu      zonal wind  [ m/s ]
    !       vvv      meridional wind  [ m/s ]
    !       Time     current time, needed for diagnostics [ time_type ]
    !       delt     physics time step [ s ]
    !
    !    intent(out) variables:
    !
    !       gwfcng_x time tendency for u eqn due to gravity-wave forcing
    !                [ m/s^2 ]
    !       gwfcng_y time tendency for v eqn due to gravity-wave forcing
    !                [ m/s^2 ]
    !
    !-------------------------------------------------------------------
    
    !-------------------------------------------------------------------
    !    local variables:
    
    
    real,    dimension (size(uuu,1), size(uuu,2), size(uuu,3))  ::  dtdz, ked_gwfc_x, ked_gwfc_y
    
    real,    dimension (size(uuu,1),size(uuu,2), 0:size(uuu,3)) ::  &
    zzchm, zu, zv, zden, zbf,    &
    gwd_xtnd, ked_xtnd, &
    gwd_ytnd, ked_ytnd
    
    integer           :: iz0
    logical           :: used
    real              :: bflim = 2.5E-5
    integer           :: ie, je
    integer           :: imax, jmax, kmax
    integer           :: i, j, k, nn
    real              :: pif = 3.14159265358979/180.
    
    !---------------------------------------------------------------------
    !    define processor extents and loop limits.
    !---------------------------------------------------------------------
    imax = size(uuu,1)
    jmax = size(uuu,2)
    kmax = size(uuu,3)
    
    ie = is + imax - 1
    je = js + jmax - 1
    
    
    if (cgdrag_alarm <= 0) then
      do j=1,jmax 
        do i=1,imax
          iz0 = source_level(i+is-1,j+js-1)
          
          !! lapse rate
          ! forward difference
          dtdz(i,j,1) = (temp(i,j,1) - temp(i,j,2))/(zfull(i,j,1) - zfull(i,j,2))
          do k=2,iz0
            ! central differnce
            dtdz(i,j,k) = (temp(i,j,k-1) - temp(i,j,k+1))/(zfull(i,j,k-1) - zfull(i,j,k+1))
          enddo
          
          !! air density
          do k=1,iz0+1
            zden(i,j,k) = pfull(i,j,k)/(temp(i,j,k)*RDGAS)
          enddo 
          
          !! bouyancy frequencu
          do k=1,iz0
            zbf(i,j,k) = (GRAV/temp(i,j,k))*(dtdz(i,j,k) + GRAV/CP_AIR)
            
            !! Buyoyancy frequency is not allowed to be less than bflim 
            if (zbf(i,j,k) < bflim) then
              zbf(i,j,k) = sqrt(bflim)
            else 
              zbf(i,j,k) = sqrt(zbf(i,j,k))
            endif
          enddo
          
          ! initialze zbf for remainig levels 
          if (id_bf_cgwd > 0) then
            zbf(i,j,iz0+1:) = 0.0
          endif
          
          !! Define local variables for the gravity wave forcing
          do k=1,iz0+1
            zzchm(i,j,k) = zfull(i,j,k)
          end do
          do k=1,iz0   
            zu(i,j,k) = uuu(i,j,k)
            zv(i,j,k) = vvv(i,j,k)
          end do
          
          ! So we have a SNEAKY level above the model top so we can calcualte 
          ! the gravity wave forcing at the top level
          ! use second order extrapolation
          zzchm(i,j,0) = 2.*zzchm(i,j,1)  - zzchm(i,j,2)!
          zu(i,j,0) = 2.*zu(i,j,1) - zu(i,j,2)
          zv(i,j,0) = 2.*zv(i,j,1) - zv(i,j,2)
          ! geometric mean
          zden(i,j,0) = zden(i,j,1)*zden(i,j,1)/zden(i,j,2)
          !constant f at model top
          zbf(i,j,0) = zbf(i,j,1)
          
          
        enddo
      enddo
      !---------------------------------------------------------------------
      !    pass the vertically-extended input arrays to gwfc. gwfc will cal-
      !    culate the gravity-wave forcing and, if desired, an effective eddy 
      !    diffusion coefficient at each level above the source level. output
      !    is returned in the vertically-extended arrays gwfcng and ked_gwfc.
      !    upon return move the output fields into model-sized arrays. 
      !---------------------------------------------------------------------
      
      call gwfc (is, ie, js, je, damp_level, source_level, source_amp, lat,   &
      zden, zu, zbf,zzchm, gwd_xtnd, ked_xtnd)
      
      gwfcng_x  (:,:,1:kmax) = gwd_xtnd(:,:,1:kmax  )
      ked_gwfc_x(:,:,1:kmax) = ked_xtnd(:,:,1:kmax  )
      
      call gwfc (is, ie, js, je, damp_level, source_level, source_amp,  lat,  &
      zden, zv, zbf,zzchm, gwd_ytnd, ked_ytnd)
      gwfcng_y  (:,:,1:kmax) = gwd_ytnd(:,:,1:kmax  )
      ked_gwfc_y(:,:,1:kmax) = ked_ytnd(:,:,1:kmax  )
      
      !--------------------------------------------------------------------
      !    store the gravity wave forcing into a processor-global array. FOR THE MODULE 
      !-------------------------------------------------------------------
      gwd_u(is:ie,js:je,:) = gwfcng_x(:,:,:)
      gwd_v(is:ie,js:je,:) = gwfcng_y(:,:,:)
      !--------------------------------------------------------------------
      !    if activated, store the effective eddy diffusivity into a 
      !    processor-global array, and if desired as a netcdf diagnostic, 
      !    send the data to diag_manager_mod.
      !--------------------------------------------------------------------
      
      if (id_kedx_cgwd > 0) then
        used = send_data (id_kedx_cgwd, ked_gwfc_x, Time, is, js, 1)
      endif
      
      if (id_kedy_cgwd > 0) then
        used = send_data (id_kedy_cgwd, ked_gwfc_y, Time, is, js, 1)
      endif
      !--------------------------------------------------------------------
      !    save any other netcdf file diagnostics that are desired.
      !--------------------------------------------------------------------
      if (id_bf_cgwd > 0) then
        used = send_data (id_bf_cgwd,  zbf(:,:,1:), Time, is, js )
      endif
      
      if (id_gwfx_cgwd > 0) then
        used = send_data (id_gwfx_cgwd, gwfcng_x, Time, is, js, 1)
      endif
      if (id_gwfy_cgwd > 0) then
        used = send_data (id_gwfy_cgwd, gwfcng_y, Time, is, js, 1)
      endif
      
    else  ! cdrag alarm is NOT READY YET 
      gwfcng_x(:,:,:) = gwd_u(is:ie,js:je,:)
      gwfcng_y(:,:,:) = gwd_v(is:ie,js:je,:)
    endif 
    
    call cg_drag_king_endts
    call cg_drag_king_time_vary(delt)
    
  end subroutine cg_drag_king_calc
  
  subroutine cg_drag_king_end   
    module_is_initialized = .false.
  end subroutine cg_drag_king_end
  
  subroutine gwfc (is, ie, js, je, damp_level, source_level, source_amp, lat, rho, u, &
    bf, z, gwf, ked)
    
    !-------------------------------------------------------------------
    !    subroutine gwfc computes the gravity wave-driven-forcing on the
    !    zonal wind given vertical profiles of wind, density, and buoyancy 
    !    frequency. 
    !    Based on version implemented in SKYHI -- 27 Oct 1998 by M.J. 
    !    Alexander and L. Bruhwiler.
    !-------------------------------------------------------------------
    !-------------------------------------------------------------------
    integer,                     intent(in)             :: is, ie, js, je
    integer, dimension(:,:),     intent(in)             :: source_level, damp_level
    real,    dimension(:,:),     intent(in)             :: source_amp, lat
    real,    dimension(:,:,0:),  intent(in)             :: rho, u, bf, z
    real,    dimension(:,:,0:),  intent(out)            :: gwf
    real,    dimension(:,:,0:),  intent(out)            :: ked
    
    !-------------------------------------------------------------------
    !  intent(in) variables:
    !
    !      is, ie, js, je   starting/ending subdomain i,j indices of data
    !                       in the physics_window being integrated
    !      source_level     k index of model level serving as gravity wave
    !                       source
    !      damp_level       k index of the lowest model level at which all drag that reaches the model top is partially dumped
    !      source_amp     amplitude of  gravity wave source [Pa]
    !                       
    !      rho              atmospheric density [ kg/m^3 ] 
    !      u                zonal wind component [ m/s ]
    !      bf               buoyancy frequency [ /s ]
    !      z                height of model levels  [ m ]
    !
    !  intent(out) variables:
    !
    !      gwf              gravity wave forcing in u equation  [ m/s^2 ]
    !
    !  intent(out), optional variables:
    !
    !      ked              eddy diffusion coefficient from gravity wave 
    !                       forcing [ m^2/s ]
    !
    !------------------------------------------------------------------
    
    !  local variables
    real,    dimension (0:size(u,3)-1 ) :: wv_frcng, diff_coeff, c0mu 
    integer, dimension (nc) ::   msk
    real   , dimension (nc) ::   c0mu0, B0
    real                    ::   fm, fe, Hb, alp2, Foc, c, test, rbh,&
    eps, Bsum, fp0, omc,dz,cw
    integer                 ::   iz0, iztop
    integer                 ::   i, j, k, ink, n
    real                    ::   ampl
    real                    :: pifinv = 180./3.14159265358979
    !------------------------------------------------------------------
    !  local variables:
    ! 
    !      wv_frcng    gravity wave forcing tendency [ m/s^2 ]
    !      diff_coeff  eddy diffusion coefficient [ m2/s ]
    !      c0mu        difference between phase speed of wave n and u 
    !                  [ m/s ]
    !      dz          delta z between model levels [ m ]
    !      fp0        factor used in determining if wave is breaking 
    !                  [ s/m ]
    !      omc         critical frequency that marks total internal 
    !                  reflection  [ /s ]
    !      msk         indicator as to whether wave n is still propagating 
    !                  upwards (msk=1), or has been removed from the 
    !                  spectrum because of breaking or reflection (msk=0)
    !      c0mu0       difference between phase speed of wave n and u at the
    !                  source level [ m/s ]
    !      B0          wave momentum flux amplitude for wave n [ (m/s)^2 ]
    !      fm          used to sum up momentum flux from all waves n 
    !                  deposited at a level [ (m/s)^2 ]
    !      fe          used to sum up contributions to diffusion coefficient
    !                  from all waves n at a level [ (m/s)^3 ]
    !      Hb          density scale height [ m ]
    !      alp2        scale height factor: 1/(2*Hb)**2  [ /m^2 ]
    !      Foc         wave breaking threshold [ s/m ]
    !      c           wave phase speed used in defining wave momentum flux
    !                  amplitude [ m/s ]
    !      test        condition defining internal reflection [ /s ]
    !      rbh         atmospheric density at half-level (geometric mean)
    !                  [ kg/m^3 ]
    !      eps         intermittency factor
    !      Bsum        total mag of gravity wave momentum flux at source 
    !                  level, divided by the density  [ m^2/s^2 ]
    !      iz0         source level vertical index for the given column
    !      i,j,k       spatial do loop indices
    !      ink         wavenumber loop index 
    !      n           phase speed loop index 
    !      ampl        gravity wave stress [Pa] 
    !
    !--------------------------------------------------------------------
    
    
    !-------------------------------------------------------------------
    !    initialize the output arrays. these will hold values at each 
    !    (i,j,k) point, summed over the wavelengths and phase speeds
    !    defining the gravity wave spectrum.
    !-------------------------------------------------------------------
    gwf = 0.0
    ked = 0.0
    
    do j=1,size(u,2)
      do i=1,size(u,1)
        
        iz0 = source_level(i+is-1,j+js-1)
        iztop = damp_level(i+is-1,j+js-1)
        
        
        ampl = source_amp(i+is-1,j+js-1)
        cw = half_width(i+is-1,j+js-1)
        Bsum = 0.
        do n=1,nc
          ! intrinistic velocity at source level
          c0mu0(n) = c0(n) - u(i,j,iz0) 
          if(c0mu0(n) ==0.0) then
            B0(n) = 0.0
          else 
            B0(n) = sign(1.0,c0mu0(n))*(Bw*exp(-log(2.0)*(c0mu0(n)/cw)**2))
          endif
          Bsum = Bsum + abs(B0(n))
        enddo  
        
        if (Bsum == 0.0) then
          call error_mesg('cg_drag_king_mod', & 
          ' zero flux input at source level',FATAL)
        endif 
        eps = (ampl)/(Bsum*nk*rho(i,j,iz0))
        do ink=1,nk ! loop over wavelengths
          ! flag that wave is still propagating upwards
          msk =1 
          do k=iz0,0,-1 ! backwards loop from source to top
            fp0 = 0.5*(rho(i,j,k)/rho(i,j,iz0))*kwv(ink)/bf(i,j,k)
            
            dz = z(i,j,k) - z(i,j,k+1)
            Hb = -dz/log(rho(i,j,k)/rho(i,j,k+1)) ! scale height
            alp2 = 0.25/(Hb*Hb) ! scale height factor
            
            ! frequency for total internal reflection if omega > omc 
            omc = sqrt(bf(i,j,k)*bf(i,j,k)*k2(ink)/(alp2 + k2(ink)))
            fm = 0.0
            fe = 0.0
            do n=1,nc ! phase speed loop
              if(msk(n) == 1) then ! at THIS level, the wave at phase speed c is still propegating
                c0mu(k) = c0(n) - u(i,j,k) ! c - u(z) not at source!!
                
                if (c0mu(k) == 0.) then
                  msk(n) = 0 ! breaking has occured at this level
                  cycle
                endif 
                ! TIR test using the critical omega from before 
                test = abs(c0mu(k))*kwv(ink) - omc
                if (test >= 0.0) then ! TIR has occured so wave not at this level anymore
                  msk(n) = 0
                  cycle
                endif 
                
                if(k==0) then! model top
                  msk(n) =0
                  if( k < iz0) then ! is this not always true?
                    fm = fm+ B0(n)
                    fe = fe + c0mu(k)*B0(n)
                  endif
                  cycle 
                endif 
                
                ! has the wave broken at this level?
                Foc = B0(n)/(c0mu(k))**3 - fp0
                if ((Foc >= 0.0) .or. (c0mu0(n)*c0mu(k) <= 0.0)) then
                  msk(n) = 0
                  if (k < iz0) then
                    fm = fm + B0(n)
                    fe = fe + c0mu(k)*B0(n)
                  endif
                endif 
                
              endif
            enddo
            ! okay we have checked all the wavenumbers, now we can calculate the gravity wave forcing at this level
            if (k < iz0) then
              rbh = sqrt(rho(i,j,k) * rho(i,j,k+1))
              wv_frcng(k) = (rho(i,j,iz0)/rbh)*fm*eps/dz
              diff_coeff(k) = (rho(i,j,iz0)/rbh)*fe*eps/(dz*bf(i,j,k)*bf(i,j,k))
              
              ! average over both layers where waves are breaking. 
              if (k==0) then
                !model top
                wv_frcng(k+1) = 0.5*wv_frcng(k+1) ! this seems a little contentious?
                diff_coeff(k+1) =  0.5*diff_coeff(k+1)
              else 
                wv_frcng(k+1) = 0.5*(wv_frcng(k+1) + wv_frcng(k))
                diff_coeff(k+1) = 0.5*(diff_coeff(k+1) + diff_coeff(k))
              endif 
            else 
              wv_frcng(iz0) = 0.0
              diff_coeff(iz0) = 0.0
            endif 
            
          enddo ! height loop 
          
          do k=1,iztop 
            !distribtue the forcing amongst all layers above iztop
            wv_frcng(k) = wv_frcng(k) + wv_frcng(0)/REAL(iztop)
            diff_coeff(k) = diff_coeff(k) + diff_coeff(0)/REAL(iztop)
          enddo 
          do k=0,iz0
            gwf(i,j,k) = gwf(i,j,k) + wv_frcng(k)
            ked(i,j,k) = ked(i,j,k) + diff_coeff(k)
          enddo 
        enddo
      enddo
    enddo 
  end subroutine gwfc
end module cg_drag_king_mod