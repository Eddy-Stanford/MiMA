
!> Utility routines for atmospheric tracers: wet and dry deposition, and interpolation of
!> emission fields.
!>
!> The deposition schemes provide consistent removal mechanisms for the tracers; they are
!> selected per tracer with `dry_deposition` and `wet_deposition` methods in the
!> `field_table`. The deposition fluxes are available as diagnostics of the module
!> `tracers`: the tracer name followed by `ddep` (dry deposition), and by `wdep_ls` and
!> `wdep_cv` (wet deposition by large-scale condensation and by convection).
!>
!> Original authors: William Cooke.
module atmos_tracer_utilities_mod

  use fms_mod, only: lowercase, &
                     write_version_number, &
                     stdlog, &
                     mpp_pe, &
                     mpp_root_pe, &
                     error_mesg, &
                     NOTE
  use time_manager_mod, only: time_type
  use diag_manager_mod, only: send_data, &
                              register_diag_field
  use tracer_manager_mod, only: query_method, &
                                get_tracer_names, &
                                get_number_tracers, &
                                MAX_TRACER_FIELDS
  use field_manager_mod, only: MODEL_ATMOS, parse
  use constants_mod, only: grav, rdgas, PI
  use horiz_interp_mod, only: horiz_interp
  use constants_mod, only: PI
  use mima_interpolator_mod, only: interpolator, &
                                   interpolate_type

  implicit none
  private
!-----------------------------------------------------------------------
!----- interfaces -------

  public wet_deposition, &
    dry_deposition, &
    interp_emiss, &
    atmos_tracer_utilities_end, &
    atmos_tracer_utilities_init

!---- version number -----
  logical :: module_is_initialized = .false.

  character(len=128) :: version = '$Id: atmos_tracer_utilities.f90,v 12.0 2005/04/14 15:52:18 fms Exp $'
  character(len=128) :: tagname = '$Name: lima $'

  character(len=7), parameter :: mod_name = 'tracers'
!-----------------------------------------------------------------------
!--- identification numbers for  diagnostic fields and axes ----
  integer, parameter :: max_tracers = MAX_TRACER_FIELDS
  integer :: id_tracer_ddep(max_tracers), id_tracer_wdep_ls(max_tracers), id_tracer_wdep_cv(max_tracers)
  character(len=32), dimension(max_tracers) :: tracer_names = ' '
  character(len=32), dimension(max_tracers) :: tracer_units = ' '
  character(len=128), dimension(max_tracers) :: tracer_longnames = ' '
  character(len=32), dimension(max_tracers) :: tracer_wdep_names = ' '
  character(len=32), dimension(max_tracers) :: tracer_wdep_units = ' '
  character(len=128), dimension(max_tracers) :: tracer_wdep_longnames = ' '
  character(len=32), dimension(max_tracers) :: tracer_ddep_names = ' '
  character(len=32), dimension(max_tracers) :: tracer_ddep_units = ' '
  character(len=128), dimension(max_tracers) :: tracer_ddep_longnames = ' '

  real, allocatable :: blon_out(:), blat_out(:)
!----------------parameter values for the diagnostic units--------------
  real, parameter :: mw_air = 0.0289644
  real, parameter :: Navo = 6.023e23

contains

!
! ######################################################################
!
  !> Initializes the module: registers the dry and wet deposition diagnostics of the tracers.
  !>
  !> The diagnostic names are the tracer name followed by `ddep` for the dry deposition and
  !> `wdep_ls`, `wdep_cv` for the wet deposition; they can be entered in the `diag_table`
  !> under the module name `tracers`. The units of the deposition fields are kg/m2/s for
  !> tracers in `mmr` or `kg/kg`, and mole/m2/s for tracers in `vmr`, `mol/mol` or `mole/mole`.
  subroutine atmos_tracer_utilities_init(lonb, latb, mass_axes, Time)

    real, dimension(:), intent(in) :: lonb, latb  !! longitudes and latitudes of the cell corners [rad]
    integer, dimension(3), intent(in) :: mass_axes  !! diagnostic axes (lon, lat, pfull)
    type(time_type), intent(in) :: Time  !! model time

    integer :: ntrace
    character(len=20) :: units = ''
!
    integer :: n, unit
    character(len=128) :: name

! Make local copies of the local domain dimensions for use
! in interp_emiss.
    allocate (blon_out(size(lonb(:))))
    allocate (blat_out(size(latb(:))))
!      allocate ( data_out(size(lonb(:))-1, size(latb(:))-1))
    blon_out = lonb
    blat_out = latb

    do n = 1, max_tracers
      write (tracer_names(n), 100) n
      write (tracer_longnames(n), 102) n
      tracer_units(n) = 'none'
    end do
100 format('tr', i2.2)
102 format('tracer ', i2.2)

    call get_number_tracers(MODEL_ATMOS, num_tracers=ntrace)
    do n = 1, ntrace
!--- set tracer tendency names where tracer names have changed ---

      call get_tracer_names(MODEL_ATMOS, n, tracer_names(n), tracer_longnames(n), tracer_units(n))
      write (name, 100) n
      if (trim(tracer_names(n)) /= name) then
        tracer_ddep_names(n) = trim(tracer_names(n))//'ddep'
        tracer_wdep_names(n) = trim(tracer_names(n))//'wdep'
      end if
      write (name, 102) n
      if (trim(tracer_longnames(n)) /= name) then
        tracer_wdep_longnames(n) = &
          trim(tracer_longnames(n))//' wet deposition for tracers'
        tracer_ddep_longnames(n) = &
          trim(tracer_longnames(n))//' dry deposition for tracers'
      end if

      select case (trim(tracer_units(n)))
      case ('mmr')
        units = 'kg/m2/s'
      case ('kg/kg')
        units = 'kg/m2/s'
      case ('vmr')
        units = 'mole/m2/s'
      case ('mol/mol')
        units = 'mole/m2/s'
      case ('mole/mole')
        units = 'mole/m2/s'
      case default
        units = trim(tracer_units(n))//' kg/m2/s'
        call error_mesg('atmos_tracer_utilities_init', &
                        ' Dry dep units set to '//trim(units)//' in atmos_tracer_utilities for '//trim(tracer_names(n)), &
                        NOTE)
      end select

      id_tracer_ddep(n) = register_diag_field(mod_name, &
                                              trim(tracer_ddep_names(n)), mass_axes(1:2), Time, &
                                              trim(tracer_ddep_longnames(n)), &
                                              trim(units), missing_value=-999.)
      id_tracer_wdep_ls(n) = register_diag_field(mod_name, &
                                                 trim(tracer_wdep_names(n))//'_ls', mass_axes(1:2), Time, &
                                                 trim(tracer_wdep_longnames(n))//' in large scale', &
                                                 trim(units), missing_value=-999.)
      id_tracer_wdep_cv(n) = register_diag_field(mod_name, &
                                                 trim(tracer_wdep_names(n))//'_cv', mass_axes(1:2), Time, &
                                                 trim(tracer_wdep_longnames(n))//' in convective scheme', &
                                                 trim(units), missing_value=-999.)
    end do

    call write_version_number(version, tagname)

    if (mpp_pe() == mpp_root_pe()) then
      call write_namelist_values(stdlog(), ntrace)
    end if

    module_is_initialized = .true.

  end subroutine atmos_tracer_utilities_init
!
!#######################################################################
!
  !> Writes the names, long names and units of the deposition diagnostics to the log file.
  subroutine write_namelist_values(unit, ntrace)
    integer, intent(in) :: unit, ntrace
    integer :: n

    write (unit, 10)
    do n = 1, ntrace
      write (unit, 11) trim(tracer_wdep_names(n)), &
        trim(tracer_wdep_longnames(n)), &
        trim(tracer_wdep_units(n))
      write (unit, 11) trim(tracer_ddep_names(n)), &
        trim(tracer_ddep_longnames(n)), &
        trim(tracer_ddep_units(n))
    end do

10  format(' &TRACER_DIAGNOSTICS_NML', &
           /, '    TRACER:  names  longnames  (units)')
11  format(a16, 2x, a, 2x, '(', a, ')')

  end subroutine write_namelist_values

!
!#######################################################################
!
  !> Computes the tendency of a tracer in the lowest model level due to dry deposition, and
  !> sends the dry deposition flux diagnostic.
  !>
  !> Two types of dry deposition are coded:
  !>
  !> 1. Wind-driven dry deposition velocity. The deposition is modelled as a resistance
  !>    problem: the total resistance is `R = Ra + Rb`, with the aerodynamic resistance
  !>    `Ra = |u|/u_star**2` and the surface resistance `Rb = surfr/u_star` (laminar layer plus
  !>    uptake; `u_star` is at least 0.1 m/s), and the deposition velocity is `Vd = 1/R`.
  !> 2. Fixed dry deposition velocity. The deposition velocity does not change, but the
  !>    variation of the depth of the surface layer implies that there is variation in the
  !>    amount deposited.
  !>
  !> To use it, add one of the following as a method for the tracer in the field table:
  !>
  !> * `"dry_deposition","wind_driven","surfr=XXX"`, where XXX is the surface resistance
  !>   coefficient (default 500);
  !> * `"dry_deposition","fixed","land=XXX, sea=YYY"`, where XXX and YYY are the dry
  !>   deposition velocities [m/s] over land and over sea.
  subroutine dry_deposition(n, is, js, u, v, T, pwt, pfull, &
                            u_star, landmask, dsinku, tracer, Time)!, dry)
    integer, intent(in)                 :: n, is, js
    !! `n`: tracer number; `is`, `js`: start indices of the arrays in the processor domain
    real, intent(in), dimension(:, :)    :: u, v, T, pwt, pfull, u_star, tracer
    !! in the lowest model level: `u`, `v`: zonal and meridional wind [m/s]; `T`: temperature
    !! [K]; `pwt`: pressure weight dp/grav [kg/m2]; `pfull`: pressure [Pa]; `u_star`: friction
    !! velocity [m/s]; `tracer`: tracer mixing ratio
    logical, intent(in), dimension(:, :) :: landmask  !! true over land
    type(time_type), intent(in)         :: Time  !! model time, for the diagnostic
!type(interpolate_type), intent(inout) :: dry
    real, intent(out), dimension(:, :)   :: dsinku
    !! amount of tracer in the lowest model level which is dry deposited per second (tracer
    !! mixing ratio per second)

    real, dimension(size(u, 1), size(u, 2)) :: hwindv, frictv, resisa, xxfm, dz, dry_data
    integer :: i, j, flagsr
    real    :: land_dry_dep_vel, sea_dry_dep_vel, surfr
    logical :: used, flag
    integer :: flag_species
    character(len=10) ::units, names
    character(len=80) :: name, control, scheme, speciesname

! Default zero
    dsinku = 0.0
    flag = query_method('dry_deposition', MODEL_ATMOS, n, name, control)

    if (.not. flag) return

! delta z = dp/(rho * grav)
! delta z = RT/g*dp/p    pwt = dp/g
    dz(:, :) = pwt(:, :)*rdgas*T(:, :)/pfull(:, :)

    call get_drydep_param(name, control, scheme, land_dry_dep_vel, sea_dry_dep_vel)

    select case (lowercase(scheme))

    case ('wind_driven')
! Calculate horizontal wind velocity and aerodynamic resistance:
!   where xxfm=(u*/u) is drag coefficient, Ra=u/(u*^2),
!   and  u*=sqrt(momentum flux)  is friction velocity.
!
!****  Compute dry sinks (loss frequency, need modification when
!****    different vdep values are to be used for species)
      flagsr = parse(control, 'surfr', surfr)
      if (flagsr == 0) surfr = 500.
      hwindv = sqrt(u**2 + v**2)
      frictv = u_star
      resisa = hwindv/(u_star*u_star)
      where (frictv .lt. 0.1) frictv = 0.1
      dsinku = (1./(surfr/frictv + resisa))/dz

    case ('fixed')
! For the moment let's try to calculate the delta-z of the bottom
! layer and using a simple dry deposition velocity times the
! timestep, idt, calculate the fraction of the lowest layer which
! deposits.
      where (landmask(:, :))
! dry dep value over the land surface divided by the height of the box.
        dsinku(:, :) = land_dry_dep_vel/dz(:, :)
      elsewhere
! dry dep value over the sea surface divided by the height of the box.
        dsinku(:, :) = sea_dry_dep_vel/dz(:, :)
      end where

    case ('file')
      flag_species = parse(control, 'name', speciesname)
      if (flag_species > 0) then
        name = trim(speciesname)
      else
        call get_tracer_names(MODEL_ATMOS, n, name)
      end if
!chemistry start
!        call interpolator(dry,Time,dry_data, trim(name),is,js)
      dsinku(:, :) = dry_data(:, :)/dz(:, :)
!chemistry end
    case ('default')
    end select

    dsinku(:, :) = max(dsinku(:, :), 0.0e+00)
    where (tracer > 0)
      dsinku = dsinku*tracer
    elsewhere
      dsinku = 0.0
    end where

! Now save the dry deposition to the diagnostic manager
! delta z = dp/(rho * grav)
! delta z *rho  = dp/g
! tracer(kgtracer/kgair) * dz(m)* rho(kgair/m3) = kgtracer/m2
! so rho drops out of the equation
    if (id_tracer_ddep(n) > 0) then
      call get_tracer_names(MODEL_ATMOS, n, names, units=units)
      select case (trim(units))
      case ('mmr')
        used = send_data(id_tracer_ddep(n), dsinku*pwt, Time, &
                         is_in=is, js_in=js)
      case ('kg/kg')
        used = send_data(id_tracer_ddep(n), dsinku*pwt, Time, &
                         is_in=is, js_in=js)
      case ('vmr')
        used = send_data(id_tracer_ddep(n), dsinku*pwt/mw_air, Time, &
                         is_in=is, js_in=js)
      case ('mol/mol')
        used = send_data(id_tracer_ddep(n), dsinku*pwt/mw_air, Time, &
                         is_in=is, js_in=js)
      case ('mole/mole')
        used = send_data(id_tracer_ddep(n), dsinku*pwt/mw_air, Time, &
                         is_in=is, js_in=js)
      case default
        used = send_data(id_tracer_ddep(n), dsinku*pwt, Time, &
                         is_in=is, js_in=js)
      end select

    end if
  end subroutine dry_deposition
!
!#######################################################################
!
  !> Computes the tendency of a tracer due to wet deposition, and sends the wet deposition
  !> flux diagnostic.
  !>
  !> Schemes allowed here are:
  !>
  !> 1. `fraction`: the tracer is removed in the same fractional amount as the modelled
  !>    precipitation rate is to a standardized precipitation rate. This scheme assumes that a
  !>    fractional area of the grid box (at most 0.5) is affected by precipitation and that
  !>    this precipitation is due to a cloud of standardized cloud liquid water content. The
  !>    removal is constant throughout the column where the specific humidity is reduced.
  !> 2. `henry`: removal according to Henry's law, which states that the ratio of the
  !>    concentration in cloud water and the partial pressure in the interstitial air is a
  !>    constant. Here the units of Henry's constant are kg/L/Pa (normally they are M/L/Pa).
  !>
  !> To use it, add one of the following as a method for the tracer in the field table:
  !>
  !> * `"wet_deposition","henry","henry=XXX, dependence=YYY"`, where XXX is Henry's constant
  !>   for the tracer and YYY is the temperature dependence of Henry's constant;
  !> * `"wet_deposition","fraction","lslwc=XXX, convlwc=YYY"`, where XXX and YYY are the
  !>   liquid water contents [kg/m3] of a standard large-scale cloud (default 0.5e-3) and of
  !>   a standard convective cloud (default 2.0e-3).
  subroutine wet_deposition(n, T, pfull, phalf, rain, snow, qdt, tracer, tracer_dt, Time, cloud_param, is, js, dt)
    integer, intent(in)                 :: n, is, js
    !! `n`: tracer number; `is`, `js`: start indices of the arrays in the processor domain
    real, intent(in), dimension(:, :, :)  :: T, pfull, phalf, qdt, tracer
    !! `T`: temperature [K]; `pfull`, `phalf`: pressure at full and half levels [Pa]; `qdt`:
    !! tendency of the specific humidity due to the cloud parametrization [kg/kg/s]; `tracer`:
    !! tracer mixing ratio
    real, intent(in), dimension(:, :)    :: rain, snow  !! rain and snow reaching the surface [kg/m2/s]
    character(len=*), intent(in)         :: cloud_param
    !! cloud parametrization: convective (`'convect'`) or large-scale (`'lscale'`)
    type(time_type), intent(in)      :: Time  !! model time, for the diagnostic
    real, intent(out), dimension(:, :, :) :: tracer_dt  !! tendency of the tracer due to wet deposition
    real, intent(in)                    :: dt  !! time step [s]
!
    real, dimension(size(T, 1), size(T, 2), size(pfull, 3))   :: wsinku
    real, dimension(size(T, 1), size(T, 2)) :: Htemp, dz, washout, scav_factor, sum_wdep
    integer, dimension(size(T, 1), size(T, 2)) :: ktopcd, kendcd
    integer :: i, j, k, kd, flaglw
    real    :: Henry_constant, Henry_variable, inv298p15, clwc, wash, premin, prenow, hwtop
!real, dimension(size(rain,1),size(rain,2)) :: prenow,hwtop
    logical :: used, flag
    character(len=80) :: name, control, scheme, units
    tracer_dt = 0.0e+00
    ktopcd = 0
    kendcd = 0
    call get_tracer_names(MODEL_ATMOS, n, name, units=units)

    flag = query_method('wet_deposition', MODEL_ATMOS, n, name, control)
    if (.not. flag) return
    call get_wetdep_param(name, control, scheme, Henry_constant, Henry_variable)
    if (lowercase(scheme) == 'henry') then
! if units = MMR
! Henry_constant = [X](aq) / Px(g)
! where [X](aq) is the concentration of tracer X in precipitation
!       Px(g) is the partial pressure of the tracer in the air
! [X](aq) = Mixing ratio (MR) in cloud / qdt
! Px(g)   = MR (non cloud) * Pfull
!
! [X](aq)/Px = MR(incloud)/qdt /(Pfull MR non cloud) = H
! => MR(in cloud) = H * qdt * Pfull* MR(non cloud)
! MR (total) = MR(incloud) + MR(noncloud)
!            = MR(noncloud) * ( 1 + H*Pfull*qdt)
! MR(incloud) = H*Pfull*qdt * MR(total)/(1+H*Pfull*qdt)
! Fraction removed = MR(incloud)/MR(total) =
!  H*Pfull*qdt/(1+H*Pfull*qdt)
!
! if units = VMR
! Henry_constant = [X](aq) / Px(g)
! where [X](aq) is the concentration of tracer X in precipitation (mole/L)
!       Px(g) is the partial pressure of the tracer in the air
! [X](aq) = Volume mixing ratio (MR) in cloud / ( qdt * dt * MW_air )
! Px(g)   = VMR (non cloud) * Pfull
!
! [X](aq)/Px = VMR(incloud)/(qdt*dt*MW_air) /(Pfull VMR non cloud) = H
! => VMR(in cloud) = H * Pfull* VMR(non cloud) * (qdt * dt * MW_air)
! VMR (total) = VMR(incloud) + VMR(noncloud)
!             = VMR(noncloud) * ( 1 + H*Pfull*qdt*dt*MW_air)
! VMR(incloud) = H*Pfull*qdt*dt*MW_air * MR(total)/(1+H*Pfull*qdt*dt*MW_air)
! Fraction removed/dt = VMR(incloud)/VMR(total)/dt =
!  H*Pfull*qdt*MW_air/(1+H*Pfull*qdt*dt*MW_air)
!
      if (Henry_constant > 0) then
        inv298p15 = 1/298.15
        kd = size(T, 3)
        do k = 1, kd
          ! Calculate the temperature dependent part of Henry's constant
          ! exp( k *(1/T - 1/298.15))
          Htemp(:, :) = exp(Henry_variable*(1/T(:, :, k) - inv298p15))
          tracer_dt(:, :, k) = 0.0
          scav_factor(:, :) = 0.0
          where (qdt(:, :, k) < 0.0)
            !qdt is -ve so need to multiply by -1.0
            scav_factor(:, :) = -1.0*Henry_constant*Htemp*pfull(:, :, k)*qdt(:, :, k)*mw_air
            tracer_dt(:, :, k) = scav_factor(:, :)/(1 + scav_factor(:, :)*dt)
          end where
        end do
      end if
    end if

    if (lowercase(scheme) == 'fraction') then
      tracer_dt = 0.0
!-----------------------------------------------------------------------
!
!     Compute areal fractions experiencing wet deposition:
!
!     Set minimum precipitation rate below which no wet removal
!     occurs to 0.01 cm/day ie 1.16e-6 mm/sec (kg/m2/s)
      premin = 1.16e-6
!
!     Large scale cloud liquid water content (kg/m3)
!     and below cloud washout efficiency (cm-1):
      flaglw = parse(control, 'lslwc', clwc)
      if (flaglw == 0) clwc = 0.5e-3
      wash = 1.0
!
!     When convective adjustment occurs, use convective cloud liquid water content:
!
      if (trim(cloud_param) .eq. 'convect') then
        flaglw = parse(control, 'convlwc', clwc)
        if (flaglw == 0) clwc = 2.0e-3
        wash = 0.3
      end if
!
      do j = 1, size(rain, 2)
        do i = 1, size(rain, 1)
          tracer_dt(i, j, :) = 0.0
          washout(i, j) = 0.0
          prenow = rain(i, j) + snow(i, j)
          if (prenow .gt. premin) then
!
! Assume that the top of the cloud is where the highest model level
! specific humidity is reduced. And the the bottom of the cloud is the
! lowest model level where specific humidity is reduced.
!
            ktopcd(i, j) = 0
            do k = size(t, 3), 1, -1
              if (qdt(i, j, k) < 0.0) ktopcd(i, j) = k
            end do
            kendcd(i, j) = 0
            do k = 1, size(t, 3)
              if (qdt(i, j, k) < 0.0) kendcd(i, j) = k
            end do
!
!     Thickness of precipitating cloud deck:
!
            if (ktopcd(i, j) .gt. 1) then
              hwtop = 0.0
              do k = ktopcd(i, j), kendcd(i, j)
                hwtop = hwtop + (phalf(i, j, k + 1) - phalf(i, j, k))*rdgas*T(i, j, k)/grav/pfull(i, j, k)
              end do
              do k = ktopcd(i, j), kendcd(i, j)
!     Areal fraction affected by precip clouds (max = 0.5):
                tracer_dt(i, j, k) = prenow/(clwc*hwtop)
              end do
            end if

            washout(i, j) = prenow*wash
          end if
        end do
      end do
    end if

! Now multiply by the tracer mixing ratio to get the actual tendency.
    tracer_dt(:, :, :) = min(max(tracer_dt(:, :, :), 0.0e+00), 0.5)
    where (tracer > 0)
      tracer_dt = tracer_dt*tracer
    elsewhere
      tracer_dt = 0.0
    end where

    sum_wdep = 0.0
    do k = 1, size(tracer_dt, 3)
! delta z = dp/(rho * grav)
! delta z = RT/g*dp/p
! tracer(kgtracer/kgair) * dz(m)* rho(kgair/m3) = kgtracer/m2
! so rho drops out of the equation
      if (units(1:3) == 'mmr') then
        sum_wdep = sum_wdep + tracer_dt(:, :, k)*(phalf(:, :, k + 1) - phalf(:, :, k))/grav
      end if
      if (units(1:3) == 'vmr') then
        sum_wdep = sum_wdep + tracer_dt(:, :, k)*(phalf(:, :, k + 1) - phalf(:, :, k))/(grav*mw_air)
      end if
    end do

    if (trim(cloud_param) .eq. 'lscale') then
      if (id_tracer_wdep_ls(n) > 0) then
        used = send_data(id_tracer_wdep_ls(n), sum_wdep, Time, &
                         is_in=is, js_in=js)
      end if
    end if
    if (trim(cloud_param) .eq. 'convect') then
      if (id_tracer_wdep_cv(n) > 0) then
        used = send_data(id_tracer_wdep_cv(n), sum_wdep, Time, &
                         is_in=is, js_in=js)
      end if
    end if
  end subroutine wet_deposition
!
!#######################################################################
!
  subroutine get_drydep_param(text_in_scheme, text_in_param, scheme, land_dry_dep_vel, sea_dry_dep_vel)
!
! Subroutine to initialiize the parameters for the dry deposition scheme.
! If the dry dep scheme is 'fixed' then the dry_deposition velocity value
! has to be set.
! If the dry dep scheme is 'wind_driven' then the dry_deposition
! velocity value will be calculated. So set to a dummy value of 0.0
! INTENT IN
!  text_in_scheme   : The text that has been parsed from tracer table as
!                     the dry deposition scheme to be used.
!  text_in_param    : The parameters that are associated with the dry
!                     deposition scheme.
! INTENT OUT
!  scheme           : The scheme that is being used.
!  land_dry_dep_vel : Dry deposition velocity over the land
!  sea_dry_dep_vel  : Dry deposition velocity over the sea
!
    character(len=*), intent(in)    :: text_in_scheme, text_in_param
    character(len=*), intent(out)   :: scheme
    real, intent(out)               :: land_dry_dep_vel, sea_dry_dep_vel

    integer :: m, m1, n, lentext, flag
    character(len=32) :: dummy

!Default
    scheme = 'None'
    land_dry_dep_vel = 0.0
    sea_dry_dep_vel = 0.0

    if (lowercase(trim(text_in_scheme(1:4))) .eq. 'wind') then
      scheme = 'Wind_driven'
      land_dry_dep_vel = 0.0
      sea_dry_dep_vel = 0.0
    end if

    if (lowercase(trim(text_in_scheme(1:5))) .eq. 'fixed') then
      scheme = 'fixed'
      flag = parse(text_in_param, 'land', land_dry_dep_vel)
      flag = parse(text_in_param, 'sea', sea_dry_dep_vel)
    end if

    if (lowercase(trim(text_in_scheme(1:4))) .eq. 'file') then
      scheme = 'file'
      land_dry_dep_vel = 0.
      sea_dry_dep_vel = 0.
    end if

  end subroutine get_drydep_param
!
!#######################################################################
!
  subroutine get_wetdep_param(text_in_scheme, text_in_param, scheme, henry_constant, henry_temp)
!
! Routine to initialize the parameters for the wet deposition scheme.
! INTENT IN
!  text_in_scheme : Text read from the tracer table which provides information on which
!                   wet deposition scheme to use.
!  text_in_param  : Parameters associated with the wet deposition scheme. These will be
!                   parsed in this routine.
! INTENT OUT
!  scheme         : Wet deposition scheme to use.
!                   Choices are None, Fraction and Henry
!  henry_constant : Henry's constant for the tracer (see wet_deposition for explanation of Henry's Law)
!  henry_temp     : The temperature dependence of the Henry's Law constant.
!
!
    character(len=*), intent(in)    :: text_in_scheme, text_in_param
    character(len=*), intent(out)   :: scheme
    real, intent(out)               :: henry_constant, henry_temp

    integer :: m, m1, n, lentext, flag
    character(len=32) :: dummy

!Default
    scheme = 'None'
    henry_constant = 0.0
    henry_temp = 0.0

    if (trim(lowercase(text_in_scheme(1:8))) .eq. 'fraction') then
      scheme = 'Fraction'
      henry_constant = 0.0
      henry_temp = 0.0
    end if

    if (trim(lowercase(text_in_scheme(1:5))) .eq. 'henry') then
      scheme = 'Henry'
      flag = parse(text_in_param, 'henry', henry_constant)
      flag = parse(text_in_param, 'dependence', henry_temp)
    end if

  end subroutine get_wetdep_param
!
!#######################################################################
!
  !> Interpolates an emission field (or any 2D field) of arbitrary resolution to the model
  !> grid, and returns the part of the global field on the local processor.
  subroutine interp_emiss(global_source, start_lon, start_lat, &
                          lon_resol, lat_resol, data_out)
    real, intent(in)  :: global_source(:, :)  !! global emission field
    real, intent(in)  :: start_lon, start_lat, lon_resol, lat_resol
    !! `start_lon`, `start_lat`: western and southern boundaries of the global field [rad];
    !! `lon_resol`, `lat_resol`: longitudinal and latitudinal resolution of the field [rad]
    real, intent(out) :: data_out(:, :)  !! interpolated field on the local processor

    real :: modydeg, modxdeg, tpi
    integer :: i, j, nlon_in, nlat_in
    real :: blon_in(size(global_source, 1) + 1)
    real :: blat_in(size(global_source, 2) + 1)
! Set up the global surface boundary condition longitude-latitude boundary values

    tpi = 2.*PI
    nlon_in = size(global_source, 1)
    nlat_in = size(global_source, 2)
! For some reason the input longitude needs to be incremented by 180 degrees.
    do i = 1, nlon_in + 1
      blon_in(i) = start_lon + float(i - 1)*lon_resol + PI
    end do
    if (abs(blon_in(nlon_in + 1) - blon_in(1) - tpi) < epsilon(blon_in)) &
      blon_in(nlon_in + 1) = blon_in(1) + tpi

    do j = 2, nlat_in
      blat_in(j) = start_lat + float(j - 1)*lat_resol
    end do
    blat_in(1) = -0.5*PI
    blat_in(nlat_in + 1) = 0.5*PI

! Now interpolate the global data to the model resolution
    call horiz_interp(global_source, blon_in, blat_in, &
                      blon_out, blat_out, data_out)

  end subroutine interp_emiss
!
!######################################################################
  !> Terminates the tracer utilities module.
  subroutine atmos_tracer_utilities_end

    deallocate (blon_out, blat_out)
    module_is_initialized = .false.

  end subroutine atmos_tracer_utilities_end

! ######################################################################
!

end module atmos_tracer_utilities_mod

