!> Black and organic carbon aerosol tracers, after Cooke et al. (1999, 2002).
!>
!> In its present implementation the black and organic carbon tracers are from the
!> combustion of fossil fuel. The annual-mean emissions are read from `INPUT/r30.bc.ann` and
!> `INPUT/r30.oc.ann` (the datasets derived in Cooke et al. 1999, on a 3.6 x 3 degree R30
!> grid) and interpolated to the model grid. The tracers are `bcphob`, `bcphil`, `ocphob` and
!> `ocphil` in the `field_table`.
!>
!> While the code should provide insights into the carbonaceous aerosol cycle, it is provided
!> more as an example of how to implement a tracer module. The parameters should be checked
!> and set to the values of previous works if a user wishes to reproduce those works.
!>
!> References:
!>
!> * Cooke, W. F., and J. J. N. Wilson, 1996: A global black carbon aerosol model.
!>   J. Geophys. Res., 101, 19395-19409.
!> * Cooke, W. F., C. Liousse, H. Cachier, and J. Feichter, 1999: Construction of a 1 x 1
!>   fossil fuel emission dataset for carbonaceous aerosol and implementation and radiative
!>   impact in the ECHAM-4 model. J. Geophys. Res., 104, 22137-22162.
!> * Cooke, W. F., V. Ramaswamy, and P. Kasibhatla, 2002: A GCM study of the global
!>   carbonaceous aerosol distribution. J. Geophys. Res., 107.
!>
!> Original authors: William Cooke.
module atmos_carbon_aerosol_mod
  use fms_mod, only: &
    mpp_pe, &
    mpp_root_pe, &
    stdlog, &
    write_version_number
  use time_manager_mod, only: time_type
  use diag_manager_mod, only: send_data, &
                              register_diag_field, &
                              register_static_field
  use tracer_manager_mod, only: get_tracer_index, &
                                set_tracer_atts
  use field_manager_mod, only: MODEL_ATMOS
  use atmos_tracer_utilities_mod, only: interp_emiss
  use constants_mod, only: PI

  implicit none
  private
!-----------------------------------------------------------------------
!----- interfaces -------

  public atmos_blackc_sourcesink, &
    atmos_organic_sourcesink, &
    atmos_carbon_aerosol_init, &
    atmos_carbon_aerosol_end

!-----------------------------------------------------------------------
!----------- namelist -------------------
!-----------------------------------------------------------------------
!
!  When initializing additional tracers, the user needs to make the
!  following changes.
!
!  Add an integer variable below for each additional tracer. This should
!  be initialized to zero.
!
!  Add id_tracername for each additional tracer. These are used in
!  initializing and outputting the tracer fields.
!
!-----------------------------------------------------------------------

! tracer number for radon
  integer :: nbcphobic = 0
  integer :: nbcphilic = 0
  integer :: nocphobic = 0
  integer :: nocphilic = 0

!--- identification numbers for  diagnostic fields and axes ----

  integer :: id_emissoc, id_emissbc

!--- Arrays to help calculate tracer sources/sinks ---
  real, allocatable, dimension(:, :) :: bcsource, ocsource

  character(len=6), parameter :: module_name = 'tracer'

  logical :: module_is_initialized = .false.
  logical :: used

!---- version number -----
  character(len=128) :: version = '$Id: atmos_carbon_aerosol.f90,v 11.0 2004/09/28 19:26:31 fms Exp $'
  character(len=128) :: tagname = '$Name: lima $'
!-----------------------------------------------------------------------

contains

!#######################################################################

  !> Computes the tendencies of hydrophobic and hydrophilic black carbon due to emission
  !> and transformation.
  !>
  !> The hydrophobic aerosol has sources from emissions (80%) and sinks from dry deposition
  !> and transformation into hydrophilic aerosol. The hydrophilic aerosol also has emission
  !> sources (20%) and has sinks of wet and dry deposition. The emissions go into the lowest
  !> model level. The deposition is computed in `atmos_tracer_utilities_mod`, not here. The
  !> transformation time used here is 1 day, which corresponds to an e-folding time of
  !> 1.44 days.
  subroutine atmos_blackc_sourcesink(lon, lat, land, pwt, &
                                     black_cphob, black_cphob_dt, &
                                     black_cphil, black_cphil_dt, &
                                     Time, is, ie, js, je, kbot)

!-----------------------------------------------------------------------
    real, intent(in), dimension(:, :)   :: lon, lat
    !! longitude and latitude of the centres of the grid cells [rad]
    real, intent(in), dimension(:, :)   :: land  !! land fraction
    real, intent(in), dimension(:, :, :) :: pwt, black_cphob, black_cphil
    !! `pwt`: pressure weight dp/grav [kg/m2]; `black_cphob`, `black_cphil`: hydrophobic and
    !! hydrophilic black carbon mixing ratio
    real, intent(out), dimension(:, :, :) :: black_cphob_dt, black_cphil_dt
    !! tendencies of the hydrophobic and hydrophilic black carbon mixing ratio
    type(time_type), intent(in)            :: Time  !! model time
    integer, intent(in)                    :: is, ie, js, je  !! local domain boundaries
    integer, intent(in), dimension(:, :), optional :: kbot  !! index of the lowest model level above the surface
!-----------------------------------------------------------------------
    real, dimension(size(black_cphob, 1), size(black_cphob, 2), size(black_cphob, 3)) :: &
      sourcephob, sinkphob, sourcephil, sinkphil
    real dtr
    integer i, j, kb, id, jd, kd, lat1
!-----------------------------------------------------------------------

    id = size(black_cphob, 1); jd = size(black_cphob, 2); kd = size(black_cphob, 3)

    dtr = PI/180.

!----------- compute black carbon source ------------

    sourcephob = 0.0
    sourcephil = 0.0

    do j = 1, jd
      sourcephob(:, j, kd) = 0.8*bcsource(:, j + js - 1)/pwt(:, j, kd)
      sourcephil(:, j, kd) = 0.2*bcsource(:, j + js - 1)/pwt(:, j, kd) + &
                             8.038e-6*black_cphob(:, j, kd)
    end do

!------- compute black carbon phobic sink --------------
!
!  BCphob has a half-life time of 1.0days
!   (corresponds to an e-folding time of 1.44 days)
!
!  sink = 1./(86400.*1.44) = 8.023e-6
!

    where (black_cphob(:, :, :) >= 0.0)
      sinkphob(:, :, :) = -8.038e-6*black_cphob(:, :, :)
    elsewhere
      sinkphob(:, :, :) = 0.0
    end where

    sinkphil(:, :, :) = 0.0

!------- tendency ------------------

    black_cphob_dt = sourcephob + sinkphob
    black_cphil_dt = sourcephil + sinkphil

!-----------------------------------------------------------------------

  end subroutine atmos_blackc_sourcesink

!#######################################################################

  !> Computes the tendency of organic carbon due to emission and transformation.
  !>
  !> The hydrophobic aerosol has sources from emissions and sinks from dry deposition and
  !> transformation into hydrophilic aerosol. The hydrophilic aerosol also has emission
  !> sources and has sinks of wet and dry deposition. The emissions go into the lowest model
  !> level. The deposition is computed in `atmos_tracer_utilities_mod`, not here. The
  !> transformation time used here is 2 days, which corresponds to an e-folding time of
  !> 2.88 days.
  subroutine atmos_organic_sourcesink(lon, lat, land, pwt, organic_carbon, organic_carbon_dt, &
                                      Time, is, ie, js, je, kbot)

!-----------------------------------------------------------------------
    real, intent(in), dimension(:, :)   :: lon, lat
    !! longitude and latitude of the centres of the grid cells [rad]
    real, intent(in), dimension(:, :)   :: land  !! land fraction
    real, intent(in), dimension(:, :, :) :: pwt, organic_carbon
    !! `pwt`: pressure weight dp/grav [kg/m2]; `organic_carbon`: organic carbon mixing ratio
    real, intent(out), dimension(:, :, :) :: organic_carbon_dt  !! tendency of the organic carbon mixing ratio
    type(time_type), intent(in) :: Time  !! model time
    integer, intent(in)                    :: is, ie, js, je  !! local domain boundaries
    integer, intent(in), dimension(:, :), optional :: kbot  !! index of the lowest model level above the surface
!-----------------------------------------------------------------------
    real, dimension(size(organic_carbon, 1), size(organic_carbon, 2), size(organic_carbon, 3)) :: &
      source, sink
    real dtr
    integer i, j, kb, id, jd, kd, lat1
!-----------------------------------------------------------------------

    id = size(organic_carbon, 1); jd = size(organic_carbon, 2); kd = size(organic_carbon, 3)

    dtr = PI/180.

!----------- compute organic carbon source ------------

    source = 0.0

    if (present(kbot)) then
      do j = 1, jd
      do i = 1, id
        kb = kbot(i, j)
        source(i, j, kb) = ocsource(i, j + js - 1)/pwt(i, j, kb)
      end do
      end do
    else
      do j = 1, je - js + 1
        source(:, j, kd) = ocsource(:, j + js - 1)/pwt(:, j, kd)
      end do
    end if

!------- compute organic carbon sink --------------
!
!  OCphob has a half-life time of 2.0days
!   (corresponds to an e-folding time of 2.88 days)
!
!  sink = 1./(86400.*2.88) = 4.019e-6
!

    where (organic_carbon(:, :, :) >= 0.0)
      sink(:, :, :) = -4.019e-6*organic_carbon(:, :, :)
    elsewhere
      sink(:, :, :) = 0.0
    end where

!------- tendency ------------------

    organic_carbon_dt = source + sink

!-----------------------------------------------------------------------

  end subroutine atmos_organic_sourcesink

!#######################################################################

  !> Initializes the carbon aerosol module: finds the indices of the carbonaceous aerosol
  !> tracers, registers the emission fields as diagnostics and reads the emissions.
  subroutine atmos_carbon_aerosol_init(lonb, latb, r, axes, Time, mask)

!-----------------------------------------------------------------------
    real, dimension(:), intent(in) :: lonb, latb  !! longitudes and latitudes of the cell corners [rad]
    real, intent(inout), dimension(:, :, :, :) :: r  !! tracer fields (nlon, nlat, nlev, ntrace)
    integer, intent(in)                        :: axes(4)  !! diagnostic axes (lon, lat, pfull, phalf)
    type(time_type), intent(in)                        :: Time  !! model time
    real, intent(in), dimension(:, :, :), optional :: mask
    !! 1. above the ground, 0. below (nlon, nlat, nlev)

    integer :: n

    if (module_is_initialized) return

!----- set initial value of carbon ------------

    n = get_tracer_index(MODEL_ATMOS, 'bcphob')
    if (n > 0) then
      nbcphobic = n
      call set_tracer_atts(MODEL_ATMOS, 'bcphob', 'hphobic_bc', 'g/g')
      if (nbcphobic > 0 .and. mpp_pe() == mpp_root_pe()) write (*, 30) 'Hydrophobic BC', nbcphobic
      if (nbcphobic > 0 .and. mpp_pe() == mpp_root_pe()) write (stdlog(), 30) 'Hydrophobic BC', nbcphobic
    end if

    n = get_tracer_index(MODEL_ATMOS, 'bcphil')
    if (n > 0) then
      nbcphilic = n
      call set_tracer_atts(MODEL_ATMOS, 'bcphil', 'hphilic_bc', 'g/g')
      if (nbcphilic > 0 .and. mpp_pe() == mpp_root_pe()) write (*, 30) 'Hydrophilic BC', nbcphilic
      if (nbcphilic > 0 .and. mpp_pe() == mpp_root_pe()) write (stdlog(), 30) 'Hydrophilic BC', nbcphilic
    end if

    n = get_tracer_index(MODEL_ATMOS, 'ocphob')
    if (n > 0) then
      nocphobic = n
      call set_tracer_atts(MODEL_ATMOS, 'ocphob', 'hphobic_oc', 'g/g')
      if (nocphobic > 0 .and. mpp_pe() == mpp_root_pe()) write (*, 30) 'Hydrophobic OC', nocphobic
      if (nocphobic > 0 .and. mpp_pe() == mpp_root_pe()) write (stdlog(), 30) 'Hydrophobic OC', nocphobic
    end if

    n = get_tracer_index(MODEL_ATMOS, 'ocphil')
    if (n > 0) then
      nocphilic = n
      call set_tracer_atts(MODEL_ATMOS, 'ocphil', 'hphilic_oc', 'g/g')
      if (nocphilic > 0 .and. mpp_pe() == mpp_root_pe()) write (*, 30) 'Hydrophilic OC', nocphilic
      if (nocphilic > 0 .and. mpp_pe() == mpp_root_pe()) write (stdlog(), 30) 'Hydrophilic OC', nocphilic
    end if

30  format(A, ' was initialized as tracer number ', i2)
    !Read in emission files
!
    id_emissbc = register_static_field('tracers', &
                                       'bcemiss', axes(1:2), &
                                       'black carbon emission', 'g/m2/s')
    id_emissoc = register_static_field('tracers', &
                                       'ocemiss', axes(1:2), &
                                       'organic carbon emission', 'g/m2/s')
!
    allocate (bcsource(size(lonb(:)) - 1, size(latb(:)) - 1))
    allocate (ocsource(size(lonb(:)) - 1, size(latb(:)) - 1))
    call tracer_input(lonb, latb, Time)

    call write_version_number(version, tagname)
    module_is_initialized = .true.

!-----------------------------------------------------------------------

  end subroutine atmos_carbon_aerosol_init

  !> Terminates the carbon aerosol module.
  subroutine atmos_carbon_aerosol_end

    module_is_initialized = .false.

  end subroutine atmos_carbon_aerosol_end

!#######################################################################
  !> Reads the black and organic carbon emissions and interpolates them to the model grid.
  subroutine tracer_input(lonb, latb, Time)
    real, dimension(:), intent(in) :: lonb, latb
    type(time_type), intent(in) :: Time

    integer      :: i, j, unit, io
    real         :: emiss
    real         :: dtr, deg_90, deg_180, deg3p6, deg3!, modxdeg, modydeg
    real         :: ZCARBONSEASON(12)
    real         :: bcsource1(100, 60)
    logical :: opened
!
! This is the Rotty seaonality for fossil fuel emissions of sulfate.
!
    data ZCARBONSEASON/1.146, 1.139, 1.081, 0.995, 0.916, 0.920, 0.910, &
      0.907, 0.934, 0.962, 1.019, 1.072/
!
    dtr = PI/180.
    deg_90 = -90.*dtr; deg_180 = -180.*dtr
    ! -90 and -180 degrees are the southwest boundaries of the
    ! emission  field you are reading in.
    deg3p6 = 3.6*dtr; deg3 = 3.*dtr; 
    ! 3.6 degrees longitude and 3 degree latitude is the resolution
    ! of the r30 emission data that I used in SKYHI.

! initialise the BC phobic and philic
! read in the emission sources here.

    do unit = 30, 100
      inquire (unit=unit, opened=opened)
      if (.not. opened) exit
    end do
    open (unit, file='INPUT/r30.bc.ann', form='formatted', action='read')
    do io = 1, 6000
      read (unit, FMT=1968, end=11) i, j, emiss
      bcsource1(i, j) = emiss
    end do
11  close (unit)
1968 format(2i3, e11.4)
1969 format(2i3, f11.3)
! Interpolate the R30 emission field to the resolution of the model.
    call interp_emiss(bcsource1, 0.0, deg_90, deg3p6, deg3, &
                      bcsource)

    if (mpp_pe() == mpp_root_pe()) write (*, *) 'Reading OC emissions'
!Now let's do the OC
!
    bcsource1 = 0.0e+00
    do unit = 30, 100
      inquire (unit=unit, opened=opened)
      if (.not. opened) exit
    end do
    open (unit, file='INPUT/r30.oc.ann', form='formatted', action='read')
    do io = 1, 6000
      read (unit, FMT=1968, end=13) i, j, emiss
      bcsource1(i, j) = emiss
    end do
13  close (unit)

! Interpolate the R30 emission field to the resolution of the model.
    call interp_emiss(bcsource1, 0.0, deg_90, deg3p6, deg3, &
                      ocsource)

! Send the emission data to the diag_manager for output.
    if (id_emissbc > 0) &
      used = send_data(id_emissbc, bcsource)
    if (id_emissoc > 0) &
      used = send_data(id_emissoc, ocsource)

  end subroutine tracer_input

end module atmos_carbon_aerosol_mod

