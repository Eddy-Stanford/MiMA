
!> Prescribed ocean heat fluxes (Q-fluxes) for the slab ocean.
!>
!> `qflux` adds the zonally symmetric meridional Q-flux of Merlis et al. (2013,
!> Part II); `warmpool` adds zonally asymmetric fluxes: a tropical warm pool and, depending
!> on `warmpool_localization_choice`, regional patterns such as the Gulf Stream, the
!> Kuroshio and the tropical Atlantic, following Garfinkel et al. (2020). Both are called
!> by `simple_surface` with `do_qflux` and `do_warmpool` in `simple_surface_nml`.
!>
!> Namelist: `qflux_nml`
!> ([namelist reference](https://eddy-stanford.github.io/MiMA/Parameters/#qflux_nml)).
!>
!> References:
!>
!> * Merlis, T. M., T. Schneider, S. Bordoni, and I. Eisenman, 2013: Hadley circulation
!>   response to orbital precession. Part II: Subtropical continent. J. Climate, 26,
!>   754-771.
!> * Garfinkel, C. I., I. White, E. P. Gerber, M. Jucker, and M. Erez, 2020:
!>   The building blocks of Northern Hemisphere wintertime stationary waves.
!>   J. Climate, 33, 5611-5633, https://doi.org/10.1175/JCLI-D-19-0181.1.
module qflux_mod

  use constants_mod, only: pi
  use fms_mod, only: input_nml_file, check_nml_error, &
                     error_mesg, FATAL

  implicit none

  real ::    qflux_amp = 26., & !! [W/m2] amplitude of the meridional Q-flux
          qflux_width = 16., & !! [deg] half-width of the meridional Q-flux
          warmpool_amp = 18., & !! [W/m2] amplitude of the warm pool
          warmpool_width = 35., & !! [deg] latitudinal width of the warm pool
          warmpool_centr = 0., & !! [deg] central latitude of the warm pool
          warmpool_phase = 140.     !! [deg] longitude phase of the warm pool

  integer :: gulf_k = 4       !! zonal wave number of the Gulf Stream perturbation (choice 2; with choice 3
  !! only in a North Atlantic term near 67N)

  real :: warmpool_k = 1.66666, & !! zonal wave number of the warm pool
          gulf_phase = 310., & !! [deg] longitude phase of the Gulf Stream perturbation (choice 2; with
          !! choice 3 only in a North Atlantic term near 67N)
          gulf_amp = 70., & !! [W/m2] Gulf Stream amplitude (choices 2 and 3; with choice 3 it scales a fixed,
          !! localized Gulf Stream pattern, and the tropical Atlantic term is only applied if `gulf_amp` > 0)
          kuroshio_amp = 40., &  !! [W/m2] Kuroshio amplitude (choices 2 and 3)
          trop_atlantic_amp = 50., &  !! [W/m2] tropical Atlantic amplitude (choices 2 and 3)
          north_sea_heat = 0., &
          !! [1] factor on `gulf_amp` for moving heat from Canada to the North Sea (choice 2 only)
          Pac_ITCZextra = 0., & !! [W/m2] extra flux in the tropical South Pacific (strengthens the local ITCZ; choice 3 only)
          Pac_SPCZextra = 0., & !! [W/m2] extra flux in the subtropical Pacific (modulates the SPCZ; choice 3 only)
          Africaextra = 0., &  !! [W/m2] extra flux near the Agulhas current (choice 3 only)
          Sampeextra = 0., & !! [W/m2] extra flux off South America (choice 3 only)
          Hawaiiextra = 30.0 !! [W/m2] extra flux near Hawaii (choice 3 only)

  integer :: warmpool_localization_choice = 3
  !! 1: cosine in longitude; 2: cosine restricted to the Indo-Pacific, plus Gulf Stream, Kuroshio and
  !! tropical Atlantic terms; 3: the localized patterns of Garfinkel et al. (2020). Which of the
  !! regional amplitudes below are used depends on this choice (see `qflux.f90`).
  logical :: qflux_initialized = .false.

  namelist /qflux_nml/ qflux_amp, qflux_width, &
    warmpool_amp, warmpool_width, warmpool_centr, &
    warmpool_k, warmpool_phase, warmpool_localization_choice, &
    gulf_k, gulf_phase, gulf_amp, kuroshio_amp, trop_atlantic_amp, &
    north_sea_heat, Pac_ITCZextra, Sampeextra, &
    Pac_SPCZextra, Africaextra, Hawaiiextra

  private

  public :: qflux_init, qflux, warmpool

contains

!########################################################
  !> Reads `qflux_nml`.
  subroutine qflux_init
    implicit none
    integer :: unit, ierr, io

    read (input_nml_file, nml=qflux_nml, iostat=io)
    ierr = check_nml_error(io, 'qflux_nml')

    qflux_initialized = .true.

  end subroutine qflux_init
!########################################################

  !> Subtracts the meridional Q-flux of Merlis et al. (2013, Part II) from `flux`.
  subroutine qflux(latb, flux)
    implicit none
    real, dimension(:), intent(in)    :: latb   !! latitudes of the cell boundaries [rad]
    real, dimension(:, :), intent(inout) :: flux   !! total ocean heat flux [W/m2]
!
    integer j
    real lat, coslat

    if (.not. qflux_initialized) then
      call error_mesg('qflux', 'qflux module not initialized', FATAL)
    end if

    do j = 1, size(latb) - 1
      lat = 0.5*(latb(j + 1) + latb(j))
      coslat = cos(lat)
      lat = lat*180./pi
      flux(:, j) = flux(:, j) - qflux_amp*(1 - 2.*lat**2/qflux_width**2)* &
                   exp(-((lat)**2/(qflux_width)**2))/coslat
    end do

  end subroutine qflux

!########################################################

  !> Adds the zonally asymmetric Q-fluxes (warm pool and regional patterns) to `flux`.
  subroutine warmpool(lonb, latb, flux)
    implicit none
    real, dimension(:), intent(in)   :: lonb, latb  !! longitudes and latitudes of the cell boundaries [rad]
    real, dimension(:, :), intent(inout):: flux       !! total ocean heat flux [W/m2]
!
    integer i, j
    real lon, lat, piphase, pigulfphase, latgulf, latgreen, latorig, africaamp

    africaamp = trop_atlantic_amp*2/4.
!    piphase = warmpool_phase/pi
    piphase = warmpool_phase*pi/180  !modified by cig, nov 15 2017
    pigulfphase = gulf_phase*pi/180  !modified by cig, nov 15 2017
    do j = 1, size(latb) - 1
      lat = 0.5*(latb(j + 1) + latb(j))*180./pi
      latgulf = (lat - 37.)/10.
      latgreen = (lat - 67.)/10.

      latorig = lat
      lat = (lat - warmpool_centr)/warmpool_width

      do i = 1, size(lonb) - 1
        lon = 0.5*(lonb(i + 1) + lonb(i))
        if (abs(lat) .le. 1.0) then
          if (warmpool_localization_choice == 1) then
            !modified by cig, nov 15 2017, note that I use 4th power not 2nd as in MJ to better match Pac warm pool
            flux(i, j) = flux(i, j) &
                 &+ (1.-lat**4.)*warmpool_amp*cos(warmpool_k*(lon - piphase))
          elseif (warmpool_localization_choice == 2 .or. warmpool_localization_choice == 3) then  !assumes k=5/3 for warmpool,
            !modified by cig, nov 15 2017
            if (lon .ge. (warmpool_phase - 54)*pi/180. .and. lon .le. (warmpool_phase + 162)*pi/180.) then
              flux(i, j) = flux(i, j) &
                 &+ (1.-lat**4.)*warmpool_amp*cos(warmpool_k*(lon - piphase))
            end if
            !modified by cig, mar 28 2019
            if (lon .ge. (warmpool_phase + 117)*pi/180. .and. lon .le. (warmpool_phase + 162)*pi/180. &
              & .and. warmpool_localization_choice == 3) then
              flux(i, j) = flux(i, j) &
                 &+ (1.-lat**4.)*warmpool_amp*sin(8*(lon - piphase - 139.5*pi/180.))
            end if
          end if
          !modified by cig, may 13 2019
          if (lon .ge. (warmpool_phase - 130)*pi/180. .and. lon .le. (warmpool_phase - 58)*pi/180. &
            & .and. warmpool_localization_choice == 3) then
            flux(i, j) = flux(i, j) &
               &+ (1.-lat**2.)*africaamp*cos(5*(lon - (piphase - 112*pi/180)))
          end if
          !modified by cig, april 30 2018
          if ((lon .ge. (gulf_phase - 22)*pi/180. .or. lon .le. (gulf_phase + 68 - 360)*pi/180.) &
            & .and. warmpool_localization_choice == 2) then
            flux(i, j) = flux(i, j) &
               &+ (1.-lat**4.)*trop_atlantic_amp*cos(gulf_k*(lon - pigulfphase))
          end if
        end if

        if (abs(latgulf) .le. 1.0 .and. warmpool_localization_choice == 2) then !add Kuroshio and Gulf streams

          if (lon .ge. (gulf_phase - 42)*pi/180. .and. lon .lt. (gulf_phase + 48)*pi/180.) then !modified by cig, june 3 2018
            flux(i, j) = flux(i, j) &
               &+ (1.-latgulf**4.)*gulf_amp*cos(gulf_k*(lon - (pigulfphase - 19.5*pi/180.)))
          end if
          !modified by cig, Nov 18 2018 to localize gulfstream more over ocean
          if (lon .ge. (gulf_phase - 42)*pi/180. .and. lon .lt. (gulf_phase + 3)*pi/180.) then
            flux(i, j) = flux(i, j) &
               &+ 0.535*(1.-latgulf**4.)*gulf_amp*sin(gulf_k*2*(lon - (pigulfphase - 19.5*pi/180.)))
          end if

          !modified by cig, june 3 2018, use k=4 for kuroshio
          if (lon .ge. (warmpool_phase - 17.5)*pi/180. .and. lon .lt. (warmpool_phase + 72.5)*pi/180.) then
            flux(i, j) = flux(i, j) &
               &+ (1.-latgulf**2.)*kuroshio_amp*cos(gulf_k*(lon - (piphase + 5*pi/180))) !modified by cig, mar 29 2019,
          end if
          !modified by cig, jan 2 2019, shift cooling kuroshio east k=6
          if (lon .ge. (warmpool_phase + 30)*pi/180. .and. lon .lt. (warmpool_phase + 90)*pi/180.) then
            flux(i, j) = flux(i, j) &
               &- (1.-latgulf**2.)*0.65*kuroshio_amp*cos((2*gulf_k - 2)*(lon - (piphase + 15*pi/180)))!modified by cig, mar 29 2019
          end if
        end if

        if (abs(latgreen) .le. 1.0 .and. warmpool_localization_choice == 2 .and. north_sea_heat .gt. 0.001) then

          !modified by cig, Nov 22 2018 to add opposite of Gulfstream further polewrd
          if (lon .ge. (gulf_phase - 52)*pi/180. .or. lon .lt. (gulf_phase + 68 - 360)*pi/180.) then
            flux(i, j) = flux(i, j) &
               & + north_sea_heat*(1.-latgreen**4.)*gulf_amp*cos((gulf_k - 1)*(lon - pigulfphase - 38*pi/180.))

          end if
          !modified by cig, Nov 22 2018 to smear out the cooling over Northern Canada over broad region
          if (lon .ge. (gulf_phase - 67)*pi/180. .and. lon .lt. (gulf_phase - 7)*pi/180.) then
            flux(i, j) = flux(i, j) &
               & + north_sea_heat*0.25*(1.-latgreen**4.)*gulf_amp*cos(2*(gulf_k - 1)*(lon - pigulfphase + 22*pi/180.))

          end if
        end if

        !add Kuroshio
        if (warmpool_localization_choice == 3 .and. kuroshio_amp .gt. 0.001 .and. latorig .le. 47. .and. latorig .ge. 5.) then
          !modified by cig, may 13 2019, Pacific sector is exp
          if (lon .ge. (warmpool_phase - 30.)*pi/180. .and. lon .lt. (warmpool_phase + 130.)*pi/180.) then
            flux(i, j) = flux(i, j) &
               &- kuroshio_amp*59.4/100.*exp(-(lon*180./pi + latorig - 268.)**2./(2*49.)) &
                 & *exp(-(lon*180./pi - latorig - 215.)**2./(2*625.))  &
               &+ kuroshio_amp*exp(-(lon*180./pi - 3*latorig - 45.)**2./(2*100.))*exp(-(lon*180./pi + latorig - 170.)**2./(2*400.))
          end if

        end if

        !add more Kurishio at expense of Southeast Asia  modified by cig, may 14 2019
        if (warmpool_localization_choice == 3 .and. warmpool_amp .gt. 0.001 .and. lon .ge. 70.*pi/180. .and. latorig .le. 60. &
          & .and. latorig .ge. -10. .and. lon .lt. 240.*pi/180.) then
          if (kuroshio_amp .gt. 0.001) then
            flux(i, j) = flux(i, j) &
            &- 27.60*exp(-(lon*180./pi - 140)**2./(2*1521))*(exp(-(latorig - 19.7)**2./(2*49)))   &
            &- (5.2)*exp(-(lon*180./pi - 140)**2./(2*64))*(exp(-(latorig - 20.)**2./(2*16)))   &
            &+ 35.4*exp(-(lon*180./pi - 160)**2./(2*400))*exp(-(latorig - 35)**2./(2*36)) &
               &+ (49.5 - Hawaiiextra*.0228)*exp(-(lon*180./pi - 3*latorig - 45.-Hawaiiextra/5)**2./(2*100.)) &
                 & *exp(-(lon*180./pi + latorig - 160.)**2./(2*400.)) &
            &+ 22.9*exp(-(lon*180./pi - 90)**2./(2*144))*exp(-(latorig - 0)**2./(2*25))
          else
            flux(i, j) = flux(i, j) &
                 &- 5.3*exp(-(lon*180./pi - 140)**2./(2*1521))*(exp(-(latorig - 19.7)**2./(2*49)))   &
                 &+ 22.9*exp(-(lon*180./pi - 90)**2./(2*144))*exp(-(latorig - 0)**2./(2*25))

          end if
        end if

        !add south pacific and SPCZ at expense of further cold tongue  modified by cig, may 28 2019, also includes Australia bit
        !from Kuroshio above
        if (warmpool_localization_choice == 3 .and. warmpool_amp .gt. 0.001 .and. latorig .le. 24. .and. latorig .ge. -78. &
          & .and. lon .ge. 129.*pi/180. .and. lon .lt. 290.*pi/180.) then
          flux(i, j) = flux(i, j) &
              &- (50.-Pac_SPCZextra*.28 - Hawaiiextra*.6)*exp(-(lon*180./pi - 270)**2./(2*81))*(exp(-(latorig + 0.)**2./(2*9)))   &
              &- (50.-Pac_SPCZextra*.28 - Hawaiiextra*.6)*exp(-(lon*180./pi - 250)**2./(2*81))*(exp(-(latorig + 1.)**2./(2*9)))   &
              &- (50.-Pac_SPCZextra*.28 - Hawaiiextra*.6)*exp(-(lon*180./pi - 230)**2./(2*81))*(exp(-(latorig + 2.)**2./(2*9)))   &
              &- (39.-Hawaiiextra*.6)*exp(-(lon*180./pi - 210)**2./(2*81))*(exp(-(latorig + 2.)**2./(2*9)))   &
              & - (36.-Hawaiiextra*.676)*exp(-(lon*180./pi - 190)**2./(2*81))*(exp(-(latorig + 0.)**2./(2*9))) &
              & - (16.-Hawaiiextra*.676)*exp(-(lon*180./pi - 170)**2./(2*81))*(exp(-(latorig + 0.)**2./(2*9))) &
              &- 40.*exp(-(lon*180./pi - 287)**2./(2*4))*(exp(-(latorig + 25.)**2./(2*81)))   &
              &- 15.*exp(-(lon*180./pi - 282)**2./(2*25))*(exp(-(latorig + 15.)**2./(2*81)))   &
              &- (25.+Pac_ITCZextra + Pac_SPCZextra - Hawaiiextra*.5)*exp(-(lon*180./pi - 240)**2./(2*1600)) &
                & *(exp(-(latorig + 21.)**2./(2*121)))   &
              &- (38.0 - Hawaiiextra*.7)*exp(-(lon*180./pi - 195)**2./(2*169))*(exp(-(latorig - 16.)**2./(2*49)))   &
              &- (51.4 - Hawaiiextra*.7)*exp(-(lon*180./pi - 225)**2./(2*169))*(exp(-(latorig - 16.)**2./(2*49)))   &
              &+ (28.2 - Hawaiiextra*25/30*.8623)*exp(-(lon*180./pi - 220)**2./(2*1600))*exp(-(latorig + 57)**2./(2*225)) &
              &+ (14.+Pac_SPCZextra*.9)*exp(-(lon*180./pi - 165)**2./(2*400))*exp(-(latorig + 20)**2./(2*25))&
              &+ (16.+Pac_SPCZextra*.9)*exp(-(lon*180./pi - 195)**2./(2*400))*exp(-(latorig + 33)**2./(2*49)) &
              &+ (50.+Pac_SPCZextra*.9)*exp(-(lon*180./pi - 155)**2./(2*9))*exp(-(latorig + 30)**2./(2*49)) &
              &+ (40.+Pac_SPCZextra*.9)*exp(-(lon*180./pi - 180)**2./(2*25))*exp(-(latorig + 40)**2./(2*25))&
              &+ (41.-Hawaiiextra*25/30)*exp(-(lon*180./pi - 240)**2./(2*900))*exp(-(latorig + 62)**2./(2*64)) &
                &+ (60.-Hawaiiextra)*exp(-(lon*180./pi - 180)**2./(2*169))*(exp(-(latorig - 6.97)**2./(2*4))) &
              &+ (47.-Hawaiiextra)*exp(-(lon*180./pi - 210)**2./(2*169))*(exp(-(latorig - 6.97)**2./(2*4))) &
              &+ (45.-Hawaiiextra)*exp(-(lon*180./pi - 240)**2./(2*169))*(exp(-(latorig - 6.97)**2./(2*4))) &
              &+ (19.5 + Pac_SPCZextra*.4 - Hawaiiextra*.5158)*exp(-(lon*180./pi - 145)**2./(2*196)) &
                & *(exp(-(latorig - 3.)**2./(2*16)))  &
              &+ (40.+Pac_SPCZextra*.35 - Hawaiiextra*.5158)*exp(-(lon*180./pi - 150)**2./(2*169))*(exp(-(latorig - 7.)**2./(2*9)))
        end if
        !add south pacific and SPCZ at expense of further cold tongue  modified by cig, may 28 2019, also includes Australia bit
        !from Kuroshio above
        if (warmpool_localization_choice == 3 .and. warmpool_amp .gt. 0.001 .and. latorig .le. 30. .and. latorig .ge. -16. &
          & .and. lon .ge. 50.*pi/180. .and. lon .lt. 220.*pi/180.) then
          flux(i, j) = flux(i, j) &
                &+ Hawaiiextra*.9*exp(-(lon*180./pi - 145)**2./(2*400))*(exp(-(latorig - 16.)**2./(2*25)))
        end if

        !add south pacific and south indian at expense of Australia  modified by cig, may 13 2019
        if (warmpool_localization_choice == 3 .and. warmpool_amp .gt. 0.001 .and. latorig .le. 10. .and. latorig .ge. -36. &
          & .and. lon .ge. 50.*pi/180. .and. lon .lt. 220.*pi/180.) then
          flux(i, j) = flux(i, j) &
              &- (qflux_amp + warmpool_amp)*1.02*exp(-(lon*180./pi - 135)**2./(2*225))*(exp(-(latorig + 20.)**2./(2*36)))   &
              &- (qflux_amp + warmpool_amp)*1.02*10/44*exp(-(lon*180./pi - 147)**2./(2*64))*(exp(-(latorig + 27.)**2./(2*49)))   &
              &+ (qflux_amp + warmpool_amp)*1.02*16.6/44*exp(-(lon*180./pi - 120)**2./(2*900))*exp(-(latorig + 20)**2./(2*36)) &
              &+ (qflux_amp + warmpool_amp)*1.02*(27.89 - Hawaiiextra*1.412)/44*exp(-(lon*180./pi - 100)**2./(2*100)) &
                & *exp(-(latorig + 10)**2./(2*16)) &
              &+ (qflux_amp + warmpool_amp)*1.02*4.9/44*exp(-(lon*180./pi - 135)**2./(2*225))*exp(-(latorig + 0)**2./(2*16))  &
              &+ (Pac_SPCZextra*.7137)*exp(-(lon*180./pi - 130)**2./(2*196))*(exp(-(latorig + 8.)**2./(2*81)))

        end if

        !add  Gulf streams
        if (warmpool_localization_choice == 3 .and. gulf_amp .gt. 0.001 .and. latorig .le. 52. .and. latorig .ge. 10.) then
          !modified by cig, may 13 2019, Atlantic sector is exp
          if (lon .ge. (warmpool_phase + 135.)*pi/180. .and. lon .lt. (warmpool_phase + 195.)*pi/180.) then
            flux(i, j) = flux(i, j) &
               &+ gulf_amp*exp(-(lon*180./pi - 2.*latorig - 220.)**2./(2*9.))*exp(-(lon*180./pi + latorig - 335.)**2./(2*625.))
          end if
          !modified by cig, may 13 2019
          if (lon .ge. (warmpool_phase + 158.)*pi/180. .and. lon .lt. (warmpool_phase + 218.)*pi/180.) then
            flux(i, j) = flux(i, j) &
               &- gulf_amp*63.9/70.*exp(-(lon*180./pi - .5*latorig - 325.)**2./(2*9.))*exp(-(latorig - 25.)**2./(2*49.))
          end if

        end if

        !add more Gulf at expense of tropical Atlantic  modified by cig, may 14 2019
        if (warmpool_localization_choice == 3 .and. trop_atlantic_amp .gt. 0.001 .and. latorig .le. 77. .and. latorig .ge. -35. &
          & .and. lon .ge. 275.*pi/180. .and. gulf_amp .gt. 0.001) then
          flux(i, j) = flux(i, j) &
              &- trop_atlantic_amp*(exp(-(lon*180./pi - 342)**2./(2*81)) + exp(-(lon*180./pi - 360.)**2./(2*64))) &
                & *(exp(-(latorig + 5.)**2./(2*25)))   &
              &- 12.6*(exp(-(lon*180./pi - 345)**2./(2*256)))*(exp(-(latorig + 16.)**2./(2*64)))   &
              &+ trop_atlantic_amp*30.65/28.*exp(-(lon*180./pi - 2*latorig - 220)**2./(2*100)) &
                & *exp(-(latorig + lon*180./pi - 375)**2./(2*900))
        end if
        if (warmpool_localization_choice == 3 .and. trop_atlantic_amp .gt. 0.001 .and. gulf_amp .gt. 0.001) then
          !second part of Gulf at expense of South America, replaced greenland perturbation
          if (lon .ge. (318.)*pi/180. .or. lon .lt. (18.)*pi/180.) then
            if (abs(latgreen) .le. 1.0) then
              flux(i, j) = flux(i, j) &
                 & + trop_atlantic_amp*36./28*(1.-latgreen**4.)*cos((gulf_k - 1)*(lon - pigulfphase - 38.*pi/180.))
            elseif (latorig .le. 15. .and. latorig .ge. -35.) then
              flux(i, j) = flux(i, j) &
                &- trop_atlantic_amp*exp(-(lon*180./pi - 0.)**2./(2*64))*(exp(-(latorig + 5.)**2./(2*25)))  &
                 &- 12.6*exp(-(lon*180./pi + 15.)**2./(2*256))*(exp(-(latorig + 16.)**2./(2*64)))
            end if
          end if
        end if

        !add Caribean and South America at expense of tropical North/South Atlantic  modified by cig, may 14 2019
        if (warmpool_localization_choice == 3 .and. trop_atlantic_amp .gt. 0.001 .and. lon .ge. 250.*pi/180. &
          & .and. lon .lt. 344.*pi/180. .and. latorig .le. 40. .and. latorig .ge. -35. .and. gulf_amp .gt. 0.001) then
          flux(i, j) = flux(i, j) &
           &- qflux_amp*.92*exp(-(lon*180./pi - 290)**2./(2*400))*(exp(-(latorig + 20)**2./(2*49)))   &
           &- 16.8*exp(-(lon*180./pi - 325)**2./(2*484))*(exp(-(latorig - 19.5)**2./(2*64)))   &
           &+ qflux_amp*1.2*exp(-(lon*180./pi - 270)**2./(2*49))*exp(-(latorig - 22)**2./(2*25)) &
           &+ qflux_amp*1.58*exp(-(lon*180./pi - 283)**2./(2*25))*(exp(-(latorig + 0)**2./(2*36)))   &
           &+ qflux_amp*1.06415*exp(-(lon*180./pi - 304)**2./(2*36))*(exp(-(latorig + 2)**2./(2*49)))  &
                &+ qflux_amp*.85*exp(-(lon*180./pi - 284)**2./(2*25))*(exp(-(latorig + 10)**2./(2*36)))   &
           &+ qflux_amp*.63*exp(-(lon*180./pi - 317)**2./(2*25))*(exp(-(latorig + 6)**2./(2*16)))  &
           &+ (42.54)*exp(-(lon*180./pi - 325)**2./(2*121))*(exp(-(latorig - 4.2)**2./(2*4)))
        end if

        !add Agulhas current add expense of Africa  modified by cig, may 13 2019
        if (warmpool_localization_choice == 3 .and. africaamp .gt. 0.001 .and. latorig .le. 35. .and. latorig .ge. -60. &
          & .and. lon .ge. 2.*pi/180. .and. lon .lt. 100.*pi/180.) then
          flux(i, j) = flux(i, j) &
              &- 30.*exp(-(lon*180./pi - 28)**2./(2*100))*(exp(-(latorig - 18.)**2./(2*50)) + exp(-(latorig + 18)**2./(2*60))) &
              &- (38.5 + Africaextra*.7709)*exp(-(lon*180./pi - 11)**2./(2*4))*exp(-(latorig + 15)**2./(2*100)) &
              &+ (83.+Africaextra)*exp(-(lon*180./pi - 50)**2./(2*625))*exp(-(latorig + 40)**2./(2*16)) &
              &- (64.22 + Africaextra*1.3)*exp(-(lon*180./pi - 50)**2./(2*400))*exp(-(latorig + 48)**2./(2*16)) &
              &+ (38.0 + Africaextra/3)*exp(-(lon*180./pi - 2./3.*latorig - 57.)**2./(2*16)) &
                & *exp(-(lon*180./pi + latorig - 10.)**2./(2*225.)) &
              &+ 20.*exp(-(lon*180./pi - 14)**2./(2*30))*(exp(-(latorig - 0.)**2./(2*50))) &
              &+ 11.*exp(-(lon*180./pi - 36)**2./(2*30))*(exp(-(latorig - 0.)**2./(2*50)))

        end if

        !add Agulhas current add expense of Africa  modified by cig, feb 19 2020
        if (warmpool_localization_choice == 3 .and. latorig .le. -30. .and. latorig .ge. -61. .and. abs(Sampeextra) .gt. 0.001) then
          flux(i, j) = flux(i, j) &
              &+ (Sampeextra*.8822)*exp(-(latorig + 40)**2./(2*16)) &
              &- (Sampeextra)*exp(-(latorig + 48)**2./(2*16))

        end if
        !add zonally symmetric heating modified by cig, June 2 2020
        if (warmpool_localization_choice == 3 .and. latorig .le. -45. .and. latorig .ge. -65. &
          & .and. abs(Pac_ITCZextra) .gt. 0.001) then
          flux(i, j) = flux(i, j) &
              &+ (Pac_ITCZextra*.74537)*exp(-(latorig + 55)**2./(2*49))

        end if

        !add dipole in south atlantic  modified by cig, august 4 2019
        if (warmpool_localization_choice == 3 .and. latorig .ge. -61. .and. latorig .le. -30. .and. lon .ge. 290.*pi/180.) then
          flux(i, j) = flux(i, j) &
              &+ (37.4 + Africaextra*0*.9349/2)*exp(-(lon*180./pi - 323)**2./(2*121))*exp(-(latorig + 36)**2./(2*16)) &
              &- (40.+0*Africaextra/2)*exp(-(lon*180./pi - 311)**2./(2*121))*exp(-(latorig + 45)**2./(2*16))

        end if

        !add Norwegian Sea at expense of Africa modified by cig, may 28 2019
        if (warmpool_localization_choice == 3 .and. trop_atlantic_amp .gt. 0.001 .and. gulf_amp .gt. 0.001) then
          if (lon .ge. 310.*pi/180. .and. latorig .ge. 10. .and. latorig .le. 35.) then
            flux(i, j) = flux(i, j) &
               &- qflux_amp*14.5/26.*exp(-(lon*180./pi - 357.)**2./(2*400))*(exp(-(latorig - 20)**2./(2*49)))
          end if
          if (lon .le. 30.*pi/180. .and. latorig .ge. 10. .and. latorig .le. 35.) then
            flux(i, j) = flux(i, j) &
             &- qflux_amp*14.5/26.*exp(-(lon*180./pi + 3)**2./(2*400))*(exp(-(latorig - 20)**2./(2*49)))
          end if
          if ((lon .le. 75.*pi/180. .or. lon .ge. (345.)*pi/180.) .and. latorig .ge. 71. .and. latorig .le. 83.) then
            flux(i, j) = flux(i, j) &
              & + qflux_amp*68./26.*(1.-((latorig - 76.)/6.5)**4.)*cos(2*(lon - 30.*pi/180.))
          end if
        end if

        !wave1 over Arctic and Hudson Bay modified by cig, may 30 2019
        if (warmpool_localization_choice == 3 .and. trop_atlantic_amp .gt. 0.001) then
          if (latorig .ge. 69. .and. latorig .le. 83.) then
            flux(i, j) = flux(i, j) &
               & + 25.*(1.-((latorig - 76.)/7.)**4.)*cos(1.*(lon - 10.*pi/180.))
          end if
          if (latorig .ge. 60. .and. latorig .le. 76.) then
            if ((lon .le. 17.*pi/180. .or. lon .ge. (347.)*pi/180.)) then
              flux(i, j) = flux(i, j) &
                   & + 68.2*(1.-((latorig - 68.)/8.)**4.)*cos(6*(lon - 2.*pi/180.))
            end if
          end if
          if (lon .ge. 260.*pi/180. .and. lon .le. 310.*pi/180. .and. latorig .ge. 55. .and. latorig .le. 85.) then
            flux(i, j) = flux(i, j) &
              &- 38*exp(-(lon*180./pi - 2.*latorig - 152.)**2./(2*100.))*exp(-(lon*180./pi + latorig - 342.)**2./(2*400.)) &
              &- 100*exp(-(lon*180./pi - 275.)**2./(2*25.))*exp(-(latorig - 58.)**2./(2*16.))
          end if
          if (lon .ge. 275.*pi/180. .and. lon .le. 335.*pi/180. .and. latorig .ge. 10. .and. latorig .le. 52.) then
            flux(i, j) = flux(i, j) &
                &+ 10.8*exp(-(lon*180./pi - 2.*latorig - 220.)**2./(2*100.))*exp(-(lon*180./pi + latorig - 335.)**2./(2*625.))
          end if
        end if

      end do
    end do
  end subroutine warmpool

!########################################################

end module qflux_mod

