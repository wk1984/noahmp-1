module noahmp_func
  use noahmp_const
  use noahmp_global
  implicit none

  public  :: noahmp_set_options
  public  :: noahmp_sflx
  public  :: frh2o

  private :: atm
  private :: phenology
  private :: energy
  private ::       thermoprop
  private ::               csnow
  private ::               tdfcnd
  private ::       radiation
  private ::               albedo
  private ::                         snowage
  private ::                         snowalb_bats
  private ::                         snowalb_class
  private ::                         groundalb
  private ::                         twostream
  private ::               surrad
  private ::       vege_flux
  private ::               sfcdif1
  private ::               sfcdif2
  private ::               stomata
  private ::               canres
  private ::               esat
  private ::               ragrb
  private ::       bare_flux
  private ::       tsnosoi
  private ::               hrt
  private ::               hstep
  private ::                         rosr12
  private ::       phasechange

  private :: water
  private ::       canwater
  private ::       snowwater
  private ::               snowfall
  private ::               combine
  private ::               divide
  private ::                         combo
  private ::               compact
  private ::               snowh2o
  private ::       soilh2o
  private ::               zwteq
  private ::               infil
  private ::               srt
  private ::                         wdfcnd1
  private ::                         wdfcnd2
  !  private ::                         INFIL
  private ::               sstep
  private ::       groundwater

  private :: carbon
  private ::       co2flux
  !  private ::       BVOCFLUX
  !  private ::       CH4FLUX

  private :: error

contains

  subroutine noahmp_sflx( &
       ILOC    , JLOC    , LAT     , YEARLEN , JULIAN  , COSZ    , & ! IN : Time/Space-related
       DT      , DX      , DZ8W    , NSOIL   , ZSOIL   , NSNOW   , & ! IN : Model configuration
       SHDFAC  , SHDMAX  , slptyp, sltyp, lutyp, ICE     , IST     , & ! IN : Vegetation/Soil characteristics
       ISC,                                         & ! IN : Vegetation/Soil characteristics
       IZ0TLND ,                                                   & ! IN : User options
       SFCTMP  , SFCPRS  , PSFC    , UU      , VV      , Q2      , & ! IN : Forcing
       QC      , SOLDN   , LWDN    , PRCP    , TBOT    , CO2AIR  , & ! IN : Forcing
       O2AIR   , FOLN    , FICEOLD , PBLH    , ZLVL    ,           & ! IN : Forcing
       ALBOLD  , SNEQVO  ,                                         & ! IN/OUT :
       STC     , soilwat    , SMC     , TAH     , EAH     , FWET    , & ! IN/OUT :
       CANLIQ  , CANICE  , TV      , TG      , QSFC    , QSNOW   , & ! IN/OUT :
       ISNOW   , ZSNSO   , snowh   , sneqv   , SNICE   , SNLIQ   , & ! IN/OUT :
       ZWT     , WA      , WT      , WSLAKE  , LFMASS  , RTMASS  , & ! IN/OUT :
       STMASS  , WOOD    , STBLCP  , FASTCP  , LAI     , SAI     , & ! IN/OUT :
       CM      , CH      , TAUSS   ,                               & ! IN/OUT :
       FSA     , FSR     , FIRA    , FSH     , SSOIL   , FCEV    , & ! OUT :
       FGEV    , FCTR    , ECAN    , ETRAN   , EDIR    , TRAD    , & ! OUT :
       TGB     , TGV     , T2MV    , T2MB    , Q2V     , Q2B     , & ! OUT :
       runsrf  , runsub  , APAR    , PSN     , SAV     , SAG     , & ! OUT :
       FSNO    , NEE     , GPP     , NPP     , fveg    , ALBEDO  , & ! OUT :
       QSNBOT  , PONDING , PONDING1, PONDING2, RSSUN   , RSSHA   , & ! OUT :
       BGAP    , WGAP    , CHV     , CHB     , EMISSI  ,           & ! OUT :
       SHG     , SHC     , SHB     , EVG     , EVB     , GHV     , & ! OUT :
       GHB     , IRG     , IRC     , IRB     , TR      , EVC     , & ! OUT :
       CHLEAF  , CHUC    , CHV2    , CHB2    , FPICE)
    ! Initial code: Guo-Yue Niu, Oct. 2007
    use noahmp_gen_param, only: KK_CSOIL
    use noahmp_gen_param, only: KK_ZBOT
    use noahmp_soil_param, only: LK_SMCMAX
    use noahmp_veg_param, only: ISURBAN
    use noahmp_veg_param, only: ISBARREN
    use noahmp_veg_param, only: LK_NROOT
    implicit none
    ! input
    integer                        , intent(in)    :: ICE    !ice (ice = 1)
    integer                        , intent(in)    :: IST    !surface type 1->soil; 2->lake
    integer                        , intent(in)    :: slptyp !slope type
    integer                        , intent(in)    :: sltyp !soil type
    integer                        , intent(in)    :: lutyp !vegetation type
    integer                        , intent(in)    :: ISC    !soil color type (1-lighest; 8-darkest)
    integer                        , intent(in)    :: NSNOW  !maximum no. of snow layers
    integer                        , intent(in)    :: NSOIL  !no. of soil layers
    integer                        , intent(in)    :: ILOC   !grid index
    integer                        , intent(in)    :: JLOC   !grid index
    real                           , intent(in)    :: DT     !time step [sec]
    real, dimension(       1:NSOIL), intent(in)    :: ZSOIL  !layer-bottom depth from soil surf (m)
    real                           , intent(in)    :: Q2     !mixing ratio (kg/kg) lowest model layer
    real                           , intent(in)    :: SFCTMP !surface air temperature [K]
    real                           , intent(in)    :: UU     !wind speed in eastward dir (m/s)
    real                           , intent(in)    :: VV     !wind speed in northward dir (m/s)
    real                           , intent(in)    :: SOLDN  !downward shortwave radiation (w/m2)
    real                           , intent(in)    :: PRCP   !precipitation rate (kg m-2 s-1)
    real                           , intent(in)    :: LWDN   !downward longwave radiation (w/m2)
    real                           , intent(in)    :: SFCPRS !pressure (pa)
    real                           , intent(inout) :: ZLVL   !reference height (m)
    real                           , intent(in)    :: COSZ   !cosine solar zenith angle [0-1]
    real                           , intent(in)    :: TBOT   !bottom condition for soil temp. [K]
    real                           , intent(in)    :: FOLN   !foliage nitrogen (%) [1-saturated]
    real                           , intent(in)    :: SHDFAC !green vegetation fraction [0.0-1.0]
    integer                        , intent(in)    :: YEARLEN!Number of days in the particular year.
    real                           , intent(in)    :: JULIAN !Julian day of year (floating point)
    real                           , intent(in)    :: LAT    !latitude (radians)
    real, dimension(-NSNOW+1:    0), intent(in)    :: FICEOLD!ice fraction at last timestep

    !jref:start; in
    integer                        , intent(in)    :: IZ0TLND
    real                           , intent(in)    :: QC     !cloud water mixing ratio
    real                           , intent(in)    :: PBLH   !planetary boundary layer height
    real                           , intent(inout) :: QSFC   !mixing ratio at lowest model layer
    real                           , intent(in)    :: PSFC   !pressure at lowest model layer
    real                           , intent(in)    :: DZ8W   !thickness of lowest layer
    real                           , intent(in)    :: DX
    real                           , intent(in)    :: SHDMAX  !yearly max vegetation fraction
    !jref:end

    ! input/output : need arbitary intial values
    real                           , intent(inout) :: QSNOW  !snowfall [mm/s]
    real                           , intent(inout) :: FWET   !wetted or snowed fraction of canopy (-)
    real                           , intent(inout) :: SNEQVO !snow mass at last time step (mm)
    real                           , intent(inout) :: EAH    !canopy air vapor pressure (pa)
    real                           , intent(inout) :: TAH    !canopy air tmeperature (k)
    real                           , intent(inout) :: ALBOLD !snow albedo at last time step (CLASS type)
    real                           , intent(inout) :: CM     !momentum drag coefficient
    real                           , intent(inout) :: CH     !sensible heat exchange coefficient
    real                           , intent(inout) :: TAUSS  !non-dimensional snow age

    ! prognostic variables
    integer                        , intent(inout) :: ISNOW  !actual no. of snow layers [-]
    real                           , intent(inout) :: CANLIQ !intercepted liquid water (mm)
    real                           , intent(inout) :: CANICE !intercepted ice mass (mm)
    real                           , intent(inout) :: sneqv  !snow water eqv. [mm]
    real, dimension(       1:NSOIL), intent(inout) :: SMC    !soil moisture (ice + liq.) [m3/m3]
    real, dimension(-NSNOW+1:NSOIL), intent(inout) :: ZSNSO  !layer-bottom depth from snow surf [m]
    real                           , intent(inout) :: snowh  !snow height [m]
    real, dimension(-NSNOW+1:    0), intent(inout) :: SNICE  !snow layer ice [mm]
    real, dimension(-NSNOW+1:    0), intent(inout) :: SNLIQ  !snow layer liquid water [mm]
    real                           , intent(inout) :: TV     !vegetation temperature (k)
    real                           , intent(inout) :: TG     !ground temperature (k)
    real, dimension(-NSNOW+1:NSOIL), intent(inout) :: STC    !snow/soil temperature [k]
    real, dimension(       1:NSOIL), intent(inout) :: soilwat   !liquid soil moisture [m3/m3]
    real                           , intent(inout) :: ZWT    !depth to water table [m]
    real                           , intent(inout) :: WA     !water storage in aquifer [mm]
    real                           , intent(inout) :: WT     !water in aquifer&saturated soil [mm]
    real                           , intent(inout) :: WSLAKE !lake water storage (can be neg.) (mm)

    ! output
    real                           , intent(out)   :: FSA    !total absorbed solar radiation (w/m2)
    real                           , intent(out)   :: FSR    !total reflected solar radiation (w/m2)
    real                           , intent(out)   :: FIRA   !total net LW rad (w/m2)  [+ to atm]
    real                           , intent(out)   :: FSH    !total sensible heat (w/m2) [+ to atm]
    real                           , intent(out)   :: FCEV   !canopy evap heat (w/m2) [+ to atm]
    real                           , intent(out)   :: FGEV   !ground evap heat (w/m2) [+ to atm]
    real                           , intent(out)   :: FCTR   !transpiration heat (w/m2) [+ to atm]
    real                           , intent(out)   :: SSOIL  !ground heat flux (w/m2)   [+ to soil]
    real                           , intent(out)   :: TRAD   !surface radiative temperature (k)
    real                                           :: TS     !surface temperature (k)
    real                           , intent(out)   :: ECAN   !evaporation of intercepted water (mm/s)
    real                           , intent(out)   :: ETRAN  !transpiration rate (mm/s)
    real                           , intent(out)   :: EDIR   !soil surface evaporation rate (mm/s]
    real                           , intent(out)   :: runsrf !surface runoff [mm/s]
    real                           , intent(out)   :: runsub !baseflow (saturation excess) [mm/s]
    real                           , intent(out)   :: PSN    !total photosynthesis (umol co2/m2/s) [+]
    real                           , intent(out)   :: APAR   !photosyn active energy by canopy (w/m2)
    real                           , intent(out)   :: SAV    !solar rad absorbed by veg. (w/m2)
    real                           , intent(out)   :: SAG    !solar rad absorbed by ground (w/m2)
    real                           , intent(out)   :: FSNO   !snow cover fraction on the ground (-)
    real                           , intent(out)   :: fveg   !green vegetation fraction [0.0-1.0]
    real                           , intent(out)   :: ALBEDO !surface albedo [-]
    real                                           :: ERRWAT !water error [kg m{-2}]
    real                           , intent(out)   :: QSNBOT !snowmelt out bottom of pack [mm/s]
    real                           , intent(out)   :: PONDING!surface ponding [mm]
    real                           , intent(out)   :: PONDING1!surface ponding [mm]
    real                           , intent(out)   :: PONDING2!surface ponding [mm]

    !jref:start; output
    real                           , intent(out)     :: T2MV   !2-m air temperature over vegetated part [k]
    real                           , intent(out)     :: T2MB   !2-m air temperature over bare ground part [k]
    real, intent(out) :: RSSUN        !sunlit leaf stomatal resistance (s/m)
    real, intent(out) :: RSSHA        !shaded leaf stomatal resistance (s/m)
    real, intent(out) :: BGAP
    real, intent(out) :: WGAP
    real, intent(out) :: TGV
    real, intent(out) :: TGB
    real              :: Q1
    real, intent(out) :: EMISSI
    !jref:end

    ! local
    integer                                        :: IZ     !do-loop index
    integer, dimension(-NSNOW+1:NSOIL)             :: IMELT  !phase change index [1-melt; 2-freeze]
    real                                           :: CMC    !intercepted water (CANICE+CANLIQ) (mm)
    real                                           :: TAUX   !wind stress: e-w (n/m2)
    real                                           :: TAUY   !wind stress: n-s (n/m2)
    real                                           :: RHOAIR !density air (kg/m3)
    !  REAL, DIMENSION(       1:    5)                :: VOCFLX !voc fluxes [ug C m-2 h-1]
    real, dimension(-NSNOW+1:NSOIL)                :: DZSNSO !snow/soil layer thickness [m]
    real                                           :: THAIR  !potential temperature (k)
    real                                           :: QAIR   !specific humidity (kg/kg) (q2/(1+q2))
    real                                           :: EAIR   !vapor pressure air (pa)
    real, dimension(       1:    2)                :: SOLAD  !incoming direct solar rad (w/m2)
    real, dimension(       1:    2)                :: SOLAI  !incoming diffuse solar rad (w/m2)
    real                                           :: QPRECC !convective precipitation (mm/s)
    real                                           :: QPRECL !large-scale precipitation (mm/s)
    real                                           :: IGS    !growing season index (0=off, 1=on)
    real                                           :: elai   !leaf area index, after burying by snow
    real                                           :: esai   !stem area index, after burying by snow
    real                                           :: BEVAP  !soil water evaporation factor (0 - 1)
    real, dimension(       1:NSOIL)                :: BTRANI !Soil water transpiration factor (0 - 1)
    real                                           :: BTRAN  !soil water transpiration factor (0 - 1)
    real                                           :: HTOP   !top of canopy layer (m)
    real                                           :: QIN    !groundwater recharge [mm/s]
    real                                           :: QDIS   !groundwater discharge [mm/s]
    real, dimension(       1:NSOIL)                :: soilice   !soil ice content (m3/m3)
    real, dimension(-NSNOW+1:    0)                :: SNICEV !partial volume ice of snow [m3/m3]
    real, dimension(-NSNOW+1:    0)                :: SNLIQV !partial volume liq of snow [m3/m3]
    real, dimension(-NSNOW+1:    0)                :: epore  !effective porosity [m3/m3]
    real                                           :: TOTSC  !total soil carbon (g/m2)
    real                                           :: TOTLB  !total living carbon (g/m2)
    real                                           :: T2M    !2-meter air temperature (k)
    real                                           :: QDEW   !ground surface dew rate [mm/s]
    real                                           :: QVAP   !ground surface evap. rate [mm/s]
    real                                           :: LATHEA !latent heat [j/kg]
    real                                           :: SWDOWN !downward solar [w/m2]
    real                                           :: QMELT  !snowmelt [mm/s]
    real                                           :: BEG_WB !water storage at begin of a step [mm]
    real,intent(out)                                              :: IRC    !canopy net LW rad. [w/m2] [+ to atm]
    real,intent(out)                                              :: IRG    !ground net LW rad. [w/m2] [+ to atm]
    real,intent(out)                                              :: SHC    !canopy sen. heat [w/m2]   [+ to atm]
    real,intent(out)                                              :: SHG    !ground sen. heat [w/m2]   [+ to atm]
    real,intent(out)                                              :: EVG    !ground evap. heat [w/m2]  [+ to atm]
    real,intent(out)                                              :: GHV    !ground heat flux [w/m2]  [+ to soil]
    real,intent(out)                                              :: IRB    !net longwave rad. [w/m2] [+ to atm]
    real,intent(out)                                              :: SHB    !sensible heat [w/m2]     [+ to atm]
    real,intent(out)                                              :: EVB    !evaporation heat [w/m2]  [+ to atm]
    real,intent(out)                                              :: GHB    !ground heat flux [w/m2] [+ to soil]
    real,intent(out)                                              :: EVC    !canopy evap. heat [w/m2]  [+ to atm]
    real,intent(out)                                              :: TR     !transpiration heat [w/m2] [+ to atm]
    real, intent(out)   :: FPICE   !snow fraction in precipitation

    !jref:start
    real                                           :: FSRV
    real                                           :: FSRG
    real,intent(out)                               :: Q2V
    real,intent(out)                               :: Q2B
    real :: Q2E
    real :: QFX
    real,intent(out)                               :: CHV    !sensible heat exchange coefficient over vegetated fraction
    real,intent(out)                               :: CHB    !sensible heat exchange coefficient over bare-ground
    real,intent(out)                               :: CHLEAF !leaf exchange coefficient
    real,intent(out)                               :: CHUC   !under canopy exchange coefficient
    real,intent(out)                               :: CHV2    !sensible heat exchange coefficient over vegetated fraction
    real,intent(out)                               :: CHB2    !sensible heat exchange coefficient over bare-ground
    !jref:end

    ! carbon
    ! inputs
    real                           , intent(in)    :: CO2AIR !atmospheric co2 concentration (pa)
    real                           , intent(in)    :: O2AIR  !atmospheric o2 concentration (pa)

    ! inputs and outputs : prognostic variables
    real                        , intent(inout)    :: LFMASS !leaf mass [g/m2]
    real                        , intent(inout)    :: RTMASS !mass of fine roots [g/m2]
    real                        , intent(inout)    :: STMASS !stem mass [g/m2]
    real                        , intent(inout)    :: WOOD   !mass of wood (incl. woody roots) [g/m2]
    real                        , intent(inout)    :: STBLCP !stable carbon in deep soil [g/m2]
    real                        , intent(inout)    :: FASTCP !short-lived carbon, shallow soil [g/m2]
    real                        , intent(inout)    :: LAI    !leaf area index [-]
    real                        , intent(inout)    :: SAI    !stem area index [-]

    ! outputs
    real                          , intent(out)    :: NEE    !net ecosys exchange (g/m2/s CO2)
    real                          , intent(out)    :: GPP    !net instantaneous assimilation [g/m2/s C]
    real                          , intent(out)    :: NPP    !net primary productivity [g/m2/s C]
    real                                           :: AUTORS !net ecosystem respiration (g/m2/s C)
    real                                           :: HETERS !organic respiration (g/m2/s C)
    real                                           :: TROOT  !root-zone averaged temperature (k)
    real                                 :: LATHEAV !latent heat vap./sublimation (j/kg)
    real                                 :: LATHEAG !latent heat vap./sublimation (j/kg)
    logical                             :: FROZEN_GROUND ! used to define latent heat pathway
    logical                             :: FROZEN_CANOPY ! used to define latent heat pathway

    ! INTENT (out) variables need to be assigned a value.  These normally get assigned values
    ! only if opt_veg == 2.
    nee = 0.0
    npp = 0.0
    gpp = 0.0

    ! re-process atmospheric forcing

    call ATM(SFCPRS, SFCTMP, Q2, PRCP, SOLDN, &
         & COSZ, THAIR, QAIR, EAIR, RHOAIR, &
         & QPRECC, QPRECL, SOLAD, SOLAI, SWDOWN)

    ! snow/soil layer thickness (m)

    do IZ = ISNOW+1, NSOIL
       if (IZ == ISNOW + 1) then
          DZSNSO(IZ) = -ZSNSO(IZ)
       else
          DZSNSO(IZ) = ZSNSO(IZ-1) - ZSNSO(IZ)
       end if
    end do

    ! root-zone temperature

    TROOT  = 0.0
    do IZ = 1, LK_NROOT(lutyp)
       TROOT = TROOT + STC(IZ) * DZSNSO(IZ) / (-ZSOIL(LK_NROOT(lutyp)))
    end do

    ! total water storage for water balance check

    if (IST == 1) then
       BEG_WB = CANLIQ + CANICE + sneqv + WA
       do IZ = 1, NSOIL
          BEG_WB = BEG_WB + SMC(IZ) * DZSNSO(IZ) * 1000.0
       end do
    end if

    ! vegetation phenology

    call phenology(lutyp, snowh, TV, LAT, &
         & YEARLEN, JULIAN, LAI, SAI, TROOT, &
         & HTOP, elai, esai, IGS)

    !input GVF should be consistent with LAI
    !     if (opt_veg == 1) THEN
    !        fveg = SHDFAC
    !        if (fveg <= 0.05) fveg = 0.05
    !     ELSE IF (opt_veg == 2 .or. opt_veg == 3) THEN
    !        fveg = 1.-EXP(-0.52*(LAI+SAI))
    !        if (fveg <= 0.05) fveg = 0.05
    !     ELSE IF (opt_veg == 4) THEN
    !        fveg = SHDMAX
    !        if (fveg <= 0.05) fveg = 0.05
    !     ELSE
    !        WRITE(*,*) "-------- FATAL CALLED IN SFLX -----------"
    !        CALL wrf_error_fatal("Namelist parameter opt_veg unknown")
    !     ENDIF
    if (opt_veg == 1) then
       fveg = SHDFAC
       if (fveg <= 0.01) fveg = 0.01
    else if (opt_veg == 2 .or. opt_veg == 3) then
       fveg = 1.0 - exp(-0.52 * (LAI + SAI))
       if (fveg <= 0.01) fveg = 0.01
    else if (opt_veg == 4 .or. opt_veg == 5) then
       fveg = SHDMAX
       if (fveg <= 0.01) fveg = 0.01
    else
       write(*,*) "-------- FATAL CALLED IN SFLX -----------"
       call wrf_error_fatal("Namelist parameter opt_veg unknown")
    end if
    if (lutyp == ISURBAN .or. lutyp == ISBARREN) fveg = 0.0
    if (elai + esai == 0.0) fveg = 0.0

    !     CALL PHENOLOGY (lutyp,IMONTH ,IDAY   ,snowh  ,TV     ,LAT   , & !in
    !                     LAI   ,SAI    ,TROOT  ,                & !in
    !                     HTOP  ,elai   ,esai   ,IGS    )   !out

    ! compute energy budget (momentum & energy fluxes and phase changes)

    call energy(sltyp, lutyp, ICE, IST, ISC, NSNOW  ,NSOIL  , & !in
         ISNOW  ,DT     ,RHOAIR ,SFCPRS ,QAIR   , & !in
         SFCTMP ,THAIR  ,LWDN   ,UU     ,VV     ,ZLVL   , & !in
         CO2AIR ,O2AIR  ,SOLAD  ,SOLAI  ,COSZ   ,IGS    , & !in
         EAIR   ,HTOP   ,TBOT   ,KK_ZBOT   ,ZSNSO  ,ZSOIL  , & !in
         elai   ,esai   ,KK_CSOIL  ,FWET   ,FOLN   ,         & !in
         fveg   ,                                         & !in
         QSNOW  ,DZSNSO ,LAT    ,CANLIQ ,CANICE ,iloc, jloc , & !in
         IMELT  ,SNICEV ,SNLIQV ,epore  ,T2M    ,FSNO   , & !out
         SAV    ,SAG    ,QMELT  ,FSA    ,FSR    ,TAUX   , & !out
         TAUY   ,FIRA   ,FSH    ,FCEV   ,FGEV   ,FCTR   , & !out
         TRAD   ,PSN    ,APAR   ,SSOIL  ,BTRANI ,BTRAN  , & !out
         PONDING,TS     ,LATHEAV , LATHEAG , frozen_canopy,frozen_ground,                         & !out
         TV     ,TG     ,STC    ,snowh  ,EAH    ,TAH    , & !inout
         SNEQVO ,sneqv  ,soilwat   ,SMC    ,SNICE  ,SNLIQ  , & !inout
         ALBOLD ,CM     ,CH     ,DX     ,DZ8W   ,Q2     , & !inout
         TAUSS  ,                                         & !inout
                                !jref:start
         QC     ,PBLH   ,QSFC   ,PSFC, IZ0TLND, & !in
         T2MV   ,T2MB  ,FSRV   , &
         FSRG   ,RSSUN   ,RSSHA ,BGAP   ,WGAP, TGV,TGB,&
         Q1     ,Q2V    ,Q2B    ,Q2E    ,CHV   ,CHB     , & !out
         EMISSI,&
         SHG,SHC,SHB,EVG,EVB,GHV,GHB,IRG,IRC,IRB,TR,EVC,CHLEAF,CHUC,CHV2,CHB2 )                                            !out
    !jref:end

    soilice(:) = max(0.0, SMC(:) - soilwat(:))
    SNEQVO  = sneqv

    QVAP = max(FGEV / LATHEAG, 0.0)       ! positive part of fgev; Barlage change to ground v3.6
    QDEW = abs(min(FGEV / LATHEAG, 0.0))  ! negative part of fgev
    EDIR = QVAP - QDEW

    ! compute water budgets (water storages, ET components, and runoff)

    call water(slptyp, sltyp, lutyp, dt, nsnow, nsoil, IMELT,UU     , & !in
         VV     ,FCEV   ,FCTR   ,QPRECC ,QPRECL ,elai   , & !in
         esai   ,SFCTMP ,QVAP   ,QDEW   ,ZSOIL  ,BTRANI , & !in
         FICEOLD,PONDING,TG     ,IST    ,fveg   ,iloc,jloc, & !in
         LATHEAV , LATHEAG , frozen_canopy,frozen_ground,                        & !in  MB
         ISNOW  ,CANLIQ ,CANICE ,TV     ,snowh  ,sneqv  , & !inout
         SNICE  ,SNLIQ  ,STC    ,ZSNSO  ,soilwat   ,SMC    , & !inout
         soilice   ,ZWT    ,WA     ,WT     ,DZSNSO ,WSLAKE , & !inout
         CMC    ,ECAN   ,ETRAN  ,FWET   ,runsrf ,runsub , & !out
         QIN    ,QDIS   ,QSNOW  ,PONDING1       ,PONDING2,&
         QSNBOT,FPICE)  !out

    !     write(*,'(a20,10F15.5)') 'SFLX:RUNOFF=',runsrf*DT,runsub*DT,EDIR*DT

    ! compute carbon budgets (carbon storages and co2 & bvoc fluxes)

    if (opt_veg == 2 .or. opt_veg == 5) then
       call carbon(NSNOW  ,NSOIL  ,lutyp, DT     ,ZSOIL  , & !in
            DZSNSO ,STC    ,SMC    ,TV     ,TG     ,PSN    , & !in
            FOLN   ,LK_SMCMAX(sltyp) ,BTRAN  ,APAR   ,fveg   ,IGS    , & !in
            TROOT  ,IST    ,LAT    ,iloc   ,jloc, & !in
            LFMASS ,RTMASS ,STMASS ,WOOD   ,STBLCP ,FASTCP , & !inout
            GPP    ,NPP    ,NEE    ,AUTORS ,HETERS ,TOTSC  , & !out
            TOTLB  ,LAI    ,SAI    )                   !out
    end if

    ! water and energy balance check

    call error(SWDOWN ,FSA    ,FSR    ,FIRA   ,FSH    ,FCEV   , & !in
         FGEV   ,FCTR   ,SSOIL  ,BEG_WB ,CANLIQ ,CANICE , & !in
         sneqv  ,WA     ,SMC    ,DZSNSO ,PRCP   ,ECAN   , & !in
         ETRAN  ,EDIR   ,runsrf ,runsub ,DT     ,NSOIL  , & !in
         NSNOW  ,IST    ,ERRWAT ,ILOC   , JLOC  ,fveg   , &
         SAV    ,SAG    ,FSRV   ,FSRG   ,ZWT  )   !in ( Except ERRWAT, which is out )

    ! urban - jref
    QFX = ETRAN + ECAN + EDIR
    if (lutyp == ISURBAN) then
       QSFC = (QFX / RHOAIR * CH) + QAIR
       Q2B = QSFC
    end if

    if (snowh <= 1.0E-6 .or. sneqv <= 1.0E-3) then
       snowh = 0.0
       sneqv = 0.0
    end if

    if (SWDOWN /= 0.0) then
       ALBEDO = FSR / SWDOWN
    else
       ALBEDO = -999.9
    end if

  end subroutine noahmp_sflx


  subroutine atm(SFCPRS, SFCTMP ,Q2     ,PRCP   ,SOLDN  ,COSZ   ,THAIR  , &
       QAIR   ,EAIR   ,RHOAIR ,QPRECC ,QPRECL ,SOLAD  ,SOLAI  , &
       SWDOWN )
    ! re-process atmospheric forcing
    implicit none
    ! inputs
    real                          , intent(in)  :: SFCPRS !pressure (pa)
    real                          , intent(in)  :: SFCTMP !surface air temperature [k]
    real                          , intent(in)  :: Q2     !mixing ratio (kg/kg)
    real                          , intent(in)  :: SOLDN  !downward shortwave radiation (w/m2)
    real                          , intent(in)  :: PRCP   !precipitation rate (kg m-2 s-1)
    real                          , intent(in)  :: COSZ   !cosine solar zenith angle [0-1]

    ! outputs
    real                          , intent(out) :: THAIR  !potential temperature (k)
    real                          , intent(out) :: QAIR   !specific humidity (kg/kg) (q2/(1+q2))
    real                          , intent(out) :: EAIR   !vapor pressure air (pa)
    real, dimension(       1:   2), intent(out) :: SOLAD  !incoming direct solar radiation (w/m2)
    real, dimension(       1:   2), intent(out) :: SOLAI  !incoming diffuse solar radiation (w/m2)
    real                          , intent(out) :: QPRECC !convective precipitation (mm/s)
    real                          , intent(out) :: QPRECL !large-scale precipitation (mm/s)
    real                          , intent(out) :: RHOAIR !density air (kg/m3)
    real                          , intent(out) :: SWDOWN !downward solar filtered by sun angle [w/m2]

    !locals

    real                                        :: PAIR   !atm bottom level pressure (pa)

    !jref: seems like PAIR should be P1000mb??
    PAIR   = SFCPRS                   ! atm bottom level pressure (pa)
    THAIR  = SFCTMP * (SFCPRS / PAIR) ** (RAIR / CPAIR)

    !       QAIR   = Q2 / (1.0+Q2)           ! mixing ratio to specific humidity [kg/kg]
    QAIR   = Q2                       ! In WRF, driver converts to specific humidity

    EAIR   = QAIR * SFCPRS / (0.622 + 0.378 * QAIR)
    RHOAIR = (SFCPRS - 0.378 * EAIR) / (RAIR * SFCTMP)

    QPRECC = 0.10 * PRCP          ! should be from the atmospheric model
    QPRECL = 0.90 * PRCP          ! should be from the atmospheric model

    if (COSZ <= 0.0) then
       SWDOWN = 0.0
    else
       SWDOWN = SOLDN
    end if

    SOLAD(1) = SWDOWN * 0.7 * 0.5     ! direct  vis
    SOLAD(2) = SWDOWN * 0.7 * 0.5     ! direct  nir
    SOLAI(1) = SWDOWN * 0.3 * 0.5     ! diffuse vis
    SOLAI(2) = SWDOWN * 0.3 * 0.5     ! diffuse nir

  end subroutine atm


  subroutine phenology(lutyp, snowh, TV, LAT, YEARLEN, JULIAN, & !in
       LAI, SAI, TROOT, HTOP, elai, esai, IGS)
    ! vegetation phenology considering vegeation canopy being buries by snow and evolution in time
    use noahmp_veg_param, only: LK_LAI12M
    use noahmp_veg_param, only: LK_SAI12M
    use noahmp_veg_param, only: ISURBAN
    use noahmp_veg_param, only: ISWATER
    use noahmp_veg_param, only: ISBARREN
    use noahmp_veg_param, only: ISICE
    use noahmp_veg_param, only: LK_HVT
    use noahmp_veg_param, only: LK_HVB
    use noahmp_veg_param, only: LK_TMIN
    implicit none

    ! inputs
    integer                , intent(in)    :: lutyp !vegetation type
    real                   , intent(in)    :: snowh  !snow height [m]
    real                   , intent(in)    :: TV     !vegetation temperature (k)
    real                   , intent(in)    :: LAT    !latitude (radians)
    integer                , intent(in)    :: YEARLEN!Number of days in the particular year
    real                   , intent(in)    :: JULIAN !Julian day of year (fractional) ( 0 <= JULIAN < YEARLEN )
    real                   , intent(in)    :: TROOT  !root-zone averaged temperature (k)
    real                   , intent(inout) :: LAI    !LAI, unadjusted for burying by snow
    real                   , intent(inout) :: SAI    !SAI, unadjusted for burying by snow

    ! outputs
    real                   , intent(out)   :: HTOP   !top of canopy layer (m)
    real                   , intent(out)   :: elai   !leaf area index, after burying by snow
    real                   , intent(out)   :: esai   !stem area index, after burying by snow
    real                   , intent(out)   :: IGS    !growing season index (0=off, 1=on)

    ! locals

    real                                   :: DB     !thickness of canopy buried by snow (m)
    real                                   :: FB     !fraction of canopy buried by snow
    real                                   :: snowhc !critical snow depth at which short vege
    !is fully covered by snow

    integer                                :: K       !index
    integer                                :: IT1,IT2 !interpolation months
    real                                   :: DAY     !current day of year ( 0 <= DAY < YEARLEN )
    real                                   :: WT1,WT2 !interpolation weights
    real                                   :: T       !current month (1.00, ..., 12.00)

    if ( opt_veg == 1 .or. opt_veg == 3 .or. opt_veg == 4 ) then

       if (LAT >= 0.0) then
          ! Northern Hemisphere
          DAY = JULIAN
       else
          ! Southern Hemisphere.  DAY is shifted by 1/2 year.
          DAY = mod(JULIAN + (0.5 * YEARLEN), real(YEARLEN))
       end if

       T = 12.0 * DAY / real(YEARLEN)
       IT1 = T + 0.5
       IT2 = IT1 + 1
       WT1 = (IT1 + 0.5) - T
       WT2 = 1.0 - WT1
       if (IT1 <  1) IT1 = 12
       if (IT2 > 12) IT2 = 1

       LAI = WT1 * LK_LAI12M(it1,lutyp) + WT2 * LK_LAI12M(it2,lutyp)
       SAI = WT1 * LK_SAI12M(it1,lutyp) + WT2 * LK_SAI12M(it2,lutyp)
    end if
    if (SAI < 0.05) SAI = 0.0                  ! MB: SAI CHECK, change to 0.05 v3.6
    if (LAI < 0.05 .or. SAI == 0.0) LAI = 0.0  ! MB: LAI CHECK

    if ((lutyp == ISWATER) .or. (lutyp == ISBARREN) .or. (lutyp == ISICE) .or. (lutyp == ISURBAN)) then
       LAI  = 0.0
       SAI  = 0.0
    end if

    !buried by snow

    DB = min(max(snowh - LK_HVB(lutyp), 0.0), LK_HVT(lutyp) - LK_HVB(lutyp))
    FB = DB / max(1.0E-06, LK_HVT(lutyp) - LK_HVB(lutyp))

    if (LK_HVT(lutyp)> 0. .and. LK_HVT(lutyp) <= 1.0) then  !MB: change to 1.0 and 0.2 to reflect
       snowhc = LK_HVT(lutyp) * exp(-snowh / 0.2)             !      changes to HVT in MPTABLE
       FB     = min(snowh, snowhc) / snowhc
    end if

    elai =  LAI * (1.0 - FB)
    esai =  SAI * (1.0 - FB)
    if (esai < 0.05) esai = 0.0                   ! MB: esai CHECK, change to 0.05 v3.6
    if (elai < 0.05 .or. esai == 0.0) elai = 0.0  ! MB: LAI CHECK

    if (TV > LK_TMIN(lutyp)) then
       IGS = 1.0
    else
       IGS = 0.0
    end if

    HTOP = LK_HVT(lutyp)

  end subroutine phenology


  subroutine error(SWDOWN ,FSA    ,FSR    ,FIRA   ,FSH    ,FCEV   , &
       FGEV   ,FCTR   ,SSOIL  ,BEG_WB ,CANLIQ ,CANICE , &
       sneqv  ,WA     ,SMC    ,DZSNSO ,PRCP   ,ECAN   , &
       ETRAN  ,EDIR   ,runsrf ,runsub ,DT     ,NSOIL  , &
       NSNOW  ,IST    ,ERRWAT, ILOC   ,JLOC   ,fveg   , &
       SAV    ,SAG    ,FSRV   ,FSRG   ,ZWT)
    ! check surface energy balance and water balance
    implicit none
    ! inputs
    integer                        , intent(in) :: NSNOW  !maximum no. of snow layers
    integer                        , intent(in) :: NSOIL  !number of soil layers
    integer                        , intent(in) :: IST    !surface type 1->soil; 2->lake
    integer                        , intent(in) :: ILOC   !grid index
    integer                        , intent(in) :: JLOC   !grid index
    real                           , intent(in) :: SWDOWN !downward solar filtered by sun angle [w/m2]
    real                           , intent(in) :: FSA    !total absorbed solar radiation (w/m2)
    real                           , intent(in) :: FSR    !total reflected solar radiation (w/m2)
    real                           , intent(in) :: FIRA   !total net longwave rad (w/m2)  [+ to atm]
    real                           , intent(in) :: FSH    !total sensible heat (w/m2)     [+ to atm]
    real                           , intent(in) :: FCEV   !canopy evaporation heat (w/m2) [+ to atm]
    real                           , intent(in) :: FGEV   !ground evaporation heat (w/m2) [+ to atm]
    real                           , intent(in) :: FCTR   !transpiration heat flux (w/m2) [+ to atm]
    real                           , intent(in) :: SSOIL  !ground heat flux (w/m2)        [+ to soil]
    real                           , intent(in) :: fveg
    real                           , intent(in) :: SAV
    real                           , intent(in) :: SAG
    real                           , intent(in) :: FSRV
    real                           , intent(in) :: FSRG
    real                           , intent(in) :: ZWT

    real                           , intent(in) :: PRCP   !precipitation rate (kg m-2 s-1)
    real                           , intent(in) :: ECAN   !evaporation of intercepted water (mm/s)
    real                           , intent(in) :: ETRAN  !transpiration rate (mm/s)
    real                           , intent(in) :: EDIR   !soil surface evaporation rate[mm/s]
    real                           , intent(in) :: runsrf !surface runoff [mm/s]
    real                           , intent(in) :: runsub !baseflow (saturation excess) [mm/s]
    real                           , intent(in) :: CANLIQ !intercepted liquid water (mm)
    real                           , intent(in) :: CANICE !intercepted ice mass (mm)
    real                           , intent(in) :: sneqv  !snow water eqv. [mm]
    real, dimension(       1:NSOIL), intent(in) :: SMC    !soil moisture (ice + liq.) [m3/m3]
    real, dimension(-NSNOW+1:NSOIL), intent(in) :: DZSNSO !snow/soil layer thickness [m]
    real                           , intent(in) :: WA     !water storage in aquifer [mm]
    real                           , intent(in) :: DT     !time step [sec]
    real                           , intent(in) :: BEG_WB !water storage at begin of a timesetp [mm]
    real                           , intent(out) :: ERRWAT !error in water balance [mm/timestep]

    integer                                     :: IZ     !do-loop index
    real                                        :: END_WB !water storage at end of a timestep [mm]
    !KWM REAL                                        :: ERRWAT !error in water balance [mm/timestep]
    real                                        :: ERRENG !error in surface energy balance [w/m2]
    real                                        :: ERRSW  !error in shortwave radiation balance [w/m2]
    real                                        :: FSRVG
    character(256)                              :: message

    !jref:start
    ERRSW   = SWDOWN - (FSA + FSR)
    !   ERRSW   = SWDOWN - (SAV+SAG + FSRV+FSRG)
    !   WRITE(*,*) "ERRSW =",ERRSW
    if (abs(ERRSW) > 0.01) then            ! w/m2
       write(*,*) "VEGETATION!"
       write(*,*) "SWDOWN*fveg =",SWDOWN*fveg
       write(*,*) "fveg*(SAV+SAG) =",fveg*SAV + SAG
       write(*,*) "fveg*(FSRV +FSRG)=",fveg*FSRV + FSRG
       write(*,*) "GROUND!"
       write(*,*) "(1-.fveg)*SWDOWN =",(1.-fveg)*SWDOWN
       write(*,*) "(1.-fveg)*SAG =",(1.-fveg)*SAG
       write(*,*) "(1.-fveg)*FSRG=",(1.-fveg)*FSRG
       write(*,*) "FSRV   =",FSRV
       write(*,*) "FSRG   =",FSRG
       write(*,*) "FSR    =",FSR
       write(*,*) "SAV    =",SAV
       write(*,*) "SAG    =",SAG
       write(*,*) "FSA    =",FSA
       !jref:end
       write(message,*) 'ERRSW =',ERRSW
       call wrf_message(trim(message))
       call wrf_error_fatal("Stop in Noah-MP")
    end if

    ERRENG = SAV+SAG-(FIRA+FSH+FCEV+FGEV+FCTR+SSOIL)
    !   ERRENG = fveg*SAV+SAG-(FIRA+FSH+FCEV+FGEV+FCTR+SSOIL)
    !   WRITE(*,*) "ERRENG =",ERRENG
    if (abs(ERRENG) > 0.01) then
       write(message,*) 'ERRENG =',ERRENG
       call wrf_message(trim(message))
       write(message,'(i6,1x,i6,1x,7F10.4)')ILOC,JLOC,FSA,FIRA,FSH,FCEV,FGEV,FCTR,SSOIL
       call wrf_message(trim(message))
       call wrf_error_fatal("energy budget problem in NOAHMP LSM")
    end if

    if (IST == 1) then                                       !soil
       END_WB = CANLIQ + CANICE + sneqv + WA
       do IZ = 1,NSOIL
          END_WB = END_WB + SMC(IZ) * DZSNSO(IZ) * 1000.0
       end do
       ERRWAT = END_WB - BEG_WB - (PRCP - ECAN - ETRAN - EDIR - runsrf - runsub) * DT
    else                 !KWM
       ERRWAT = 0.0      !KWM
    end if
  end subroutine error


  subroutine energy(sltyp, lutyp, ICE, IST, ISC, nsnow, nsoil, & !in
       ISNOW, DT     ,RHOAIR ,SFCPRS ,QAIR   , & !in
       SFCTMP ,THAIR  ,LWDN   ,UU     ,VV     ,ZREF   , & !in
       CO2AIR ,O2AIR  ,SOLAD  ,SOLAI  ,COSZ   ,IGS    , & !in
       EAIR   ,HTOP   ,TBOT   ,ZBOT   ,ZSNSO  ,ZSOIL  , & !in
       elai   ,esai   ,CSOIL  ,FWET   ,FOLN   ,         & !in
       fveg   ,                                         & !in
       QSNOW  ,DZSNSO ,LAT    ,CANLIQ ,CANICE ,ILOC   , JLOC, & !in
       IMELT  ,SNICEV ,SNLIQV ,epore  ,T2M    ,FSNO   , & !out
       SAV    ,SAG    ,QMELT  ,FSA    ,FSR    ,TAUX   , & !out
       TAUY   ,FIRA   ,FSH    ,FCEV   ,FGEV   ,FCTR   , & !out
       TRAD   ,PSN    ,APAR   ,SSOIL  ,BTRANI ,BTRAN  , & !out
       PONDING,TS     ,LATHEAV , LATHEAG , frozen_canopy,frozen_ground,                       & !out
       TV     ,TG     ,STC    ,snowh  ,EAH    ,TAH    , & !inout
       SNEQVO ,sneqv  ,soilwat   ,SMC    ,SNICE  ,SNLIQ  , & !inout
       ALBOLD ,CM     ,CH     ,DX     ,DZ8W   ,Q2     , &   !inout
       TAUSS  ,                                         & !inout
                                !jref:start
       QC     ,PBLH   ,QSFC   ,PSFC, IZ0TLND, & !in
       T2MV   ,T2MB   ,FSRV   , &
       FSRG   ,RSSUN  ,RSSHA  ,BGAP   ,WGAP,TGV,TGB,&
       Q1     ,Q2V    ,Q2B    ,Q2E    ,CHV  ,CHB, EMISSI,&
       SHG,SHC,SHB,EVG,EVB,GHV,GHB,IRG,IRC,IRB,TR,EVC,CHLEAF,CHUC,CHV2,CHB2 )   !out
    ! we use different approaches to deal with subgrid features of radiation transfer and turbulent
    ! transfer. We use 'tile' approach to compute turbulent fluxes, while we use modified two-
    ! stream to compute radiation transfer. Tile approach, assemblying vegetation canopies together,
    ! may expose too much ground surfaces (either covered by snow or grass) to solar radiation. The
    ! modified two-stream assumes vegetation covers fully the gridcell but with gaps between tree
    ! crowns.
    !
    ! turbulence transfer : 'tile' approach to compute energy fluxes in vegetated fraction and
    !                         bare fraction separately and then sum them up weighted by fraction
    !                     --------------------------------------
    !                    / O  O  O  O  O  O  O  O  /          /
    !                   /  |  |  |  |  |  |  |  | /          /
    !                  / O  O  O  O  O  O  O  O  /          /
    !                 /  |  |  |tile1|  |  |  | /  tile2   /
    !                / O  O  O  O  O  O  O  O  /  bare    /
    !               /  |  |  | vegetated |  | /          /
    !              / O  O  O  O  O  O  O  O  /          /
    !             /  |  |  |  |  |  |  |  | /          /
    !            --------------------------------------
    !
    ! radiation transfer : modified two-stream (Yang and Friedl, 2003, JGR; Niu ang Yang, 2004, JGR)
    !                     --------------------------------------  two-stream treats leaves as
    !                    /   O   O   O   O   O   O   O   O    /  cloud over the entire grid-cell,
    !                   /    |   |   |   |   |   |   |   |   / while the modified two-stream
    !                  /   O   O   O   O   O   O   O   O    / aggregates cloudy leaves into
    !                 /    |   |   |   |   |   |   |   |   / tree crowns with gaps (as shown in
    !                /   O   O   O   O   O   O   O   O    / the left figure). We assume these
    !               /    |   |   |   |   |   |   |   |   / tree crowns are evenly distributed
    !              /   O   O   O   O   O   O   O   O    / within the gridcell with 100% veg
    !             /    |   |   |   |   |   |   |   |   / fraction, but with gaps. The 'tile'
    !            -------------------------------------- approach overlaps too much shadows.
    !
    use noahmp_gen_param, only: KK_MLTFCT
    use noahmp_gen_param, only: KK_Z0SNO
    use noahmp_gen_param, only: KK_EMSSOIL
    use noahmp_gen_param, only: KK_EMSLAKE
    use noahmp_soil_param, only: LK_BEXP
    use noahmp_soil_param, only: LK_PSISAT
    use noahmp_soil_param, only: LK_SMCMAX
    use noahmp_soil_param, only: LK_SMCREF
    use noahmp_soil_param, only: LK_SMCWLT
    use noahmp_veg_param, only: ISURBAN
    use noahmp_veg_param, only: LK_NROOT
    use noahmp_veg_param, only: LK_Z0MVT
    use noahmp_veg_param, only: LK_CWPVT
    implicit none
    ! inputs
    integer                           , intent(in)    :: ILOC
    integer                           , intent(in)    :: JLOC
    integer                           , intent(in)    :: ICE    !ice (ice = 1)
    integer, intent(in) :: sltyp
    integer                           , intent(in)    :: lutyp !vegetation physiology type
    integer                           , intent(in)    :: IST    !surface type: 1->soil; 2->lake
    integer                           , intent(in)    :: ISC    !soil color type (1-lighest; 8-darkest)
    integer                           , intent(in)    :: NSNOW  !maximum no. of snow layers
    integer                           , intent(in)    :: NSOIL  !number of soil layers
    integer                           , intent(in)    :: ISNOW  !actual no. of snow layers
    real                              , intent(in)    :: DT     !time step [sec]
    real                              , intent(in)    :: QSNOW  !snowfall on the ground (mm/s)
    real                              , intent(in)    :: RHOAIR !density air (kg/m3)
    real                              , intent(in)    :: EAIR   !vapor pressure air (pa)
    real                              , intent(in)    :: SFCPRS !pressure (pa)
    real                              , intent(in)    :: QAIR   !specific humidity (kg/kg)
    real                              , intent(in)    :: SFCTMP !air temperature (k)
    real                              , intent(in)    :: THAIR  !potential temperature (k)
    real                              , intent(in)    :: LWDN   !downward longwave radiation (w/m2)
    real                              , intent(in)    :: UU     !wind speed in e-w dir (m/s)
    real                              , intent(in)    :: VV     !wind speed in n-s dir (m/s)
    real   , dimension(       1:    2), intent(in)    :: SOLAD  !incoming direct solar rad. (w/m2)
    real   , dimension(       1:    2), intent(in)    :: SOLAI  !incoming diffuse solar rad. (w/m2)
    real                              , intent(in)    :: COSZ   !cosine solar zenith angle (0-1)
    real                              , intent(in)    :: elai   !LAI adjusted for burying by snow
    real                              , intent(in)    :: esai   !LAI adjusted for burying by snow
    real                              , intent(in)    :: CSOIL  !vol. soil heat capacity [j/m3/k]
    real                              , intent(in)    :: FWET   !fraction of canopy that is wet [-]
    real                              , intent(in)    :: HTOP   !top of canopy layer (m)
    real                              , intent(in)    :: fveg   !greeness vegetation fraction (-)
    real                              , intent(in)    :: LAT    !latitude (radians)
    real                              , intent(in)    :: CANLIQ !canopy-intercepted liquid water (mm)
    real                              , intent(in)    :: CANICE !canopy-intercepted ice mass (mm)
    real                              , intent(in)    :: FOLN   !foliage nitrogen (%)
    real                              , intent(in)    :: CO2AIR !atmospheric co2 concentration (pa)
    real                              , intent(in)    :: O2AIR  !atmospheric o2 concentration (pa)
    real                              , intent(in)    :: IGS    !growing season index (0=off, 1=on)

    real                              , intent(in)    :: ZREF   !reference height (m)
    real                              , intent(in)    :: TBOT   !bottom condition for soil temp. (k)
    real                              , intent(in)    :: ZBOT   !depth for TBOT [m]
    real   , dimension(-NSNOW+1:NSOIL), intent(in)    :: ZSNSO  !layer-bottom depth from snow surf [m]
    real   , dimension(       1:NSOIL), intent(in)    :: ZSOIL  !layer-bottom depth from soil surf [m]
    real   , dimension(-NSNOW+1:NSOIL), intent(in)    :: DZSNSO !depth of snow & soil layer-bottom [m]

    !jref:start; in
    integer                           , intent(in)    :: IZ0TLND
    real                              , intent(in)    :: QC     !cloud water mixing ratio
    real                              , intent(in)    :: PBLH   !planetary boundary layer height
    real                              , intent(inout) :: QSFC   !mixing ratio at lowest model layer
    real                              , intent(in)    :: PSFC   !pressure at lowest model layer
    real                              , intent(in)    :: DX     !horisontal resolution
    real                              , intent(in)    :: DZ8W   !thickness of lowest layer
    real                              , intent(in)    :: Q2     !mixing ratio (kg/kg)
    !jref:end

    ! outputs
    integer, dimension(-NSNOW+1:NSOIL), intent(out)   :: IMELT  !phase change index [1-melt; 2-freeze]
    real   , dimension(-NSNOW+1:    0), intent(out)   :: SNICEV !partial volume ice [m3/m3]
    real   , dimension(-NSNOW+1:    0), intent(out)   :: SNLIQV !partial volume liq. water [m3/m3]
    real   , dimension(-NSNOW+1:    0), intent(out)   :: epore  !effective porosity [m3/m3]
    real                              , intent(out)   :: FSNO   !snow cover fraction (-)
    real                              , intent(out)   :: QMELT  !snowmelt [mm/s]
    real                              , intent(out)   :: PONDING!pounding at ground [mm]
    real                              , intent(out)   :: SAV    !solar rad. absorbed by veg. (w/m2)
    real                              , intent(out)   :: SAG    !solar rad. absorbed by ground (w/m2)
    real                              , intent(out)   :: FSA    !tot. absorbed solar radiation (w/m2)
    real                              , intent(out)   :: FSR    !tot. reflected solar radiation (w/m2)
    real                              , intent(out)   :: TAUX   !wind stress: e-w (n/m2)
    real                              , intent(out)   :: TAUY   !wind stress: n-s (n/m2)
    real                              , intent(out)   :: FIRA   !total net LW. rad (w/m2)   [+ to atm]
    real                              , intent(out)   :: FSH    !total sensible heat (w/m2) [+ to atm]
    real                              , intent(out)   :: FCEV   !canopy evaporation (w/m2)  [+ to atm]
    real                              , intent(out)   :: FGEV   !ground evaporation (w/m2)  [+ to atm]
    real                              , intent(out)   :: FCTR   !transpiration (w/m2)       [+ to atm]
    real                              , intent(out)   :: TRAD   !radiative temperature (k)
    real                              , intent(out)   :: T2M    !2 m height air temperature (k)
    real                              , intent(out)   :: PSN    !total photosyn. (umolco2/m2/s) [+]
    real                              , intent(out)   :: APAR   !total photosyn. active energy (w/m2)
    real                              , intent(out)   :: SSOIL  !ground heat flux (w/m2)   [+ to soil]
    real   , dimension(       1:NSOIL), intent(out)   :: BTRANI !soil water transpiration factor (0-1)
    real                              , intent(out)   :: BTRAN  !soil water transpiration factor (0-1)
    !  REAL                              , INTENT(out)   :: LATHEA !latent heat vap./sublimation (j/kg)
    real                              , intent(out)   :: LATHEAV !latent heat vap./sublimation (j/kg)
    real                              , intent(out)   :: LATHEAG !latent heat vap./sublimation (j/kg)
    logical                           , intent(out)   :: FROZEN_GROUND ! used to define latent heat pathway
    logical                           , intent(out)   :: FROZEN_CANOPY ! used to define latent heat pathway

    !jref:start
    real                              , intent(out)   :: FSRV    !veg. reflected solar radiation (w/m2)
    real                              , intent(out)   :: FSRG    !ground reflected solar radiation (w/m2)
    real, intent(out) :: RSSUN        !sunlit leaf stomatal resistance (s/m)
    real, intent(out) :: RSSHA        !shaded leaf stomatal resistance (s/m)
    !jref:end - out for debug

    !jref:start; output
    real                              , intent(out)   :: T2MV   !2-m air temperature over vegetated part [k]
    real                              , intent(out)   :: T2MB   !2-m air temperature over bare ground part [k]
    real                              , intent(out)   :: BGAP
    real                              , intent(out)   :: WGAP
    !jref:end

    ! input & output
    real                              , intent(inout) :: TS     !surface temperature (k)
    real                              , intent(inout) :: TV     !vegetation temperature (k)
    real                              , intent(inout) :: TG     !ground temperature (k)
    real   , dimension(-NSNOW+1:NSOIL), intent(inout) :: STC    !snow/soil temperature [k]
    real                              , intent(inout) :: snowh  !snow height [m]
    real                              , intent(inout) :: sneqv  !snow mass (mm)
    real                              , intent(inout) :: SNEQVO !snow mass at last time step (mm)
    real   , dimension(       1:NSOIL), intent(inout) :: soilwat   !liquid soil moisture [m3/m3]
    real   , dimension(       1:NSOIL), intent(inout) :: SMC    !soil moisture (ice + liq.) [m3/m3]
    real   , dimension(-NSNOW+1:    0), intent(inout) :: SNICE  !snow ice mass (kg/m2)
    real   , dimension(-NSNOW+1:    0), intent(inout) :: SNLIQ  !snow liq mass (kg/m2)
    real                              , intent(inout) :: EAH    !canopy air vapor pressure (pa)
    real                              , intent(inout) :: TAH    !canopy air temperature (k)
    real                              , intent(inout) :: ALBOLD !snow albedo at last time step(CLASS type)
    real                              , intent(inout) :: TAUSS  !non-dimensional snow age
    real                              , intent(inout) :: CM     !momentum drag coefficient
    real                              , intent(inout) :: CH     !sensible heat exchange coefficient
    real                              , intent(inout) :: Q1
    !  REAL                                              :: Q2E
    real,                               intent(out)   :: EMISSI

    ! local
    integer                                           :: IZ     !do-loop index
    logical                                           :: VEG    !true if vegetated surface
    real                                              :: UR     !wind speed at height ZLVL (m/s)
    real                                              :: ZLVL   !reference height (m)
    real                                              :: FSUN   !sunlit fraction of canopy [-]
    real                                              :: RB     !leaf boundary layer resistance (s/m)
    real                                              :: RSURF  !ground surface resistance (s/m)
    real                                              :: L_RSURF!Dry-layer thickness for computing RSURF (Sakaguchi and Zeng, 2009)
    real                                              :: D_RSURF!Reduced vapor diffusivity in soil for computing RSURF (SZ09)
    real                                              :: BEVAP  !soil water evaporation factor (0- 1)
    real                                              :: MOL    !Monin-Obukhov length (m)
    real                                              :: VAI    !sum of LAI  + stem area index [m2/m2]
    real                                              :: CWP    !canopy wind extinction parameter
    real                                              :: ZPD    !zero plane displacement (m)
    real                                              :: Z0M    !z0 momentum (m)
    real                                              :: ZPDG   !zero plane displacement (m)
    real                                              :: Z0MG   !z0 momentum, ground (m)
    real                                              :: EMV    !vegetation emissivity
    real                                              :: EMG    !ground emissivity
    real                                              :: FIRE   !emitted IR (w/m2)

    real                                              :: LAISUN !sunlit leaf area index (m2/m2)
    real                                              :: LAISHA !shaded leaf area index (m2/m2)
    real                                              :: PSNSUN !sunlit photosynthesis (umolco2/m2/s)
    real                                              :: PSNSHA !shaded photosynthesis (umolco2/m2/s)
    !jref:start - for debug
    !  REAL                                              :: RSSUN  !sunlit stomatal resistance (s/m)
    !  REAL                                              :: RSSHA  !shaded stomatal resistance (s/m)
    !jref:end - for debug
    real                                              :: PARSUN !par absorbed per sunlit LAI (w/m2)
    real                                              :: PARSHA !par absorbed per shaded LAI (w/m2)

    real, dimension(-NSNOW+1:NSOIL)                   :: FACT   !temporary used in phase change
    real, dimension(-NSNOW+1:NSOIL)                   :: DF     !thermal conductivity [w/m/k]
    real, dimension(-NSNOW+1:NSOIL)                   :: HCPCT  !heat capacity [j/m3/k]
    real                                              :: BDSNO  !bulk density of snow (kg/m3)
    real                                              :: FMELT  !melting factor for snow cover frac
    real                                              :: GX     !temporary variable
    real, dimension(-NSNOW+1:NSOIL)                   :: PHI    !light through water (w/m2)
    !  REAL                                              :: GAMMA  !psychrometric constant (pa/k)
    real                                              :: GAMMAV  !psychrometric constant (pa/k)
    real                                              :: GAMMAG  !psychrometric constant (pa/k)
    real                                              :: PSI    !surface layer soil matrix potential (m)
    real                                              :: RHSUR  !raltive humidity in surface soil/snow air space (-)

    ! temperature and fluxes over vegetated fraction

    real                                              :: TAUXV  !wind stress: e-w dir [n/m2]
    real                                              :: TAUYV  !wind stress: n-s dir [n/m2]
    real,intent(out)                                              :: IRC    !canopy net LW rad. [w/m2] [+ to atm]
    real,intent(out)                                              :: IRG    !ground net LW rad. [w/m2] [+ to atm]
    real,intent(out)                                              :: SHC    !canopy sen. heat [w/m2]   [+ to atm]
    real,intent(out)                                              :: SHG    !ground sen. heat [w/m2]   [+ to atm]
    !jref:start
    real,intent(out)                                  :: Q2V
    real,intent(out)                                  :: Q2B
    real,intent(out)                                  :: Q2E
    !jref:end
    real,intent(out)                                              :: EVC    !canopy evap. heat [w/m2]  [+ to atm]
    real,intent(out)                                              :: EVG    !ground evap. heat [w/m2]  [+ to atm]
    real,intent(out)                                              :: TR     !transpiration heat [w/m2] [+ to atm]
    real,intent(out)                                              :: GHV    !ground heat flux [w/m2]  [+ to soil]
    real,intent(out)                                  :: TGV    !ground surface temp. [k]
    real                                              :: CMV    !momentum drag coefficient
    real,intent(out)                                  :: CHV    !sensible heat exchange coefficient

    ! temperature and fluxes over bare soil fraction

    real                                              :: TAUXB  !wind stress: e-w dir [n/m2]
    real                                              :: TAUYB  !wind stress: n-s dir [n/m2]
    real,intent(out)                                              :: IRB    !net longwave rad. [w/m2] [+ to atm]
    real,intent(out)                                              :: SHB    !sensible heat [w/m2]     [+ to atm]
    real,intent(out)                                              :: EVB    !evaporation heat [w/m2]  [+ to atm]
    real,intent(out)                                              :: GHB    !ground heat flux [w/m2] [+ to soil]
    real,intent(out)                                  :: TGB    !ground surface temp. [k]
    real                                              :: CMB    !momentum drag coefficient
    real,intent(out)                                  :: CHB    !sensible heat exchange coefficient
    real,intent(out)                                  :: CHLEAF !leaf exchange coefficient
    real,intent(out)                                  :: CHUC   !under canopy exchange coefficient
    !jref:start
    real,intent(out)                                  :: CHV2    !sensible heat conductance, canopy air to ZLVL air (m/s)
    real,intent(out)                                  :: CHB2    !sensible heat conductance, canopy air to ZLVL air (m/s)
    real                                  :: noahmpres

    !jref:end

    real, parameter                   :: MPE    = 1.E-6
    real, parameter                   :: PSIWLT = -150.  !metric potential for wilting point (m)
    real, parameter                   :: Z0     = 0.01   ! Bare-soil roughness length (m) (i.e., under the canopy)

    ! initialize fluxes from veg. fraction
    TAUXV     = 0.0
    TAUYV     = 0.0
    IRC       = 0.0
    SHC       = 0.0
    IRG       = 0.0
    SHG       = 0.0
    EVG       = 0.0
    EVC       = 0.0
    TR        = 0.0
    GHV       = 0.0
    PSNSUN    = 0.0
    PSNSHA    = 0.0
    T2MV      = 0.0
    Q2V       = 0.0
    CHV       = 0.0
    CHLEAF    = 0.0
    CHUC      = 0.0
    CHV2      = 0.0

    ! wind speed at reference height: ur >= 1
    ur = max(sqrt(uu * uu + vv * vv), 1.0)

    ! vegetated or non-vegetated
    VAI = elai + esai
    VEG = .false.
    if (VAI > 0.0) VEG = .true.

    ! ground snow cover fraction [Niu and Yang, 2007, JGR]
    FSNO = 0.0
    if (snowh > 0.0)  then
       BDSNO    = sneqv / snowh
       FMELT    = (BDSNO / 100.0) ** KK_MLTFCT
       FSNO     = tanh(snowh /(2.5 * Z0 * FMELT))
    end if

    ! ground roughness length
    if (IST == 2) then
       if (TG <= TFRZ) then
          Z0MG = 0.01 * (1.0 - FSNO) + FSNO * KK_Z0SNO
       else
          Z0MG = 0.01
       end if
    else
       Z0MG = Z0 * (1.0-FSNO) + FSNO * KK_Z0SNO
    end if

    ! roughness length and displacement height

    ZPDG  = snowh
    if (VEG) then
       Z0M  = LK_Z0MVT(lutyp)
       ZPD  = 0.65 * HTOP
       if (snowh > ZPD) ZPD  = snowh
    else
       Z0M  = Z0MG
       ZPD  = ZPDG
    end if

    ZLVL = max(ZPD, HTOP) + ZREF
    if (ZPDG >= ZLVL) ZLVL = ZPDG + ZREF
    !     UR   = UR*LOG(ZLVL/Z0M)/LOG(10./Z0M)       !input UR is at 10m

    ! canopy wind absorption coeffcient
    CWP = LK_CWPVT(lutyp)

    ! Thermal properties of soil, snow, lake, and frozen soil
    call thermoprop(sltyp, NSOIL   ,NSNOW   ,ISNOW   ,IST     ,DZSNSO  , & !in
         DT      ,snowh   ,SNICE   ,SNLIQ   ,CSOIL   , & !in
         SMC     ,soilwat    ,TG      ,STC     ,UR      , & !in
         LAT     ,Z0M     ,ZLVL    ,lutyp, & !in
         DF      ,HCPCT   ,SNICEV  ,SNLIQV  ,epore   , & !out
         FACT    )                              !out

    ! Solar radiation: absorbed & reflected by the ground and canopy
    call  radiation(lutyp  ,IST     ,ISC     ,ICE     ,NSOIL   , & !in
         SNEQVO  ,sneqv   ,DT      ,COSZ    ,snowh   , & !in
         TG      ,TV      ,FSNO    ,QSNOW   ,FWET    , & !in
         elai    ,esai    ,SMC     ,SOLAD   ,SOLAI   , & !in
         fveg    ,ILOC    ,JLOC    ,                   & !in
         ALBOLD  ,TAUSS   ,                            & !inout
         FSUN    ,LAISUN  ,LAISHA  ,PARSUN  ,PARSHA  , & !out
         SAV     ,SAG     ,FSR     ,FSA     ,FSRV    , &
         FSRG    ,BGAP    ,WGAP    )            !out

    ! vegetation and ground emissivity
    EMV = 1.0 - exp(-(elai + esai) / 1.0)
    if (ICE == 1) then
       EMG = 0.98 * (1.0 - FSNO) + 1.0 * FSNO
    elseif (ist == 1) then
       EMG = KK_EMSSOIL * (1.0 - FSNO) + 1.0 * FSNO
    else
       EMG = KK_EMSLAKE * (1.0 - FSNO) + 1.0 * FSNO
    end if

    ! soil moisture factor controlling stomatal resistance

    BTRAN = 0.0

    if (IST ==1 ) then
       do IZ = 1, LK_NROOT(lutyp)
          if (opt_btr == 1) then                  ! Noah
             GX    = (soilwat(IZ) - LK_SMCWLT(sltyp)) / (LK_SMCREF(sltyp) - (LK_SMCWLT(sltyp)))
          end if
          if (opt_btr == 2) then                  ! CLM
             PSI   = max(PSIWLT, -LK_PSISAT(sltyp) * (max(0.01, soilwat(IZ)) / LK_SMCMAX(sltyp)) ** (-LK_BEXP(sltyp)))
             GX    = (1.0 - PSI / PSIWLT) / (1.0 + LK_PSISAT(sltyp) / PSIWLT)
          end if
          if (opt_btr == 3) then                  ! SSiB
             PSI   = max(PSIWLT, -LK_PSISAT(sltyp) * (max(0.01, soilwat(IZ)) / LK_SMCMAX(sltyp)) ** (-LK_BEXP(sltyp)))
             GX    = 1.0 - exp(-5.8 * (log(PSIWLT / PSI)))
          end if

          GX = min(1.0, max(0.0, GX))
          BTRANI(IZ) = max(MPE, DZSNSO(IZ) / (-ZSOIL(LK_NROOT(lutyp))) * GX)
          BTRAN      = BTRAN + BTRANI(IZ)
       end do
       BTRAN = max(MPE, BTRAN)

       BTRANI(1:LK_NROOT(lutyp)) = BTRANI(1:LK_NROOT(lutyp)) / BTRAN
    end if

    ! soil surface resistance for ground evap.
    BEVAP = max(0.0, soilwat(1) / LK_SMCMAX(sltyp))
    if (IST == 2) then
       RSURF = 1.0         ! avoid being divided by 0
       RHSUR = 1.0
    else

       ! RSURF based on Sakaguchi and Zeng, 2009
       ! taking the "residual water content" to be the wilting point,
       ! and correcting the exponent on the D term (typo in SZ09 ?)
       L_RSURF = (-ZSOIL(1)) * (exp((1.0 - min(1.0, soilwat(1) / LK_SMCMAX(sltyp))) ** 5) - 1.0 ) / (2.71828 - 1.0)
       D_RSURF = 2.2E-5 * LK_SMCMAX(sltyp) * LK_SMCMAX(sltyp) &
            & * (1.0 - LK_SMCWLT(sltyp) / LK_SMCMAX(sltyp)) ** (2.0 + 3.0 / LK_BEXP(sltyp))
       RSURF = L_RSURF / D_RSURF

       ! Older RSURF computations:
       !    RSURF = FSNO * 1. + (1.-FSNO)* EXP(8.25-4.225*BEVAP) !Sellers (1992)
       !    RSURF = FSNO * 1. + (1.-FSNO)* EXP(8.25-6.0  *BEVAP) !adjusted to decrease RSURF for wet soil

       if (soilwat(1) < 0.01 .and. snowh == 0.0) RSURF = 1.0E6
       PSI   = -LK_PSISAT(sltyp) * (max(0.01, soilwat(1)) / LK_SMCMAX(sltyp)) ** (-LK_BEXP(sltyp))
       RHSUR = FSNO + (1.0 - FSNO) * exp(PSI * GRAV / (RVAP * TG))
    end if

    ! urban - jref
    if (lutyp == ISURBAN .and. snowh == 0.0) then
       RSURF = 1.0E6
    end if

    ! set psychrometric constant

    if (TV > TFRZ) then           ! Barlage: add distinction between ground and
       LATHEAV = HVAP                ! vegetation in v3.6
       frozen_canopy = .false.
    else
       LATHEAV = HSUB
       frozen_canopy = .true.
    end if
    GAMMAV = CPAIR * SFCPRS / (0.622 * LATHEAV)

    if (TG > TFRZ) then
       LATHEAG = HVAP
       frozen_ground = .false.
    else
       LATHEAG = HSUB
       frozen_ground = .true.
    end if
    GAMMAG = CPAIR * SFCPRS / (0.622 * LATHEAG)

    !     IF (SFCTMP > TFRZ) THEN
    !        LATHEA = HVAP
    !     ELSE
    !        LATHEA = HSUB
    !     END IF
    !     GAMMA = CPAIR*SFCPRS/(0.622*LATHEA)

    ! Surface temperatures of the ground and canopy and energy fluxes

    if (VEG .and. fveg > 0) then
       TGV = TG
       CMV = CM
       CHV = CH
       call vege_flux(NSNOW   ,NSOIL   ,ISNOW   ,lutyp  ,VEG     , & !in
            DT      ,SAV     ,SAG     ,LWDN    ,UR      , & !in
            UU      ,VV      ,SFCTMP  ,THAIR   ,QAIR    , & !in
            EAIR    ,RHOAIR  ,snowh   ,VAI     ,GAMMAV   ,GAMMAG   , & !in
            FWET    ,LAISUN  ,LAISHA  ,CWP     ,DZSNSO  , & !in
            HTOP    ,ZLVL    ,ZPD     ,Z0M     ,fveg    , & !in
            Z0MG    ,EMV     ,EMG     ,CANLIQ           , & !in
            CANICE  ,STC     ,DF      ,RSSUN   ,RSSHA   , & !in
            RSURF   ,LATHEAV ,LATHEAG ,PARSUN  ,PARSHA  ,IGS     , & !in
            FOLN    ,CO2AIR  ,O2AIR   ,BTRAN   ,SFCPRS  , & !in
            RHSUR   ,ILOC    ,JLOC    ,Q2      , & !in
            EAH     ,TAH     ,TV      ,TGV     ,CMV     , & !inout
            CHV     ,DX      ,DZ8W    ,                   & !inout
            TAUXV   ,TAUYV   ,IRG     ,IRC     ,SHG     , & !out
            SHC     ,EVG     ,EVC     ,TR      ,GHV     , & !out
            T2MV    ,PSNSUN  ,PSNSHA  ,                   & !out
                                !jref:start
            QC      ,PBLH    ,QSFC    ,PSFC, & !in
            IZ0TLND ,Q2V     ,CHV2, CHLEAF, CHUC)               !inout
       !jref:end
    end if

    TGB = TG
    CMB = CM
    CHB = CH
    call bare_flux(lutyp, dt, NSNOW, nsoil, ISNOW, SAG     , & !in
         LWDN    ,UR      ,UU      ,VV      ,SFCTMP  , & !in
         THAIR   ,QAIR    ,EAIR    ,RHOAIR  ,snowh   , & !in
         DZSNSO  ,ZLVL    ,ZPDG    ,Z0MG    ,          & !in
         EMG     ,STC     ,DF      ,RSURF   ,LATHEAG  , & !in
         GAMMAG   ,RHSUR   ,ILOC    ,JLOC    ,Q2      , & !in
         TGB     ,CMB     ,CHB     ,                   & !inout
         TAUXB   ,TAUYB   ,IRB     ,SHB     ,EVB     , & !out
         GHB     ,T2MB    ,DX      ,DZ8W, & !out
         QC      ,PBLH    ,QSFC    ,PSFC, & !in
         IZ0TLND ,SFCPRS  ,Q2B,   CHB2)                          !in
    !jref:end

    !energy balance at vege canopy: SAV          =(IRC+SHC+EVC+TR)     *fveg  at   fveg
    !energy balance at vege ground: SAG*    fveg =(IRG+SHG+EVG+GHV)    *fveg  at   fveg
    !energy balance at bare ground: SAG*(1.-fveg)=(IRB+SHB+EVB+GHB)*(1.-fveg) at 1-fveg

    if (VEG .and. fveg > 0) then
       TAUX  = fveg * TAUXV     + (1.0 - fveg) * TAUXB
       TAUY  = fveg * TAUYV     + (1.0 - fveg) * TAUYB
       FIRA  = fveg * IRG       + (1.0 - fveg) * IRB       + IRC
       FSH   = fveg * SHG       + (1.0 - fveg) * SHB       + SHC
       FGEV  = fveg * EVG       + (1.0 - fveg) * EVB
       SSOIL = fveg * GHV       + (1.0 - fveg) * GHB
       FCEV  = EVC
       FCTR  = TR
       TG    = fveg * TGV       + (1.0 - fveg) * TGB
       T2M   = fveg * T2MV      + (1.0 - fveg) * T2MB
       TS    = fveg * TV        + (1.0 - fveg) * TGB
       CM    = fveg * CMV       + (1.0 - fveg) * CMB      ! better way to average?
       CH    = fveg * CHV       + (1.0 - fveg) * CHB
       Q1    = fveg * (EAH*0.622/(SFCPRS - 0.378*EAH)) + (1.0 - fveg)*QSFC
       Q2E   = fveg * Q2V       + (1.0 - fveg) * Q2B
    else
       TAUX  = TAUXB
       TAUY  = TAUYB
       FIRA  = IRB
       FSH   = SHB
       FGEV  = EVB
       SSOIL = GHB
       TG    = TGB
       T2M   = T2MB
       FCEV  = 0.
       FCTR  = 0.
       TS    = TG
       CM    = CMB
       CH    = CHB
       Q1    = QSFC
       Q2E   = Q2B
       RSSUN = 0.0
       RSSHA = 0.0
       TGV   = TGB
       CHV   = CHB
    end if

    FIRE = LWDN + FIRA

    if (FIRE <= 0.0) then
       write(6,*) 'emitted longwave <0; skin T may be wrong due to inconsistent'
       write(6,*) 'input of SHDFAC with LAI'
       write(6,*) ILOC, JLOC, 'SHDFAC=',fveg,'VAI=',VAI,'TV=',TV,'TG=',TG
       write(6,*) 'LWDN=',LWDN,'FIRA=',FIRA,'SNOWH=',snowh
       call wrf_error_fatal("STOP in Noah-MP")
    end if

    ! Compute a net emissivity
    EMISSI = fveg * ( EMG*(1-EMV) + EMV + EMV*(1-EMV)*(1-EMG) ) + &
         (1-fveg) * EMG

    ! When we're computing a TRAD, subtract from the emitted IR the
    ! reflected portion of the incoming LWDN, so we're just
    ! considering the IR originating in the canopy/ground system.

    TRAD = ((FIRE - (1.0 - EMISSI) * LWDN ) / (EMISSI * SB) ) ** 0.25

    ! Old TRAD calculation not taking into account Emissivity:
    ! TRAD = (FIRE/SB)**0.25

    APAR = PARSUN * LAISUN + PARSHA * LAISHA
    PSN  = PSNSUN * LAISUN + PSNSHA * LAISHA

    ! 3L snow & 4L soil temperatures
    call tsnosoi(ICE     ,NSOIL   ,NSNOW   ,ISNOW   ,IST     , & !in
         TBOT    ,ZSNSO   ,SSOIL   ,DF      ,HCPCT   , & !in
         ZBOT    ,SAG     ,DT      ,snowh   ,DZSNSO  , & !in
         TG      ,ILOC    ,JLOC    ,                   & !in
         STC     )                                       !inout

    ! adjusting snow surface temperature
    if (opt_stc == 2) then
       if (snowh > 0.05 .and. TG > TFRZ) then
          TGV = TFRZ
          TGB = TFRZ
          if (VEG .and. fveg > 0) then
             TG    = fveg * TGV       + (1.0 - fveg) * TGB
             TS    = fveg * TV        + (1.0 - fveg) * TGB
          else
             TG    = TGB
             TS    = TGB
          end if
       end if
    end if

    ! energy released or consumed by snow & frozen soil
    call phasechange(sltyp, NSNOW   ,NSOIL   ,ISNOW   ,DT      ,FACT    , & !in
         DZSNSO  ,HCPCT   ,IST     ,ILOC    ,JLOC    , & !in
         STC     ,SNICE   ,SNLIQ   ,sneqv   ,snowh   , & !inout
         SMC     ,soilwat    ,                            & !inout
         QMELT   ,IMELT   ,PONDING )                     !out
  end subroutine energy


  subroutine thermoprop(sltyp, NSOIL   ,NSNOW   ,ISNOW   ,IST     ,DZSNSO  , & !in
       DT      ,snowh   ,SNICE   ,SNLIQ   ,CSOIL   , & !in
       SMC     ,soilwat    ,TG      ,STC     ,UR      , & !in
       LAT     ,Z0M     ,ZLVL    ,lutyp, & !in
       DF      ,HCPCT   ,SNICEV  ,SNLIQV  ,epore   , & !out
       FACT    )                                       !out
    use noahmp_soil_param, only: LK_SMCMAX
    use noahmp_veg_param, only: ISURBAN
    implicit none
    ! inputs
    integer, intent(in) :: sltyp
    integer                        , intent(in)  :: NSOIL   !number of soil layers
    integer                        , intent(in)  :: NSNOW   !maximum no. of snow layers
    integer                        , intent(in)  :: ISNOW   !actual no. of snow layers
    integer                        , intent(in)  :: IST     !surface type
    real                           , intent(in)  :: DT      !time step [s]
    real, dimension(-NSNOW+1:    0), intent(in)  :: SNICE   !snow ice mass (kg/m2)
    real, dimension(-NSNOW+1:    0), intent(in)  :: SNLIQ   !snow liq mass (kg/m2)
    real, dimension(-NSNOW+1:NSOIL), intent(in)  :: DZSNSO  !thickness of snow/soil layers [m]
    real, dimension(       1:NSOIL), intent(in)  :: SMC     !soil moisture (ice + liq.) [m3/m3]
    real, dimension(       1:NSOIL), intent(in)  :: soilwat    !liquid soil moisture [m3/m3]
    real                           , intent(in)  :: snowh   !snow height [m]
    real                           , intent(in)  :: CSOIL   !vol. soil heat capacity [j/m3/k]
    real,                            intent(in)  :: TG      !surface temperature (k)
    real, dimension(-NSNOW+1:NSOIL), intent(in)  :: STC     !snow/soil/lake temp. (k)
    real,                            intent(in)  :: UR      !wind speed at ZLVL (m/s)
    real,                            intent(in)  :: LAT     !latitude (radians)
    real,                            intent(in)  :: Z0M     !roughness length (m)
    real,                            intent(in)  :: ZLVL    !reference height (m)
    integer                        , intent(in)  :: lutyp  !lutyp type

    ! outputs
    real, dimension(-NSNOW+1:NSOIL), intent(out) :: DF      !thermal conductivity [w/m/k]
    real, dimension(-NSNOW+1:NSOIL), intent(out) :: HCPCT   !heat capacity [j/m3/k]
    real, dimension(-NSNOW+1:    0), intent(out) :: SNICEV  !partial volume of ice [m3/m3]
    real, dimension(-NSNOW+1:    0), intent(out) :: SNLIQV  !partial volume of liquid water [m3/m3]
    real, dimension(-NSNOW+1:    0), intent(out) :: epore   !effective porosity [m3/m3]
    real, dimension(-NSNOW+1:NSOIL), intent(out) :: FACT    !computing energy for phase change

    ! locals
    integer :: IZ
    real, dimension(-NSNOW+1:    0)              :: CVSNO   !volumetric specific heat (j/m3/k)
    real, dimension(-NSNOW+1:    0)              :: TKSNO   !snow thermal conductivity (j/m3/k)
    real, dimension(       1:NSOIL)              :: soilice    !soil ice content

    ! compute snow thermal conductivity and heat capacity

    call CSNOW (ISNOW   ,NSNOW   ,NSOIL   ,SNICE   ,SNLIQ   ,DZSNSO  , & !in
         TKSNO   ,CVSNO   ,SNICEV  ,SNLIQV  ,epore   )   !out

    do IZ = ISNOW+1, 0
       DF   (IZ) = TKSNO(IZ)
       HCPCT(IZ) = CVSNO(IZ)
    end do

    ! compute soil thermal properties

    do  IZ = 1, NSOIL
       soilice(IZ)  = SMC(IZ) - soilwat(IZ)
       HCPCT(IZ) = soilwat(IZ) * CWAT + (1.0 - LK_SMCMAX(sltyp)) * CSOIL &
            + (LK_SMCMAX(sltyp) - SMC(IZ)) * CPAIR + soilice(IZ) * CICE
       call tdfcnd(sltyp, DF(IZ), SMC(IZ), soilwat(IZ))
    end do

    if (lutyp == ISURBAN) then
       do IZ = 1, NSOIL
          DF(IZ) = 3.24
       end do
    end if

    ! heat flux reduction effect from the overlying green canopy, adapted from
    ! section 2.1.2 of Peters-Lidard et al. (1997, JGR, VOL 102(D4)).
    ! not in use because of the separation of the canopy layer from the ground.
    ! but this may represent the effects of leaf litter (Niu comments)
    !       DF1 = DF1 * EXP (SBETA * SHDFAC)

    ! compute lake thermal properties
    ! (no consideration of turbulent mixing for this version)

    if (IST == 2) then
       do IZ = 1, NSOIL
          if (STC(IZ) > TFRZ) then
             HCPCT(IZ) = CWAT
             DF(IZ)    = TKWAT  !+ KEDDY * CWAT
          else
             HCPCT(IZ) = CICE
             DF(IZ)    = TKICE
          end if
       end do
    end if

    ! combine a temporary variable used for melting/freezing of snow and frozen soil

    do IZ = ISNOW+1, NSOIL
       FACT(IZ) = DT / (HCPCT(IZ) * DZSNSO(IZ))
    end do

    ! snow/soil interface

    if (ISNOW == 0) then
       DF(1) = (DF(1) * DZSNSO(1) + 0.35 * snowh) / (snowh + DZSNSO(1))
    else
       DF(1) = (DF(1) * DZSNSO(1) + DF(0) * DZSNSO(0)) / (DZSNSO(0) + DZSNSO(1))
    end if
  end subroutine thermoprop


  subroutine csnow(ISNOW   ,NSNOW   ,NSOIL   ,SNICE   ,SNLIQ   ,DZSNSO  , & !in
       TKSNO   ,CVSNO   ,SNICEV  ,SNLIQV  ,epore   )   !out
    ! Snow bulk density,volumetric capacity, and thermal conductivity
    implicit none

    ! inputs
    integer,                          intent(in) :: ISNOW  !number of snow layers (-)
    integer                        ,  intent(in) :: NSNOW  !maximum no. of snow layers
    integer                        ,  intent(in) :: NSOIL  !number of soil layers
    real, dimension(-NSNOW+1:    0),  intent(in) :: SNICE  !snow ice mass (kg/m2)
    real, dimension(-NSNOW+1:    0),  intent(in) :: SNLIQ  !snow liq mass (kg/m2)
    real, dimension(-NSNOW+1:NSOIL),  intent(in) :: DZSNSO !snow/soil layer thickness [m]

    ! outputs
    real, dimension(-NSNOW+1:    0), intent(out) :: CVSNO  !volumetric specific heat (j/m3/k)
    real, dimension(-NSNOW+1:    0), intent(out) :: TKSNO  !thermal conductivity (w/m/k)
    real, dimension(-NSNOW+1:    0), intent(out) :: SNICEV !partial volume of ice [m3/m3]
    real, dimension(-NSNOW+1:    0), intent(out) :: SNLIQV !partial volume of liquid water [m3/m3]
    real, dimension(-NSNOW+1:    0), intent(out) :: epore  !effective porosity [m3/m3]

    ! locals
    integer :: IZ
    real, dimension(-NSNOW+1:    0) :: BDSNOI  !bulk density of snow(kg/m3)

    !---------------------------------------------------------------------------------------------------
    ! thermal capacity of snow

    do IZ = ISNOW+1, 0
       SNICEV(IZ)   = min(1.0, SNICE(IZ) / (DZSNSO(IZ) * DENICE) )
       epore(IZ)    = 1.0 - SNICEV(IZ)
       SNLIQV(IZ)   = min(epore(IZ), SNLIQ(IZ) / (DZSNSO(IZ) * DENWAT))
    end do

    do IZ = ISNOW+1, 0
       BDSNOI(IZ) = (SNICE(IZ) + SNLIQ(IZ)) / DZSNSO(IZ)
       CVSNO(IZ) = CICE * SNICEV(IZ) + CWAT * SNLIQV(IZ)
       !      CVSNO(IZ) = 0.525E06                          ! constant
    end do

    ! thermal conductivity of snow

    do IZ = ISNOW+1, 0
       TKSNO(IZ) = 3.2217E-6 * BDSNOI(IZ) ** 2.0          ! Stieglitz(yen,1965)
       !    TKSNO(IZ) = 2E-2+2.5E-6*BDSNOI(IZ)*BDSNOI(IZ)   ! Anderson, 1976
       !    TKSNO(IZ) = 0.35                                ! constant
       !    TKSNO(IZ) = 2.576E-6*BDSNOI(IZ)**2. + 0.074    ! Verseghy (1991)
       !    TKSNO(IZ) = 2.22*(BDSNOI(IZ)/1000.)**1.88      ! Douvill(Yen, 1981)
    end do

  end subroutine csnow


  subroutine tdfcnd(sltyp, DF, SMC, soilwat)
    ! Calculate thermal diffusivity and conductivity of the soil.
    ! Peters-Lidard approach (Peters-Lidard et al., 1998)
    ! Code history:
    ! June 2001 changes: frozen soil condition.
    use noahmp_soil_param, only: LK_SMCMAX
    use noahmp_soil_param, only: LK_QUARTZ
    implicit none
    integer, intent(in) :: sltyp
    real, intent(in)    :: SMC    ! total soil water
    real, intent(in)    :: soilwat   ! liq. soil water
    real, intent(out)   :: DF     ! thermal diffusivity

    ! local variables
    real  :: AKE
    real  :: GAMMD
    real  :: THKDRY
    real  :: THKO     ! thermal conductivity for other soil components
    real  :: THKQTZ   ! thermal conductivity for quartz
    real  :: THKSAT   !
    real  :: THKS     ! thermal conductivity for the solids
    real  :: THKW     ! water thermal conductivity
    real  :: SATRATIO
    real  :: XU
    real  :: XUNFROZ
    ! We now get quartz as an input argument:
    !      DATA QUARTZ /0.82, 0.10, 0.25, 0.60, 0.52,
    !     &             0.35, 0.60, 0.40, 0.82/
    !
    ! If the soil has any moisture content compute a partial sum/product
    ! otherwise use a constant value which works well with most soils
    !
    !  QUARTZ ....QUARTZ CONTENT (SOIL TYPE DEPENDENT)
    !
    ! USE AS IN PETERS-LIDARD, 1998 (MODIF. FROM JOHANSEN, 1975).

    !                                  PABLO GRUNMANN, 08/17/98
    ! Refs.:
    !      Farouki, O.T.,1986: Thermal properties of soils. Series on Rock
    !              and Soil Mechanics, Vol. 11, Trans Tech, 136 pp.
    !      Johansen, O., 1975: Thermal conductivity of soils. PH.D. Thesis,
    !              University of Trondheim,
    !      Peters-Lidard, C. D., et al., 1998: The effect of soil thermal
    !              conductivity parameterization on surface energy fluxes
    !              and temperatures. Journal of The Atmospheric Sciences,
    !              Vol. 55, pp. 1209-1224.
    !
    ! NEEDS PARAMETERS
    ! POROSITY(SOIL TYPE):
    !      POROS = SMCMAX
    ! SATURATION RATIO:
    ! PARAMETERS  W/(M.K)
    SATRATIO = SMC / LK_SMCMAX(sltyp)
    THKW = 0.57
    !      IF (QUARTZ .LE. 0.2) THKO = 3.0
    THKO = 2.0
    ! SOLIDS' CONDUCTIVITY
    ! QUARTZ' CONDUCTIVITY
    THKQTZ = 7.7

    ! UNFROZEN FRACTION (FROM 1., i.e., 100%LIQUID, TO 0. (100% FROZEN))
    THKS = (THKQTZ ** LK_QUARTZ(sltyp))* (THKO ** (1.0 - LK_QUARTZ(sltyp)))

    ! UNFROZEN VOLUME FOR SATURATION (POROSITY*XUNFROZ)
    XUNFROZ = soilwat / SMC
    ! SATURATED THERMAL CONDUCTIVITY
    XU = XUNFROZ * LK_SMCMAX(sltyp)

    ! DRY DENSITY IN KG/M3
    THKSAT = THKS ** (1.0 - LK_SMCMAX(sltyp))* TKICE ** (LK_SMCMAX(sltyp) - XU)* THKW ** XU

    ! DRY THERMAL CONDUCTIVITY IN W.M-1.K-1
    GAMMD = (1.0 - LK_SMCMAX(sltyp)) * 2700.0

    THKDRY = (0.135 * GAMMD + 64.7)/ (2700.0 - 0.947 * GAMMD)
    ! FROZEN
    if ((soilwat + 0.0005) <  SMC) then
       AKE = SATRATIO
       ! UNFROZEN
       ! RANGE OF VALIDITY FOR THE KERSTEN NUMBER (AKE)
    else
       ! KERSTEN NUMBER (USING "FINE" FORMULA, VALID FOR SOILS CONTAINING AT
       ! LEAST 5% OF PARTICLES WITH DIAMETER LESS THAN 2.E-6 METERS.)
       ! (FOR "COARSE" FORMULA, SEE PETERS-LIDARD ET AL., 1998).
       if (SATRATIO > 0.1) then
          AKE = log10(SATRATIO) + 1.0
          ! USE K = KDRY
       else
          AKE = 0.0
       end if
       !  THERMAL CONDUCTIVITY
    end if

    DF = AKE * (THKSAT - THKDRY) + THKDRY

  end subroutine tdfcnd


  subroutine radiation(lutyp  ,IST     ,ISC     ,ICE     ,NSOIL   , & !in
       SNEQVO  ,sneqv   ,DT      ,COSZ    ,snowh   , & !in
       TG      ,TV      ,FSNO    ,QSNOW   ,FWET    , & !in
       elai    ,esai    ,SMC     ,SOLAD   ,SOLAI   , & !in
       fveg    ,ILOC    ,JLOC    ,                   & !in
       ALBOLD  ,TAUSS   ,                            & !inout
       FSUN    ,LAISUN  ,LAISHA  ,PARSUN  ,PARSHA  , & !out
       SAV     ,SAG     ,FSR     ,FSA     ,FSRV    , &
       FSRG    ,BGAP    ,WGAP)            !out
    implicit none
    ! input
    integer, intent(in)                  :: ILOC
    integer, intent(in)                  :: JLOC
    integer, intent(in)                  :: lutyp !vegetation type
    integer, intent(in)                  :: IST    !surface type
    integer, intent(in)                  :: ISC    !soil color type (1-lighest; 8-darkest)
    integer, intent(in)                  :: ICE    !ice (ice = 1)
    integer, intent(in)                  :: NSOIL  !number of soil layers

    real, intent(in)                     :: DT     !time step [s]
    real, intent(in)                     :: QSNOW  !snowfall (mm/s)
    real, intent(in)                     :: SNEQVO !snow mass at last time step(mm)
    real, intent(in)                     :: sneqv  !snow mass (mm)
    real, intent(in)                     :: snowh  !snow height (mm)
    real, intent(in)                     :: COSZ   !cosine solar zenith angle (0-1)
    real, intent(in)                     :: TG     !ground temperature (k)
    real, intent(in)                     :: TV     !vegetation temperature (k)
    real, intent(in)                     :: elai   !LAI, one-sided, adjusted for burying by snow
    real, intent(in)                     :: esai   !SAI, one-sided, adjusted for burying by snow
    real, intent(in)                     :: FWET   !fraction of canopy that is wet
    real, dimension(1:NSOIL), intent(in) :: SMC    !volumetric soil water [m3/m3]
    real, dimension(1:2)    , intent(in) :: SOLAD  !incoming direct solar radiation (w/m2)
    real, dimension(1:2)    , intent(in) :: SOLAI  !incoming diffuse solar radiation (w/m2)
    real, intent(in)                     :: FSNO   !snow cover fraction (-)
    real, intent(in)                     :: fveg   !green vegetation fraction [0.0-1.0]

    ! inout
    real,                  intent(inout) :: ALBOLD !snow albedo at last time step (CLASS type)
    real,                  intent(inout) :: TAUSS  !non-dimensional snow age.

    ! output
    real, intent(out)                    :: FSUN   !sunlit fraction of canopy (-)
    real, intent(out)                    :: LAISUN !sunlit leaf area (-)
    real, intent(out)                    :: LAISHA !shaded leaf area (-)
    real, intent(out)                    :: PARSUN !average absorbed par for sunlit leaves (w/m2)
    real, intent(out)                    :: PARSHA !average absorbed par for shaded leaves (w/m2)
    real, intent(out)                    :: SAV    !solar radiation absorbed by vegetation (w/m2)
    real, intent(out)                    :: SAG    !solar radiation absorbed by ground (w/m2)
    real, intent(out)                    :: FSA    !total absorbed solar radiation (w/m2)
    real, intent(out)                    :: FSR    !total reflected solar radiation (w/m2)

    !jref:start
    real, intent(out)                    :: FSRV    !veg. reflected solar radiation (w/m2)
    real, intent(out)                    :: FSRG    !ground reflected solar radiation (w/m2)
    real, intent(out)                    :: BGAP
    real, intent(out)                    :: WGAP
    !jref:end

    ! local
    real                                 :: FAGE   !snow age function (0 - new snow)
    real, dimension(1:2)                 :: ALBGRD !ground albedo (direct)
    real, dimension(1:2)                 :: ALBGRI !ground albedo (diffuse)
    real, dimension(1:2)                 :: ALBD   !surface albedo (direct)
    real, dimension(1:2)                 :: ALBI   !surface albedo (diffuse)
    real, dimension(1:2)                 :: FABD   !flux abs by veg (per unit direct flux)
    real, dimension(1:2)                 :: FABI   !flux abs by veg (per unit diffuse flux)
    real, dimension(1:2)                 :: FTDD   !down direct flux below veg (per unit dir flux)
    real, dimension(1:2)                 :: FTID   !down diffuse flux below veg (per unit dir flux)
    real, dimension(1:2)                 :: FTII   !down diffuse flux below veg (per unit dif flux)
    !jref:start
    real, dimension(1:2)                 :: FREVI
    real, dimension(1:2)                 :: FREVD
    real, dimension(1:2)                 :: FREGI
    real, dimension(1:2)                 :: FREGD
    !jref:end

    real                                 :: FSHA   !shaded fraction of canopy
    real                                 :: VAI    !total LAI + stem area index, one sided

    real, parameter :: MPE = 1.0E-6
    logical :: VEG  !true: vegetated for surface temperature calculation

    ! surface abeldo

    call ALBEDO (lutyp ,IST    ,ISC    ,ICE    ,NSOIL  , & !in
         DT     ,COSZ   ,FAGE   ,elai   ,esai   , & !in
         TG     ,TV     ,snowh  ,FSNO   ,FWET   , & !in
         SMC    ,SNEQVO ,sneqv  ,QSNOW  ,fveg   , & !in
         ILOC   ,JLOC   ,                         & !in
         ALBOLD ,TAUSS                          , & !inout
         ALBGRD ,ALBGRI ,ALBD   ,ALBI   ,FABD   , & !out
         FABI   ,FTDD   ,FTID   ,FTII   ,FSUN   , & !)   !out
         FREVI  ,FREVD   ,FREGD ,FREGI  ,BGAP   , & !inout
         WGAP)

    ! surface radiation

    FSHA = 1.0 - FSUN
    LAISUN = elai * FSUN
    LAISHA = elai * FSHA
    VAI = elai+ esai
    if (VAI > 0.0) then
       VEG = .true.
    else
       VEG = .false.
    end if

    call SURRAD(MPE    ,FSUN   ,FSHA   ,elai   ,VAI    , & !in
         LAISUN ,LAISHA ,SOLAD  ,SOLAI  ,FABD   , & !in
         FABI   ,FTDD   ,FTID   ,FTII   ,ALBGRD , & !in
         ALBGRI ,ALBD   ,ALBI   ,ILOC   ,JLOC   , & !in
         PARSUN ,PARSHA ,SAV    ,SAG    ,FSA    , & !out
         FSR    ,                                 & !out
         FREVI  ,FREVD  ,FREGD  ,FREGI  ,FSRV   , & !inout
         FSRG)

  end subroutine radiation


  subroutine albedo(lutyp ,IST    ,ISC    ,ICE    ,NSOIL  , & !in
       DT     ,COSZ   ,FAGE   ,elai   ,esai   , & !in
       TG     ,TV     ,snowh  ,FSNO   ,FWET   , & !in
       SMC    ,SNEQVO ,sneqv  ,QSNOW  ,fveg   , & !in
       ILOC   ,JLOC   ,                         & !in
       ALBOLD ,TAUSS                          , & !inout
       ALBGRD ,ALBGRI ,ALBD   ,ALBI   ,FABD   , & !out
       FABI   ,FTDD   ,FTID   ,FTII   ,FSUN   , & !out
       FREVI  ,FREVD  ,FREGD  ,FREGI  ,BGAP   , & !out
       WGAP)
    ! surface albedos. also fluxes (per unit incoming direct and diffuse
    ! radiation) reflected, transmitted, and absorbed by vegetation.
    ! also sunlit fraction of the canopy.
    use noahmp_global, only: nband
    use noahmp_veg_param, only: LK_RHOL
    use noahmp_veg_param, only: LK_RHOS
    use noahmp_veg_param, only: LK_TAUL
    use noahmp_veg_param, only: LK_TAUS
    implicit none
    ! input
    integer,                  intent(in)  :: ILOC
    integer,                  intent(in)  :: JLOC
    integer,                  intent(in)  :: NSOIL  !number of soil layers
    integer,                  intent(in)  :: lutyp !vegetation type
    integer,                  intent(in)  :: IST    !surface type
    integer,                  intent(in)  :: ISC    !soil color type (1-lighest; 8-darkest)
    integer,                  intent(in)  :: ICE    !ice (ice = 1)

    real,                     intent(in)  :: DT     !time step [sec]
    real,                     intent(in)  :: QSNOW  !snowfall
    real,                     intent(in)  :: COSZ   !cosine solar zenith angle for next time step
    real,                     intent(in)  :: snowh  !snow height (mm)
    real,                     intent(in)  :: TG     !ground temperature (k)
    real,                     intent(in)  :: TV     !vegetation temperature (k)
    real,                     intent(in)  :: elai   !LAI, one-sided, adjusted for burying by snow
    real,                     intent(in)  :: esai   !SAI, one-sided, adjusted for burying by snow
    real,                     intent(in)  :: FSNO   !fraction of grid covered by snow
    real,                     intent(in)  :: FWET   !fraction of canopy that is wet
    real,                     intent(in)  :: SNEQVO !snow mass at last time step(mm)
    real,                     intent(in)  :: sneqv  !snow mass (mm)
    real,                     intent(in)  :: fveg   !green vegetation fraction [0.0-1.0]
    real, dimension(1:NSOIL), intent(in)  :: SMC    !volumetric soil water (m3/m3)

    ! inout
    real,                  intent(inout)  :: ALBOLD !snow albedo at last time step (CLASS type)
    real,                  intent(inout)  :: TAUSS  !non-dimensional snow age

    ! output
    real, dimension(1:    2), intent(out) :: ALBGRD !ground albedo (direct)
    real, dimension(1:    2), intent(out) :: ALBGRI !ground albedo (diffuse)
    real, dimension(1:    2), intent(out) :: ALBD   !surface albedo (direct)
    real, dimension(1:    2), intent(out) :: ALBI   !surface albedo (diffuse)
    real, dimension(1:    2), intent(out) :: FABD   !flux abs by veg (per unit direct flux)
    real, dimension(1:    2), intent(out) :: FABI   !flux abs by veg (per unit diffuse flux)
    real, dimension(1:    2), intent(out) :: FTDD   !down direct flux below veg (per unit dir flux)
    real, dimension(1:    2), intent(out) :: FTID   !down diffuse flux below veg (per unit dir flux)
    real, dimension(1:    2), intent(out) :: FTII   !down diffuse flux below veg (per unit dif flux)
    real,                     intent(out) :: FSUN   !sunlit fraction of canopy (-)
    !jref:start
    real, dimension(1:    2), intent(out) :: FREVD
    real, dimension(1:    2), intent(out) :: FREVI
    real, dimension(1:    2), intent(out) :: FREGD
    real, dimension(1:    2), intent(out) :: FREGI
    real, intent(out) :: BGAP
    real, intent(out) :: WGAP
    !jref:end

    ! local
    real                 :: FAGE     !snow age function
    real                 :: ALB
    integer              :: IB       !indices
    integer              :: IC       !direct beam: ic=0; diffuse: ic=1

    real                 :: WL       !fraction of LAI+SAI that is LAI
    real                 :: WS       !fraction of LAI+SAI that is SAI
    real                 :: MPE      !prevents overflow for division by zero

    real, dimension(1:2) :: RHO      !leaf/stem reflectance weighted by fraction LAI and SAI
    real, dimension(1:2) :: TAU      !leaf/stem transmittance weighted by fraction LAI and SAI
    real, dimension(1:2) :: FTDI     !down direct flux below veg per unit dif flux = 0
    real, dimension(1:2) :: ALBSND   !snow albedo (direct)
    real, dimension(1:2) :: ALBSNI   !snow albedo (diffuse)

    real                 :: VAI      !elai+esai
    real                 :: GDIR     !average projected leaf/stem area in solar direction
    real                 :: EXT      !optical depth direct beam per unit leaf + stem area

    MPE = 1.0E-06
    BGAP = 0.0
    WGAP = 0.0

    ! initialize output because solar radiation only done if COSZ > 0

    do IB = 1, nband
       ALBD(IB) = 0.0
       ALBI(IB) = 0.0
       ALBGRD(IB) = 0.0
       ALBGRI(IB) = 0.0
       FABD(IB) = 0.0
       FABI(IB) = 0.0
       FTDD(IB) = 0.0
       FTID(IB) = 0.0
       FTII(IB) = 0.0
       if (IB == 1) FSUN = 0.0
    end do

    if (COSZ <= 0) return

    ! weight reflectance/transmittance by LAI and SAI

    do IB = 1, nband
       VAI = elai + esai
       WL  = elai / max(VAI,MPE)
       WS  = esai / max(VAI,MPE)
       RHO(IB) = max(LK_RHOL(ib,lutyp) * WL + LK_RHOS(ib,lutyp) * WS, MPE)
       TAU(IB) = max(LK_TAUL(ib,lutyp) * WL + LK_TAUS(ib,lutyp) * WS, MPE)
    end do

    ! snow age

    call snowage(DT, TG, SNEQVO, sneqv, TAUSS, FAGE)

    ! snow albedos: only if COSZ > 0 and FSNO > 0

    if (opt_alb == 1) &
         call SNOWALB_BATS(FSNO,COSZ,FAGE,ALBSND,ALBSNI)
    if (opt_alb == 2) then
       call SNOWALB_CLASS(QSNOW, DT, ALB, ALBOLD, ALBSND, ALBSNI, ILOC, JLOC)
       ALBOLD = ALB
    end if

    ! ground surface albedo

    call GROUNDALB(NSOIL, ICE     ,IST     ,ISC     , & !in
         FSNO    ,SMC     ,ALBSND  ,ALBSNI  ,COSZ    , & !in
         TG      ,ILOC    ,JLOC    ,                   & !in
         ALBGRD  ,ALBGRI  )                              !out

    ! loop over nband wavebands to calculate surface albedos and solar
    ! fluxes for unit incoming direct (IC=0) and diffuse flux (IC=1)

    do ib = 1, nband
       ic = 0      ! direct
       call TWOSTREAM(lutyp, ib, ic, cosz, vai, & !in
            & FWET   ,TV      ,ALBGRD  ,ALBGRI  ,RHO    , & !in
            & TAU    ,fveg    ,IST     ,ILOC    ,JLOC   , & !in
            & FABD   ,ALBD    ,FTDD    ,FTID    ,GDIR   , &!)   !out
            & FREVD  ,FREGD   ,BGAP    ,WGAP)

       ic = 1      ! diffuse
       call TWOSTREAM(lutyp, ib, ic, cosz, vai, & !in
            & FWET   ,TV      ,ALBGRD  ,ALBGRI  ,RHO    , & !in
            & TAU    ,fveg    ,IST     ,ILOC    ,JLOC   , & !in
            & FABI   ,ALBI    ,FTDI    ,FTII    ,GDIR   , & !)   !out
            & FREVI  ,FREGI   ,BGAP    ,WGAP)

    end do

    ! sunlit fraction of canopy. set FSUN = 0 if FSUN < 0.01.

    EXT = GDIR / COSZ * sqrt(1.0 - RHO(1) - TAU(1))
    FSUN = (1.0 - exp(-EXT * VAI)) / max(EXT * VAI, MPE)
    EXT = FSUN

    if (EXT < 0.01) then
       WL = 0.
    else
       WL = EXT
    end if
    FSUN = WL
  end subroutine albedo


  subroutine surrad(MPE     ,FSUN    ,FSHA    ,elai    ,VAI     , & !in
       LAISUN  ,LAISHA  ,SOLAD   ,SOLAI   ,FABD    , & !in
       FABI    ,FTDD    ,FTID    ,FTII    ,ALBGRD  , & !in
       ALBGRI  ,ALBD    ,ALBI    ,ILOC    ,JLOC    , & !in
       PARSUN  ,PARSHA  ,SAV     ,SAG     ,FSA     , & !out
       FSR     , & !)                                       !out
       FREVI   ,FREVD   ,FREGD   ,FREGI   ,FSRV    , &
       FSRG) !inout
    use noahmp_global, only: nband
    implicit none
    ! input

    integer, intent(in)              :: ILOC
    integer, intent(in)              :: JLOC
    real, intent(in)                 :: MPE     !prevents underflow errors if division by zero

    real, intent(in)                 :: FSUN    !sunlit fraction of canopy
    real, intent(in)                 :: FSHA    !shaded fraction of canopy
    real, intent(in)                 :: elai    !leaf area, one-sided
    real, intent(in)                 :: VAI     !leaf + stem area, one-sided
    real, intent(in)                 :: LAISUN  !sunlit leaf area index, one-sided
    real, intent(in)                 :: LAISHA  !shaded leaf area index, one-sided

    real, dimension(1:2), intent(in) :: SOLAD   !incoming direct solar radiation (w/m2)
    real, dimension(1:2), intent(in) :: SOLAI   !incoming diffuse solar radiation (w/m2)
    real, dimension(1:2), intent(in) :: FABD    !flux abs by veg (per unit incoming direct flux)
    real, dimension(1:2), intent(in) :: FABI    !flux abs by veg (per unit incoming diffuse flux)
    real, dimension(1:2), intent(in) :: FTDD    !down dir flux below veg (per incoming dir flux)
    real, dimension(1:2), intent(in) :: FTID    !down dif flux below veg (per incoming dir flux)
    real, dimension(1:2), intent(in) :: FTII    !down dif flux below veg (per incoming dif flux)
    real, dimension(1:2), intent(in) :: ALBGRD  !ground albedo (direct)
    real, dimension(1:2), intent(in) :: ALBGRI  !ground albedo (diffuse)
    real, dimension(1:2), intent(in) :: ALBD    !overall surface albedo (direct)
    real, dimension(1:2), intent(in) :: ALBI    !overall surface albedo (diffuse)

    real, dimension(1:2), intent(in) :: FREVD    !overall surface albedo veg (direct)
    real, dimension(1:2), intent(in) :: FREVI    !overall surface albedo veg (diffuse)
    real, dimension(1:2), intent(in) :: FREGD    !overall surface albedo grd (direct)
    real, dimension(1:2), intent(in) :: FREGI    !overall surface albedo grd (diffuse)

    ! output

    real, intent(out)                :: PARSUN  !average absorbed par for sunlit leaves (w/m2)
    real, intent(out)                :: PARSHA  !average absorbed par for shaded leaves (w/m2)
    real, intent(out)                :: SAV     !solar radiation absorbed by vegetation (w/m2)
    real, intent(out)                :: SAG     !solar radiation absorbed by ground (w/m2)
    real, intent(out)                :: FSA     !total absorbed solar radiation (w/m2)
    real, intent(out)                :: FSR     !total reflected solar radiation (w/m2)
    real, intent(out)                :: FSRV    !reflected solar radiation by vegetation
    real, intent(out)                :: FSRG    !reflected solar radiation by ground

    ! locals
    integer                          :: IB      !waveband number (1=vis, 2=nir)

    real                             :: ABS     !absorbed solar radiation (w/m2)
    real                             :: RNIR    !reflected solar radiation [nir] (w/m2)
    real                             :: RVIS    !reflected solar radiation [vis] (w/m2)
    real                             :: LAIFRA  !leaf area fraction of canopy
    real                             :: TRD     !transmitted solar radiation: direct (w/m2)
    real                             :: TRI     !transmitted solar radiation: diffuse (w/m2)
    real, dimension(1:2)             :: CAD     !direct beam absorbed by canopy (w/m2)
    real, dimension(1:2)             :: CAI     !diffuse radiation absorbed by canopy (w/m2)

    ! zero summed solar fluxes

    SAG = 0.0
    SAV = 0.0
    FSA = 0.0

    ! loop over nband wavebands

    do IB = 1, nband

       ! absorbed by canopy

       CAD(IB) = SOLAD(IB) * FABD(IB)
       CAI(IB) = SOLAI(IB) * FABI(IB)
       SAV     = SAV + CAD(IB) + CAI(IB)
       FSA     = FSA + CAD(IB) + CAI(IB)

       ! transmitted solar fluxes incident on ground

       TRD = SOLAD(IB) * FTDD(IB)
       TRI = SOLAD(IB) * FTID(IB) + SOLAI(IB) * FTII(IB)

       ! solar radiation absorbed by ground surface

       ABS = TRD * (1.0 - ALBGRD(IB)) + TRI * (1.0 - ALBGRI(IB))
       SAG = SAG + ABS
       FSA = FSA + ABS
    end do

    ! partition visible canopy absorption to sunlit and shaded fractions
    ! to get average absorbed par for sunlit and shaded leaves

    LAIFRA = elai / max(VAI, MPE)
    if (FSUN > 0.0) then
       PARSUN = (CAD(1) + FSUN * CAI(1)) * LAIFRA / max(LAISUN, MPE)
       PARSHA = (FSHA * CAI(1)) * LAIFRA / max(LAISHA, MPE)
    else
       PARSUN = 0.0
       PARSHA = (CAD(1) + CAI(1)) * LAIFRA /max(LAISHA, MPE)
    end if

    ! reflected solar radiation

    RVIS = ALBD(1) * SOLAD(1) + ALBI(1) * SOLAI(1)
    RNIR = ALBD(2) * SOLAD(2) + ALBI(2) * SOLAI(2)
    FSR  = RVIS + RNIR

    ! reflected solar radiation of veg. and ground (combined ground)
    FSRV = FREVD(1) * SOLAD(1) + FREVI(1) * SOLAI(1) + FREVD(2) * SOLAD(2) + FREVI(2) * SOLAI(2)
    FSRG = FREGD(1) * SOLAD(1) + FREGI(1) * SOLAI(1) + FREGD(2) * SOLAD(2) + FREGI(2) * SOLAI(2)


  end subroutine surrad


  subroutine snowage(DT, TG, SNEQVO, sneqv, TAUSS, FAGE)
    use noahmp_gen_param, only: KK_SWEMAX
    implicit none
    ! snowage from BATS
    !input
    real, intent(in) :: DT        !main time step (s)
    real, intent(in) :: TG        !ground temperature (k)
    real, intent(in) :: SNEQVO    !snow mass at last time step(mm)
    real, intent(in) :: sneqv     !snow water per unit ground area (mm)

    !output
    real, intent(out) :: FAGE     !snow age

    !input/output
    real, intent(inout) :: TAUSS      !non-dimensional snow age
    !local
    real            :: TAGE       !total aging effects
    real            :: AGE1       !effects of grain growth due to vapor diffusion
    real            :: AGE2       !effects of grain growth at freezing of melt water
    real            :: AGE3       !effects of soot
    real            :: DELA       !temporary variable
    real            :: SGE        !temporary variable
    real            :: DELS       !temporary variable
    real            :: DELA0      !temporary variable
    real            :: ARG        !temporary variable
    ! See Yang et al. (1997) J.of Climate for detail.

    if (sneqv <= 0.0) then
       TAUSS = 0.0
    else if (sneqv > 800.0) then
       TAUSS = 0.0
    else
       DELA0 = 1.0E-6 * DT
       ARG   = 5.0E3 * (1.0 / TFRZ - 1.0 / TG)
       AGE1  = exp(ARG)
       AGE2  = exp(AMIN1(0.0, 10.0 * ARG))
       AGE3  = 0.3
       TAGE  = AGE1 + AGE2 + AGE3
       DELA  = DELA0 * TAGE
       DELS  = AMAX1(0.0, sneqv - SNEQVO) / KK_SWEMAX
       SGE   = (TAUSS + DELA) * (1.0 - DELS)
       TAUSS = AMAX1(0.0, SGE)
    end if

    FAGE= TAUSS / (TAUSS + 1.0)

  end subroutine snowage


  subroutine snowalb_bats(FSNO, COSZ, FAGE, ALBSND, ALBSNI)
    use noahmp_global, only: nband
    implicit none
    ! input

    real,intent(in) :: COSZ    !cosine solar zenith angle
    real,intent(in) :: FSNO    !snow cover fraction (-)
    real,intent(in) :: FAGE    !snow age correction

    ! output

    real, dimension(1:2),intent(out) :: ALBSND !snow albedo for direct(1=vis, 2=nir)
    real, dimension(1:2),intent(out) :: ALBSNI !snow albedo for diffuse

    ! local
    integer :: IB          !waveband class

    real :: FZEN                 !zenith angle correction
    real :: CF1                  !temperary variable
    real :: SL2                  !2.*SL
    real :: SL1                  !1/SL
    real :: SL                   !adjustable parameter
    real, parameter :: C1 = 0.2  !default in BATS
    real, parameter :: C2 = 0.5  !default in BATS
    !  REAL, PARAMETER :: C1 = 0.2 * 2. ! double the default to match Sleepers River's
    !  REAL, PARAMETER :: C2 = 0.5 * 2. ! snow surface albedo (double aging effects)
    ! zero albedos for all points

    ALBSND(1:nband) = 0.0
    ALBSNI(1:nband) = 0.0

    ! when cosz > 0

    SL = 2.0
    SL1 = 1.0 / SL
    SL2 = 2.0 * SL
    CF1 = ((1.0 + SL1) / (1.0 + SL2 * COSZ) - SL1)
    FZEN = AMAX1(CF1, 0.0)

    ALBSNI(1) = 0.95 * (1.0 - C1 * FAGE)
    ALBSNI(2) = 0.65 * (1.0 - C2 * FAGE)

    ALBSND(1) = ALBSNI(1) + 0.4 * FZEN * (1.0 - ALBSNI(1))    !  vis direct
    ALBSND(2) = ALBSNI(2) + 0.4 * FZEN * (1.0 - ALBSNI(2))    !  nir direct

  end subroutine snowalb_bats


  subroutine snowalb_class(QSNOW, DT, ALB, ALBOLD, ALBSND, ALBSNI, &
       & ILOC, JLOC)
    use noahmp_global, only: nband
    use noahmp_gen_param, only: KK_SWEMAX
    implicit none

    ! input
    integer,intent(in) :: ILOC !grid index
    integer,intent(in) :: JLOC !grid index

    real,intent(in) :: QSNOW     !snowfall (mm/s)
    real,intent(in) :: DT        !time step (sec)
    real,intent(in) :: ALBOLD    !snow albedo at last time step

    ! in & out

    real,                intent(inout) :: ALB        !
    ! output

    real, dimension(1:2),intent(out) :: ALBSND !snow albedo for direct(1=vis, 2=nir)
    real, dimension(1:2),intent(out) :: ALBSNI !snow albedo for diffuse

    ! local
    integer :: IB          !waveband class

    ! zero albedos for all points

    ALBSND(1:nband) = 0.0
    ALBSNI(1:nband) = 0.0

    ! when cosz > 0

    ALB = 0.55 + (ALBOLD - 0.55) * exp(-0.01 * DT / 3600.0)

    ! 1 mm fresh snow(SWE) -- 10mm snow depth, assumed the fresh snow density 100kg/m3
    ! here assume 1cm snow depth will fully cover the old snow

    if (QSNOW > 0.0) then
       ALB = ALB + min(QSNOW * DT, KK_SWEMAX) * (0.84 - ALB) / KK_SWEMAX
    end if

    ALBSNI(1) = ALB         ! vis diffuse
    ALBSNI(2) = ALB         ! nir diffuse
    ALBSND(1) = ALB         ! vis direct
    ALBSND(2) = ALB         ! nir direct

  end subroutine snowalb_class


  subroutine groundalb(NSOIL, ICE, IST, ISC, FSNO, SMC, & !in
       & ALBSND, ALBSNI, COSZ, TG, ILOC, JLOC, &           !in
       & ALBGRD, ALBGRI)                                   !out
    use noahmp_global, only: nband
    use noahmp_gen_param, only: KK_ALBLAKE
    use noahmp_soil_param, only: LK_ALBSAT
    use noahmp_soil_param, only: LK_ALBDRY
    implicit none

    ! inputs
    integer,                  intent(in)  :: ILOC   !grid index
    integer,                  intent(in)  :: JLOC   !grid index
    integer,                  intent(in)  :: NSOIL  !number of soil layers
    integer,                  intent(in)  :: ICE    !value of ist for land ice
    integer,                  intent(in)  :: IST    !surface type
    integer,                  intent(in)  :: ISC    !soil color class (1-lighest; 8-darkest)
    real,                     intent(in)  :: FSNO   !fraction of surface covered with snow (-)
    real,                     intent(in)  :: TG     !ground temperature (k)
    real,                     intent(in)  :: COSZ   !cosine solar zenith angle (0-1)
    real, dimension(1:NSOIL), intent(in)  :: SMC    !volumetric soil water content (m3/m3)
    real, dimension(1:    2), intent(in)  :: ALBSND !direct beam snow albedo (vis, nir)
    real, dimension(1:    2), intent(in)  :: ALBSNI !diffuse snow albedo (vis, nir)

    !output

    real, dimension(1:    2), intent(out) :: ALBGRD !ground albedo (direct beam: vis, nir)
    real, dimension(1:    2), intent(out) :: ALBGRI !ground albedo (diffuse: vis, nir)

    !local
    integer                               :: IB     !waveband number (1=vis, 2=nir)
    real                                  :: INC    !soil water correction factor for soil albedo
    real                                  :: ALBSOD !soil albedo (direct)
    real                                  :: ALBSOI !soil albedo (diffuse)

    do IB = 1, nband
       INC = max(0.11 - 0.40 * SMC(1), 0.0)
       if (IST == 1)  then                     !soil
          ALBSOD = min(LK_ALBSAT(ib,isc) + INC, LK_ALBDRY(ib,isc))
          ALBSOI = ALBSOD
       else if (TG > TFRZ) then               !unfrozen lake, wetland
          ALBSOD = 0.06 / (max(0.01, COSZ) ** 1.7 + 0.15)
          ALBSOI = 0.06
       else                                      !frozen lake, wetland
          ALBSOD = KK_ALBLAKE(ib)
          ALBSOI = ALBSOD
       end if

       ! increase desert and semi-desert albedos

       if (IST == 1 .and. ISC == 9) then
          ALBSOD = ALBSOD + 0.10
          ALBSOI = ALBSOI + 0.10
       end if

       ALBGRD(IB) = ALBSOD*(1.0 - FSNO) + ALBSND(IB) * FSNO
       ALBGRI(IB) = ALBSOI*(1.0 - FSNO) + ALBSNI(IB) * FSNO
    end do

  end subroutine groundalb


  subroutine twostream(lutyp, ib, ic, COSZ    ,VAI    , & !in
       FWET   ,T       ,ALBGRD  ,ALBGRI  ,RHO    , & !in
       TAU    ,fveg    ,IST     ,ILOC    ,JLOC   , & !in
       FAB    ,FRE     ,FTD     ,FTI     ,GDIR   , & !)   !out
       FREV   ,FREG    ,BGAP    ,WGAP)
    ! use two-stream approximation of Dickinson (1983) Adv Geophysics
    ! 25:305-353 and Sellers (1985) Int J Remote Sensing 6:1335-1372
    ! to calculate fluxes absorbed by vegetation, reflected by vegetation,
    ! and transmitted through vegetation for unit incoming direct or diffuse
    ! flux given an underlying surface with known albedo.
    use noahmp_gen_param, only: KK_OMEGAS
    use noahmp_gen_param, only: KK_BETADS
    use noahmp_gen_param, only: KK_BETAIS
    use noahmp_veg_param, only: LK_XL
    use noahmp_veg_param, only: LK_HVT
    use noahmp_veg_param, only: LK_HVB
    use noahmp_veg_param, only: LK_RCROWN
    implicit none
    ! input

    integer,              intent(in)  :: ILOC    !grid index
    integer,              intent(in)  :: JLOC    !grid index
    integer,              intent(in)  :: IST     !surface type
    integer,              intent(in)  :: IB      !waveband number
    integer,              intent(in)  :: IC      !0=unit incoming direct; 1=unit incoming diffuse
    integer,              intent(in)  :: lutyp  !vegetation type

    real,                 intent(in)  :: COSZ    !cosine of direct zenith angle (0-1)
    real,                 intent(in)  :: VAI     !one-sided leaf+stem area index (m2/m2)
    real,                 intent(in)  :: FWET    !fraction of lai, sai that is wetted (-)
    real,                 intent(in)  :: T       !surface temperature (k)

    real, dimension(1:2), intent(in)  :: ALBGRD  !direct  albedo of underlying surface (-)
    real, dimension(1:2), intent(in)  :: ALBGRI  !diffuse albedo of underlying surface (-)
    real, dimension(1:2), intent(in)  :: RHO     !leaf+stem reflectance
    real, dimension(1:2), intent(in)  :: TAU     !leaf+stem transmittance
    real,                 intent(in)  :: fveg    !green vegetation fraction [0.0-1.0]

    ! output

    real, dimension(1:2), intent(out) :: FAB     !flux abs by veg layer (per unit incoming flux)
    real, dimension(1:2), intent(out) :: FRE     !flux refl above veg layer (per unit incoming flux)
    real, dimension(1:2), intent(out) :: FTD     !down dir flux below veg layer (per unit in flux)
    real, dimension(1:2), intent(out) :: FTI     !down dif flux below veg layer (per unit in flux)
    real,                 intent(out) :: GDIR    !projected leaf+stem area in solar direction
    real, dimension(1:2), intent(out) :: FREV    !flux reflected by veg layer   (per unit incoming flux)
    real, dimension(1:2), intent(out) :: FREG    !flux reflected by ground (per unit incoming flux)

    ! local
    real                              :: OMEGA   !fraction of intercepted radiation that is scattered
    real                              :: OMEGAL  !omega for leaves
    real                              :: BETAI   !upscatter parameter for diffuse radiation
    real                              :: BETAIL  !betai for leaves
    real                              :: BETAD   !upscatter parameter for direct beam radiation
    real                              :: BETADL  !betad for leaves
    real                              :: EXT     !optical depth of direct beam per unit leaf area
    real                              :: AVMU    !average diffuse optical depth

    real                              :: COSZI   !0.001 <= cosz <= 1.000
    real                              :: ASU     !single scattering albedo
    real                              :: CHIL    ! -0.4 <= xl <= 0.6

    real                              :: TMP0,TMP1,TMP2,TMP3,TMP4,TMP5,TMP6,TMP7,TMP8,TMP9
    real                              :: P1,P2,P3,P4,S1,S2,U1,U2,U3
    real                              :: B,C,D,D1,D2,F,H,H1,H2,H3,H4,H5,H6,H7,H8,H9,H10
    real                              :: PHI1,PHI2,SIGMA
    real                              :: FTDS,FTIS,FRES
    real                              :: DENfveg
    real                              :: VAI_SPREAD
    !jref:start
    real                              :: FREVEG,FREBAR,FTopt_veg,FTIVEG,FTDBAR,FTIBAR
    real                              :: THETAZ
    !jref:end

    !  variables for the modified two-stream scheme
    !  Niu and Yang (2004), JGR

    real, parameter :: PAI = 3.14159265
    real :: HD       !crown depth (m)
    real :: BB       !vertical crown radius (m)
    real :: THETAP   !angle conversion from SZA
    real :: FA       !foliage volume density (m-1)
    real :: NEWVAI   !effective LSAI (-)

    real,intent(inout) :: BGAP     !between canopy gap fraction for beam (-)
    real,intent(inout) :: WGAP     !within canopy gap fraction for beam (-)

    real :: KOPEN    !gap fraction for diffue light (-)
    real :: GAP      !total gap fraction for beam ( <=1-shafac )

    ! compute within and between gaps
    VAI_SPREAD = VAI
    if (VAI == 0.0) then
       GAP     = 1.0
       KOPEN   = 1.0
    else
       if (opt_rad == 1) then
          DENfveg = -log(max(1.0 - fveg, 0.01)) / (PAI * LK_RCROWN(lutyp) ** 2)
          HD      = LK_HVT(lutyp) - LK_HVB(lutyp)
          BB      = 0.5 * HD
          THETAP  = atan(BB / LK_RCROWN(lutyp) * tan(acos(max(0.01,COSZ))) )
          ! BGAP    = EXP(-DEN(lutyp) * PAI * RC(lutyp)**2/COS(THETAP) )
          BGAP    = exp(-DENfveg * PAI * LK_RCROWN(lutyp) ** 2 / cos(THETAP) )
          FA      = VAI / (1.33 * PAI * LK_RCROWN(lutyp) ** 3.0 * (BB / LK_RCROWN(lutyp)) * DENfveg)
          NEWVAI  = HD * FA
          WGAP    = (1.0 - BGAP) * exp(-0.5 * NEWVAI / COSZ)
          GAP     = min(1.0 - fveg, BGAP + WGAP)

          KOPEN   = 0.05
       end if

       if (opt_rad == 2) then
          GAP     = 0.0
          KOPEN   = 0.0
       end if

       if (opt_rad == 3) then
          GAP     = 1.0-fveg
          KOPEN   = 1.0-fveg
       end if
    end if

    ! calculate two-stream parameters OMEGA, BETAD, BETAI, AVMU, GDIR, EXT.
    ! OMEGA, BETAD, BETAI are adjusted for snow. values for OMEGA*BETAD
    ! and OMEGA*BETAI are calculated and then divided by the new OMEGA
    ! because the product OMEGA*BETAI, OMEGA*BETAD is used in solution.
    ! also, the transmittances and reflectances (TAU, RHO) are linear
    ! weights of leaf and stem values.

    COSZI  = max(0.001, COSZ)
    CHIL   = min(max(LK_XL(lutyp), -0.4), 0.6)
    if (abs(CHIL) <= 0.01) CHIL = 0.01
    PHI1   = 0.5 - 0.633*CHIL - 0.330*CHIL*CHIL
    PHI2   = 0.877 * (1.-2.*PHI1)
    GDIR   = PHI1 + PHI2*COSZI
    EXT    = GDIR/COSZI
    AVMU   = ( 1. - PHI1/PHI2 * log((PHI1+PHI2)/PHI1) ) / PHI2
    OMEGAL = RHO(IB) + TAU(IB)
    TMP0   = GDIR + PHI2*COSZI
    TMP1   = PHI1*COSZI
    ASU    = 0.5*OMEGAL*GDIR/TMP0 * ( 1.-TMP1/TMP0*log((TMP1+TMP0)/TMP1) )
    BETADL = (1.+AVMU*EXT)/(OMEGAL*AVMU*EXT)*ASU
    BETAIL = 0.5 * ( RHO(IB)+TAU(IB) + (RHO(IB)-TAU(IB))   &
         * ((1.+CHIL)/2.)**2 ) / OMEGAL

    ! adjust omega, betad, and betai for intercepted snow

    if (T > TFRZ) then                                !no snow
       TMP0 = OMEGAL
       TMP1 = BETADL
       TMP2 = BETAIL
    else
       TMP0 =  (1.0 - FWET) * OMEGAL          + FWET * KK_OMEGAS(ib)
       TMP1 = ((1.0 - FWET) * OMEGAL * BETADL + FWET * KK_OMEGAS(ib) * KK_BETADS ) / TMP0
       TMP2 = ((1.0 - FWET) * OMEGAL * BETAIL + FWET * KK_OMEGAS(ib) * KK_BETAIS ) / TMP0
    end if

    OMEGA = TMP0
    BETAD = TMP1
    BETAI = TMP2

    ! absorbed, reflected, transmitted fluxes per unit incoming radiation

    B = 1. - OMEGA + OMEGA*BETAI
    C = OMEGA*BETAI
    TMP0 = AVMU*EXT
    D = TMP0 * OMEGA*BETAD
    F = TMP0 * OMEGA*(1.-BETAD)
    TMP1 = B*B - C*C
    H = sqrt(TMP1) / AVMU
    SIGMA = TMP0*TMP0 - TMP1
    if ( abs (SIGMA) < 1.e-6 ) SIGMA = sign(1.e-6,SIGMA)
    P1 = B + AVMU*H
    P2 = B - AVMU*H
    P3 = B + TMP0
    P4 = B - TMP0
    S1 = exp(-H*VAI)
    S2 = exp(-EXT*VAI)
    if (IC == 0) then
       U1 = B - C/ALBGRD(IB)
       U2 = B - C*ALBGRD(IB)
       U3 = F + C*ALBGRD(IB)
    else
       U1 = B - C/ALBGRI(IB)
       U2 = B - C*ALBGRI(IB)
       U3 = F + C*ALBGRI(IB)
    end if
    TMP2 = U1 - AVMU*H
    TMP3 = U1 + AVMU*H
    D1 = P1*TMP2/S1 - P2*TMP3*S1
    TMP4 = U2 + AVMU*H
    TMP5 = U2 - AVMU*H
    D2 = TMP4/S1 - TMP5*S1
    H1 = -D*P4 - C*F
    TMP6 = D - H1*P3/SIGMA
    TMP7 = ( D - C - H1/SIGMA*(U1+TMP0) ) * S2
    H2 = ( TMP6*TMP2/S1 - P2*TMP7 ) / D1
    H3 = - ( TMP6*TMP3*S1 - P1*TMP7 ) / D1
    H4 = -F*P3 - C*D
    TMP8 = H4/SIGMA
    TMP9 = ( U3 - TMP8*(U2-TMP0) ) * S2
    H5 = - ( TMP8*TMP4/S1 + TMP9 ) / D2
    H6 = ( TMP8*TMP5*S1 + TMP9 ) / D2
    H7 = (C*TMP2) / (D1*S1)
    H8 = (-C*TMP3*S1) / D1
    H9 = TMP4 / (D2*S1)
    H10 = (-TMP5*S1) / D2

    ! downward direct and diffuse fluxes below vegetation
    ! Niu and Yang (2004), JGR.

    if (IC == 0) then
       FTDS = S2                           *(1.0-GAP) + GAP
       FTIS = (H4*S2/SIGMA + H5*S1 + H6/S1)*(1.0-GAP)
    else
       FTDS = 0.
       FTIS = (H9*S1 + H10/S1)*(1.0-KOPEN) + KOPEN
    end if
    FTD(IB) = FTDS
    FTI(IB) = FTIS

    ! flux reflected by the surface (veg. and ground)

    if (IC == 0) then
       FRES   = (H1/SIGMA + H2 + H3)*(1.0-GAP  ) + ALBGRD(IB)*GAP
       FREVEG = (H1/SIGMA + H2 + H3)*(1.0-GAP  )
       FREBAR = ALBGRD(IB)*GAP                   !jref - separate veg. and ground reflection
    else
       FRES   = (H7 + H8) *(1.0-KOPEN) + ALBGRI(IB)*KOPEN
       FREVEG = (H7 + H8) *(1.0-KOPEN) + ALBGRI(IB)*KOPEN
       FREBAR = 0                                !jref - separate veg. and ground reflection
    end if
    FRE(IB) = FRES

    FREV(IB) = FREVEG
    FREG(IB) = FREBAR

    ! flux absorbed by vegetation

    FAB(IB) = 1.0 - FRE(IB) - (1.0 - ALBGRD(IB)) * FTD(IB) &
         - (1.0 - ALBGRI(IB)) * FTI(IB)

    !if (iloc == 1 .and. jloc ==  2) then
    !  write(*,'(a7,2i2,5(a6,f8.4),2(a9,f8.4))') "ib,ic: ",ib,ic," GAP: ",GAP," FTD: ",FTD(IB)," FTI: ",FTI(IB)," FRE: ", &
    !         FRE(IB)," FAB: ",FAB(IB)," ALBGRD: ",ALBGRD(IB)," ALBGRI: ",ALBGRI(IB)
    !end if

  end subroutine twostream


  subroutine vege_flux(NSNOW   ,NSOIL   ,ISNOW   ,lutyp  ,VEG     , & !in
       DT      ,SAV     ,SAG     ,LWDN    ,UR      , & !in
       UU      ,VV      ,SFCTMP  ,THAIR   ,QAIR    , & !in
       EAIR    ,RHOAIR  ,snowh   ,VAI     ,GAMMAV   ,GAMMAG,  & !in
       FWET    ,LAISUN  ,LAISHA  ,CWP     ,DZSNSO  , & !in
       HTOP    ,ZLVL    ,ZPD     ,Z0M     ,fveg    , & !in
       Z0MG    ,EMV     ,EMG     ,CANLIQ  ,          & !in
       CANICE  ,STC     ,DF      ,RSSUN   ,RSSHA   , & !in
       RSURF   ,LATHEAV ,LATHEAG  ,PARSUN  ,PARSHA  ,IGS     , & !in
       FOLN    ,CO2AIR  ,O2AIR   ,BTRAN   ,SFCPRS  , & !in
       RHSUR   ,ILOC    ,JLOC    ,Q2      , & !in
       EAH     ,TAH     ,TV      ,TG      ,CM      , & !inout
       CH      ,DX      ,DZ8W    ,                   & !
       TAUXV   ,TAUYV   ,IRG     ,IRC     ,SHG     , & !out
       SHC     ,EVG     ,EVC     ,TR      ,GH      , & !out
       T2MV    ,PSNSUN  ,PSNSHA  ,                   & !out
       QC      ,PBLH    ,QSFC    ,PSFC, & !in
       IZ0TLND ,Q2V     ,CAH2,CHLEAF,CHUC)              !inout
    ! use newton-raphson iteration to solve for vegetation (tv) and
    ! ground (tg) temperatures that balance the surface energy budgets
    !
    ! vegetated:
    ! -SAV + IRC[TV] + SHC[TV] + EVC[TV] + TR[TV] = 0
    ! -SAG + IRG[TG] + SHG[TG] + EVG[TG] + GH[TG] = 0
    use noahmp_gen_param, only: KK_CZIL
    use noahmp_veg_param, only: ISURBAN
    implicit none
    ! input
    integer,                         intent(in) :: ILOC   !grid index
    integer,                         intent(in) :: JLOC   !grid index
    logical,                         intent(in) :: VEG    !true if vegetated surface
    integer,                         intent(in) :: NSNOW  !maximum no. of snow layers
    integer,                         intent(in) :: NSOIL  !number of soil layers
    integer,                         intent(in) :: ISNOW  !actual no. of snow layers
    integer,                         intent(in) :: lutyp !vegetation physiology type
    real,                            intent(in) :: fveg   !greeness vegetation fraction (-)
    real,                            intent(in) :: SAV    !solar rad absorbed by veg (w/m2)
    real,                            intent(in) :: SAG    !solar rad absorbed by ground (w/m2)
    real,                            intent(in) :: LWDN   !atmospheric longwave radiation (w/m2)
    real,                            intent(in) :: UR     !wind speed at height zlvl (m/s)
    real,                            intent(in) :: UU     !wind speed in eastward dir (m/s)
    real,                            intent(in) :: VV     !wind speed in northward dir (m/s)
    real,                            intent(in) :: SFCTMP !air temperature at reference height (k)
    real,                            intent(in) :: THAIR  !potential temp at reference height (k)
    real,                            intent(in) :: EAIR   !vapor pressure air at zlvl (pa)
    real,                            intent(in) :: QAIR   !specific humidity at zlvl (kg/kg)
    real,                            intent(in) :: RHOAIR !density air (kg/m**3)
    real,                            intent(in) :: DT     !time step (s)

    real,                            intent(in) :: snowh  !actual snow depth [m]
    real,                            intent(in) :: FWET   !wetted fraction of canopy
    real,                            intent(in) :: HTOP   !top of canopy layer (m)
    real,                            intent(in) :: CWP    !canopy wind parameter

    real,                            intent(in) :: VAI    !total leaf area index + stem area index
    real,                            intent(in) :: LAISUN !sunlit leaf area index, one-sided (m2/m2)
    real,                            intent(in) :: LAISHA !shaded leaf area index, one-sided (m2/m2)
    real,                            intent(in) :: ZLVL   !reference height (m)
    real,                            intent(in) :: ZPD    !zero plane displacement (m)
    real,                            intent(in) :: Z0M    !roughness length, momentum (m)
    real,                            intent(in) :: Z0MG   !roughness length, momentum, ground (m)
    real,                            intent(in) :: EMV    !vegetation emissivity
    real,                            intent(in) :: EMG    !ground emissivity

    real, dimension(-NSNOW+1:NSOIL), intent(in) :: STC    !soil/snow temperature (k)
    real, dimension(-NSNOW+1:NSOIL), intent(in) :: DF     !thermal conductivity of snow/soil (w/m/k)
    real, dimension(-NSNOW+1:NSOIL), intent(in) :: DZSNSO !thinkness of snow/soil layers (m)
    real,                            intent(in) :: CANLIQ !intercepted liquid water (mm)
    real,                            intent(in) :: CANICE !intercepted ice mass (mm)
    real,                            intent(in) :: RSURF  !ground surface resistance (s/m)
    !  REAL,                            INTENT(in) :: GAMMA  !psychrometric constant (pa/K)
    !  REAL,                            INTENT(in) :: LATHEA !latent heat of vaporization/subli (j/kg)
    real,                            intent(in) :: GAMMAV  !psychrometric constant (pa/K)
    real,                            intent(in) :: LATHEAV !latent heat of vaporization/subli (j/kg)
    real,                            intent(in) :: GAMMAG  !psychrometric constant (pa/K)
    real,                            intent(in) :: LATHEAG !latent heat of vaporization/subli (j/kg)
    real,                            intent(in) :: PARSUN !par absorbed per unit sunlit lai (w/m2)
    real,                            intent(in) :: PARSHA !par absorbed per unit shaded lai (w/m2)
    real,                            intent(in) :: FOLN   !foliage nitrogen (%)
    real,                            intent(in) :: CO2AIR !atmospheric co2 concentration (pa)
    real,                            intent(in) :: O2AIR  !atmospheric o2 concentration (pa)
    real,                            intent(in) :: IGS    !growing season index (0=off, 1=on)
    real,                            intent(in) :: SFCPRS !pressure (pa)
    real,                            intent(in) :: BTRAN  !soil water transpiration factor (0 to 1)
    real,                            intent(in) :: RHSUR  !raltive humidity in surface soil/snow air space (-)

    integer                        , intent(in) :: IZ0TLND
    real                           , intent(in) :: QC     !cloud water mixing ratio
    real                           , intent(in) :: PBLH   !planetary boundary layer height
    real                           , intent(in) :: PSFC   !pressure at lowest model layer
    real                           , intent(in) :: DX     !grid spacing
    real                           , intent(in) :: Q2     !mixing ratio (kg/kg)
    real                           , intent(in) :: DZ8W   !thickness of lowest layer
    real                           , intent(inout) :: QSFC   !mixing ratio at lowest model layer

    ! input/output
    real,                         intent(inout) :: EAH    !canopy air vapor pressure (pa)
    real,                         intent(inout) :: TAH    !canopy air temperature (k)
    real,                         intent(inout) :: TV     !vegetation temperature (k)
    real,                         intent(inout) :: TG     !ground temperature (k)
    real,                         intent(inout) :: CM     !momentum drag coefficient
    real,                         intent(inout) :: CH     !sensible heat exchange coefficient

    ! output
    ! -FSA + FIRA + FSH + (FCEV + FCTR + FGEV) + FCST + SSOIL = 0
    real,                           intent(out) :: TAUXV  !wind stress: e-w (n/m2)
    real,                           intent(out) :: TAUYV  !wind stress: n-s (n/m2)
    real,                           intent(out) :: IRC    !net longwave radiation (w/m2) [+= to atm]
    real,                           intent(out) :: SHC    !sensible heat flux (w/m2)     [+= to atm]
    real,                           intent(out) :: EVC    !evaporation heat flux (w/m2)  [+= to atm]
    real,                           intent(out) :: IRG    !net longwave radiation (w/m2) [+= to atm]
    real,                           intent(out) :: SHG    !sensible heat flux (w/m2)     [+= to atm]
    real,                           intent(out) :: EVG    !evaporation heat flux (w/m2)  [+= to atm]
    real,                           intent(out) :: TR     !transpiration heat flux (w/m2)[+= to atm]
    real,                           intent(out) :: GH     !ground heat (w/m2) [+ = to soil]
    real,                           intent(out) :: T2MV   !2 m height air temperature (k)
    real,                           intent(out) :: PSNSUN !sunlit leaf photosynthesis (umolco2/m2/s)
    real,                           intent(out) :: PSNSHA !shaded leaf photosynthesis (umolco2/m2/s)
    real,                           intent(out) :: CHLEAF !leaf exchange coefficient
    real,                           intent(out) :: CHUC   !under canopy exchange coefficient

    real,                           intent(out) :: Q2V
    real :: CAH    !sensible heat conductance, canopy air to ZLVL air (m/s)
    real :: U10V    !10 m wind speed in eastward dir (m/s)
    real :: V10V    !10 m wind speed in eastward dir (m/s)
    real :: WSPD

    ! locals
    real :: CW           !water vapor exchange coefficient
    real :: FV           !friction velocity (m/s)
    real :: WSTAR        !friction velocity n vertical direction (m/s) (only for SFCDIF2)
    real :: Z0H          !roughness length, sensible heat (m)
    real :: Z0HG         !roughness length, sensible heat (m)
    real :: RB           !bulk leaf boundary layer resistance (s/m)
    real :: RAMC         !aerodynamic resistance for momentum (s/m)
    real :: RAHC         !aerodynamic resistance for sensible heat (s/m)
    real :: RAWC         !aerodynamic resistance for water vapor (s/m)
    real :: RAMG         !aerodynamic resistance for momentum (s/m)
    real :: RAHG         !aerodynamic resistance for sensible heat (s/m)
    real :: RAWG         !aerodynamic resistance for water vapor (s/m)

    real, intent(out) :: RSSUN        !sunlit leaf stomatal resistance (s/m)
    real, intent(out) :: RSSHA        !shaded leaf stomatal resistance (s/m)

    real :: MOL          !Monin-Obukhov length (m)
    real :: DTV          !change in tv, last iteration (k)
    real :: DTG          !change in tg, last iteration (k)

    real :: AIR,CIR      !coefficients for ir as function of ts**4
    real :: CSH          !coefficients for sh as function of ts
    real :: CEV          !coefficients for ev as function of esat[ts]
    real :: CGH          !coefficients for st as function of ts
    real :: ATR,CTR      !coefficients for tr as function of esat[ts]
    real :: ATA,BTA      !coefficients for tah as function of ts
    real :: AEA,BEA      !coefficients for eah as function of esat[ts]

    real :: ESTV         !saturation vapor pressure at tv (pa)
    real :: ESTG         !saturation vapor pressure at tg (pa)
    real :: DESTV        !d(es)/dt at ts (pa/k)
    real :: DESTG        !d(es)/dt at tg (pa/k)
    real :: ESATW        !es for water
    real :: ESATI        !es for ice
    real :: DSATW        !d(es)/dt at tg (pa/k) for water
    real :: DSATI        !d(es)/dt at tg (pa/k) for ice

    real :: FM           !momentum stability correction, weighted by prior iters
    real :: FH           !sen heat stability correction, weighted by prior iters
    real :: FHG          !sen heat stability correction, ground
    real :: HCAN         !canopy height (m) [note: hcan >= z0mg]

    real :: A            !temporary calculation
    real :: B            !temporary calculation
    real :: CVH          !sensible heat conductance, leaf surface to canopy air (m/s)
    real :: CAW          !latent heat conductance, canopy air ZLVL air (m/s)
    real :: CTW          !transpiration conductance, leaf to canopy air (m/s)
    real :: CEW          !evaporation conductance, leaf to canopy air (m/s)
    real :: CGW          !latent heat conductance, ground to canopy air (m/s)
    real :: COND         !sum of conductances (s/m)
    real :: UC           !wind speed at top of canopy (m/s)
    real :: KH           !turbulent transfer coefficient, sensible heat, (m2/s)
    real :: H            !temporary sensible heat flux (w/m2)
    real :: HG           !temporary sensible heat flux (w/m2)
    real :: MOZ          !Monin-Obukhov stability parameter
    real :: MOZG         !Monin-Obukhov stability parameter
    real :: MOZOLD       !Monin-Obukhov stability parameter from prior iteration
    real :: FM2          !Monin-Obukhov momentum adjustment at 2m
    real :: FH2          !Monin-Obukhov heat adjustment at 2m
    real :: CH2          !Surface exchange at 2m
    real :: THSTAR          !Surface exchange at 2m

    real :: THVAIR
    real :: THAH
    real :: RAHC2        !aerodynamic resistance for sensible heat (s/m)
    real :: RAWC2        !aerodynamic resistance for water vapor (s/m)
    real, intent(out):: CAH2         !sensible heat conductance for diagnostics
    real :: CH2V         !exchange coefficient for 2m over vegetation.
    real :: CQ2V         !exchange coefficient for 2m over vegetation.
    real :: EAH2         !2m vapor pressure over canopy
    real :: QFX        !moisture flux
    real :: E1


    real :: VAIE         !total leaf area index + stem area index,effective
    real :: LAISUNE      !sunlit leaf area index, one-sided (m2/m2),effective
    real :: LAISHAE      !shaded leaf area index, one-sided (m2/m2),effective

    integer :: K         !index
    integer :: iter      !iteration index

    !jref - NITERC test from 5 to 20
    integer, parameter :: NITERC = 20   !number of iterations for surface temperature
    !jref - NITERG test from 3-5
    integer, parameter :: NITERG = 5   !number of iterations for ground temperature
    integer :: MOZSGN    !number of times MOZ changes sign
    real    :: MPE       !prevents overflow error if division by zero

    integer :: LITER     !Last iteration


    real :: T, TDC       !Kelvin to degree Celsius with limit -50 to +50

    character(len=80) ::  message

    TDC(T)   = min(50.0, max(-50.0, (T-TFRZ)) )

    MPE = 1E-6
    LITER = 0
    FV = 0.1

    ! initialization variables that do not depend on stability iteration
    DTV = 0.0
    DTG = 0.0
    MOZSGN = 0
    MOZOLD = 0.0
    HG = 0.0
    H = 0.0
    QFX = 0.0

    ! convert grid-cell LAI to the fractional vegetated area (fveg)
    VAIE    = min(6.0, VAI / fveg)
    LAISUNE = min(6.0, LAISUN / fveg)
    LAISHAE = min(6.0, LAISHA / fveg)

    ! saturation vapor pressure at ground temperature

    T = TDC(TG)
    call ESAT(T, ESATW, ESATI, DSATW, DSATI)
    if (T > 0.0) then
       ESTG = ESATW
    else
       ESTG = ESATI
    end if

    !jref - consistent surface specific humidity for sfcdif3 and sfcdif4

    QSFC = 0.622 * EAIR / (PSFC - 0.378 * EAIR)

    ! canopy height

    HCAN = HTOP
    UC = UR * log(HCAN / Z0M) / log(ZLVL / Z0M)
    if ((HCAN - ZPD) <= 0.0) then
       write(message,*) "CRITICAL PROBLEM: HCAN <= ZPD"
       call wrf_message ( message )
       write(message,*) 'i,j point=',ILOC, JLOC
       call wrf_message ( message )
       write(message,*) 'HCAN  =',HCAN
       call wrf_message ( message )
       write(message,*) 'ZPD   =',ZPD
       call wrf_message ( message )
       write (message, *) 'SNOWH =',snowh
       call wrf_message ( message )
       call wrf_error_fatal ( "CRITICAL PROBLEM IN MODULE_SF_NOAHMPLSM:VEGEFLUX" )
    end if

    ! prepare for longwave rad.
    AIR = -EMV*(1.0 + (1.0 - EMV) * (1.0 - EMG)) * LWDN - EMV * EMG * SB * TG ** 4
    CIR = (2.0 - EMV * (1.0 - EMG)) * EMV * SB

    loop1: do iter = 1, NITERC    !  begin stability iteration

       if (iter == 1) then
          Z0H  = Z0M
          Z0HG = Z0MG
       else
          Z0H  = Z0M    !* EXP(-CZIL*0.4*258.2*SQRT(FV*Z0M))
          Z0HG = Z0MG   !* EXP(-CZIL*0.4*258.2*SQRT(FV*Z0MG))
       end if

       ! aerodyn resistances between heights zlvl and d+z0v
       if (opt_sfc == 1) then
          call SFCDIF1(iter   ,SFCTMP ,RHOAIR ,H      ,QAIR   , & !in
               ZLVL   ,ZPD    ,Z0M    ,Z0H    ,UR     , & !in
               MPE    ,ILOC   ,JLOC   ,                 & !in
               MOZ    ,MOZSGN ,FM     ,FH     ,FM2,FH2, & !inout
               CM     ,CH     ,FV     ,CH2     )          !out
       end if

       if (opt_sfc == 2) then
          call SFCDIF2(iter   ,Z0M    ,TAH    ,THAIR  ,UR     , & !in
               KK_CZIL   ,ZLVL   ,ILOC   ,JLOC   ,         & !in
               CM     ,CH     ,MOZ    ,WSTAR  ,         & !in
               FV     )                                   !out
          ! Undo the multiplication by windspeed that SFCDIF2
          ! applies to exchange coefficients CH and CM:
          CH = CH / UR
          CM = CM / UR
       end if

       RAMC = max(1.0, 1.0 / (CM * UR))
       RAHC = max(1.0, 1.0 / (CH * UR))
       RAWC = RAHC

       ! aerodyn resistance between heights z0g and d+z0v, RAG, and leaf
       ! boundary layer resistance, RB
       call RAGRB(iter   ,VAIE   ,RHOAIR ,HG     ,TAH    , & !in
            ZPD    ,Z0MG   ,Z0HG   ,HCAN   ,UC     , & !in
            Z0H    ,FV     ,CWP    ,lutyp ,MPE    , & !in
            TV     ,MOZG   ,FHG    ,ILOC   ,JLOC   , & !inout
            RAMG   ,RAHG   ,RAWG   ,RB     )           !out

       ! es and d(es)/dt evaluated at tv
       T = TDC(TV)
       call ESAT(T, ESATW, ESATI, DSATW, DSATI)
       if (T > 0.0) then
          ESTV  = ESATW
          DESTV = DSATW
       else
          ESTV  = ESATI
          DESTV = DSATI
       end if

       ! stomatal resistance
       if (iter == 1) then
          if (opt_crs == 1) then  ! Ball-Berry
             call stomata(lutyp, igs, sfcprs, sfctmp, parsun, tv, &
                  &       eah, estv, o2air, co2air, foln, btran, &
                  &       rb, rssun, psnsun)
            
             call stomata(lutyp, igs, sfcprs, sfctmp, parsha, tv, &
                  &       eah, estv, o2air, co2air, foln, btran, &
                  &       rb, rssha, psnsha)
          end if

          if (opt_crs == 2) then  ! Jarvis
             call canres(lutyp, sfcprs, tv, parsun, eah, btran, rssun, psnsun)

             call canres(lutyp, sfcprs, tv, parsha, eah, btran, rssha, psnsha)
          end if
       end if

       ! prepare for sensible heat flux above veg.
       CAH  = 1.0 / RAHC
       CVH  = 2.0 * VAIE / RB
       CGH  = 1.0 / RAHG
       COND = CAH + CVH + CGH
       ATA  = (SFCTMP * CAH + TG * CGH) / COND
       BTA  = CVH / COND
       CSH  = (1.0 - BTA) * RHOAIR * CPAIR * CVH

       ! prepare for latent heat flux above veg.
       CAW  = 1.0 / RAWC
       CEW  = FWET * VAIE / RB
       CTW  = (1.0 - FWET) * (LAISUNE / (RB + RSSUN) + LAISHAE / (RB + RSSHA))
       CGW  = 1.0 / (RAWG + RSURF)
       COND = CAW + CEW + CTW + CGW
       AEA  = (EAIR * CAW + ESTG * CGW) / COND
       BEA  = (CEW + CTW) / COND
       CEV  = (1.0 - BEA) * CEW * RHOAIR * CPAIR / GAMMAV   ! Barlage: change to vegetation v3.6
       CTR  = (1.0 - BEA) * CTW * RHOAIR * CPAIR / GAMMAV

       ! evaluate surface fluxes with current temperature and solve for dts
       TAH = ATA + BTA * TV               ! canopy air T.
       EAH = AEA + BEA * ESTV             ! canopy air e

       IRC = fveg * (AIR + CIR * TV ** 4)
       SHC = fveg * RHOAIR * CPAIR * CVH * (  TV - TAH)
       EVC = fveg * RHOAIR * CPAIR * CEW * (ESTV - EAH) / GAMMAV ! Barlage: change to v in v3.6
       TR  = fveg * RHOAIR * CPAIR * CTW * (ESTV - EAH) / GAMMAV
       if (TV > TFRZ) then
          EVC = min(CANLIQ * LATHEAV / DT, EVC)    ! Barlage: add if block for canice in v3.6
       else
          EVC = min(CANICE * LATHEAV / DT, EVC)
       end if

       B   = SAV - IRC - SHC - EVC - TR                          !additional w/m2
       A   = fveg * (4.0 * CIR * TV ** 3 + CSH + (CEV + CTR) * DESTV) !volumetric heat capacity
       DTV = B / A

       IRC = IRC + fveg * 4.0 * CIR * TV ** 3 * DTV
       SHC = SHC + fveg * CSH * DTV
       EVC = EVC + fveg * CEV * DESTV * DTV
       TR  = TR  + fveg * CTR * DESTV * DTV

       ! update vegetation surface temperature
       TV  = TV + DTV
       !        TAH = ATA + BTA*TV               ! canopy air T; update here for consistency

       ! for computing M-O length in the next iteration
       H  = RHOAIR * CPAIR * (TAH - SFCTMP) /RAHC
       HG = RHOAIR * CPAIR * (TG - TAH) /RAHG

       ! consistent specific humidity from canopy air vapor pressure
       QSFC = (0.622 * EAH) / (SFCPRS - 0.378 * EAH)

       if (LITER == 1) then
          exit loop1
       end if
       if (iter >= 5 .and. abs(DTV) <= 0.01 .and. LITER == 0) then
          LITER = 1
       end if

    end do loop1 ! end stability iteration

    ! under-canopy fluxes and tg

    AIR = -EMG * (1.0 - EMV) * LWDN - EMG * EMV * SB * TV ** 4
    CIR = EMG * SB
    CSH = RHOAIR * CPAIR / RAHG
    CEV = RHOAIR * CPAIR / (GAMMAG * (RAWG + RSURF))  ! Barlage: change to ground v3.6
    CGH = 2.0 * DF(ISNOW+1) / DZSNSO(ISNOW+1)

    loop2: do iter = 1, NITERG

       T = TDC(TG)
       call ESAT(T, ESATW, ESATI, DSATW, DSATI)
       if (T > 0.0) then
          ESTG  = ESATW
          DESTG = DSATW
       else
          ESTG  = ESATI
          DESTG = DSATI
       end if

       IRG = CIR * TG ** 4 + AIR
       SHG = CSH * (TG - TAH)
       EVG = CEV * (ESTG * RHSUR - EAH)
       GH  = CGH * (TG - STC(ISNOW+1))

       B = SAG - IRG - SHG - EVG - GH
       A = 4.0 * CIR * TG ** 3 + CSH + CEV * DESTG + CGH
       DTG = B / A

       IRG = IRG + 4.0 * CIR * TG ** 3 * DTG
       SHG = SHG + CSH * DTG
       EVG = EVG + CEV * DESTG * DTG
       GH  = GH  + CGH * DTG
       TG  = TG  + DTG

    end do loop2

    !     TAH = (CAH*SFCTMP + CVH*TV + CGH*TG)/(CAH + CVH + CGH)

    ! if snow on ground and TG > TFRZ: reset TG = TFRZ. reevaluate ground fluxes.

    if (opt_stc == 1) then
       if (snowh > 0.05 .and. TG > TFRZ) then
          TG  = TFRZ
          IRG = CIR * TG ** 4 - EMG * (1.0 - EMV) * LWDN - EMG * EMV * SB * TV ** 4
          SHG = CSH * (TG         - TAH)
          EVG = CEV * (ESTG * RHSUR - EAH)
          GH  = SAG - (IRG + SHG + EVG)
       end if
    end if

    ! wind stresses

    TAUXV = -RHOAIR * CM * UR * UU
    TAUYV = -RHOAIR * CM * UR * VV

    ! consistent vegetation air temperature and vapor pressure since TG is not consistent with the TAH/EAH
    ! calculation.
    !     TAH = SFCTMP + (SHG+SHC)/(RHOAIR*CPAIR*CAH)
    !     TAH = SFCTMP + (SHG*fveg+SHC)/(RHOAIR*CPAIR*CAH) ! ground flux need fveg
    !     EAH = EAIR + (EVC+fveg*(TR+EVG))/(RHOAIR*CAW*CPAIR/GAMMAG )
    !     QFX = (QSFC-QAIR)*RHOAIR*CAW !*CPAIR/GAMMAG

    ! 2m temperature over vegetation ( corrected for low CQ2V values )
    if (opt_sfc == 1 .or. opt_sfc == 2) then
       !      CAH2 = FV*1./KARMAN*LOG((2.+Z0H)/Z0H)
       CAH2 = FV * KARMAN / log((2.0 + Z0H) / Z0H)
       CAH2 = FV * KARMAN / (log((2.0 + Z0H) / Z0H) - FH2)
       CQ2V = CAH2
       if (CAH2 < 1.E-5) then
          T2MV = TAH
          !         Q2V  = (EAH*0.622/(SFCPRS - 0.378*EAH))
          Q2V  = QSFC
       else
          T2MV = TAH - (SHG + SHC / fveg) / (RHOAIR * CPAIR) * 1.0 / CAH2
          !         Q2V = (EAH*0.622/(SFCPRS - 0.378*EAH))- QFX/(RHOAIR*FV)* 1./KARMAN * LOG((2.+Z0H)/Z0H)
          Q2V = QSFC - ((EVC + TR) / fveg + EVG) / (LATHEAV * RHOAIR) * 1.0 / CQ2V
       end if
    end if

    ! update CH for output
    CH = CAH
    CHLEAF = CVH
    CHUC = 1.0 / RAHG

  end subroutine vege_flux


  subroutine bare_flux(lutyp, dt, nsnow, nsoil, ISNOW, SAG     , & !in
       LWDN    ,UR      ,UU      ,VV      ,SFCTMP  , & !in
       THAIR   ,QAIR    ,EAIR    ,RHOAIR  ,snowh   , & !in
       DZSNSO  ,ZLVL    ,ZPD     ,Z0M     ,          & !in
       EMG     ,STC     ,DF      ,RSURF   ,LATHEA  , & !in
       GAMMA   ,RHSUR   ,ILOC    ,JLOC    ,Q2      , & !in
       TGB     ,CM      ,CH      ,          & !inout
       TAUXB   ,TAUYB   ,IRB     ,SHB     ,EVB     , & !out
       GHB     ,T2MB    ,DX      ,DZ8W, & !out
       QC      ,PBLH    ,QSFC    ,PSFC, & !in
       IZ0TLND ,SFCPRS  ,Q2B     ,EHB2)   !in
    ! use newton-raphson iteration to solve ground (tg) temperature
    ! that balances the surface energy budgets for bare soil fraction.
    !
    ! bare soil:
    ! -SAB + IRB[TG] + SHB[TG] + EVB[TG] + GHB[TG] = 0
    use noahmp_gen_param, only: KK_CZIL
    use noahmp_veg_param, only: ISURBAN
    use noahmp_veg_param, only: ISBARREN
    implicit none
    ! input
    integer                        , intent(in) :: lutyp != ISBARREN
    integer                        , intent(in) :: ILOC   !grid index
    integer                        , intent(in) :: JLOC   !grid index
    integer,                         intent(in) :: NSNOW  !maximum no. of snow layers
    integer,                         intent(in) :: NSOIL  !number of soil layers
    integer,                         intent(in) :: ISNOW  !actual no. of snow layers
    real,                            intent(in) :: DT     !time step (s)
    real,                            intent(in) :: SAG    !solar radiation absorbed by ground (w/m2)
    real,                            intent(in) :: LWDN   !atmospheric longwave radiation (w/m2)
    real,                            intent(in) :: UR     !wind speed at height zlvl (m/s)
    real,                            intent(in) :: UU     !wind speed in eastward dir (m/s)
    real,                            intent(in) :: VV     !wind speed in northward dir (m/s)
    real,                            intent(in) :: SFCTMP !air temperature at reference height (k)
    real,                            intent(in) :: THAIR  !potential temperature at height zlvl (k)
    real,                            intent(in) :: QAIR   !specific humidity at height zlvl (kg/kg)
    real,                            intent(in) :: EAIR   !vapor pressure air at height (pa)
    real,                            intent(in) :: RHOAIR !density air (kg/m3)
    real,                            intent(in) :: snowh  !actual snow depth [m]
    real, dimension(-NSNOW+1:NSOIL), intent(in) :: DZSNSO !thickness of snow/soil layers (m)
    real,                            intent(in) :: ZLVL   !reference height (m)
    real,                            intent(in) :: ZPD    !zero plane displacement (m)
    real,                            intent(in) :: Z0M    !roughness length, momentum, ground (m)
    real,                            intent(in) :: EMG    !ground emissivity
    real, dimension(-NSNOW+1:NSOIL), intent(in) :: STC    !soil/snow temperature (k)
    real, dimension(-NSNOW+1:NSOIL), intent(in) :: DF     !thermal conductivity of snow/soil (w/m/k)
    real,                            intent(in) :: RSURF  !ground surface resistance (s/m)
    real,                            intent(in) :: LATHEA !latent heat of vaporization/subli (j/kg)
    real,                            intent(in) :: GAMMA  !psychrometric constant (pa/k)
    real,                            intent(in) :: RHSUR  !raltive humidity in surface soil/snow air space (-)

    !jref:start; in
    integer                        , intent(in) :: IZ0TLND
    real                           , intent(in) :: QC     !cloud water mixing ratio
    real                           , intent(in) :: PBLH   !planetary boundary layer height
    real                           , intent(inout) :: QSFC   !mixing ratio at lowest model layer
    real                           , intent(in) :: PSFC   !pressure at lowest model layer
    real                           , intent(in) :: SFCPRS !pressure at lowest model layer
    real                           , intent(in) :: DX     !horisontal grid spacing
    real                           , intent(in) :: Q2     !mixing ratio (kg/kg)
    real                           , intent(in) :: DZ8W   !thickness of lowest layer
    !jref:end

    ! input/output
    real,                         intent(inout) :: TGB    !ground temperature (k)
    real,                         intent(inout) :: CM     !momentum drag coefficient
    real,                         intent(inout) :: CH     !sensible heat exchange coefficient

    ! output
    ! -SAB + IRB[TG] + SHB[TG] + EVB[TG] + GHB[TG] = 0

    real,                           intent(out) :: TAUXB  !wind stress: e-w (n/m2)
    real,                           intent(out) :: TAUYB  !wind stress: n-s (n/m2)
    real,                           intent(out) :: IRB    !net longwave rad (w/m2)   [+ to atm]
    real,                           intent(out) :: SHB    !sensible heat flux (w/m2) [+ to atm]
    real,                           intent(out) :: EVB    !latent heat flux (w/m2)   [+ to atm]
    real,                           intent(out) :: GHB    !ground heat flux (w/m2)  [+ to soil]
    real,                           intent(out) :: T2MB   !2 m height air temperature (k)
    !jref:start
    real,                           intent(out) :: Q2B    !bare ground heat conductance
    real :: EHB    !bare ground heat conductance
    real :: U10B    !10 m wind speed in eastward dir (m/s)
    real :: V10B    !10 m wind speed in eastward dir (m/s)
    real :: WSPD
    !jref:end

    ! local variables

    real :: TAUX       !wind stress: e-w (n/m2)
    real :: TAUY       !wind stress: n-s (n/m2)
    real :: FIRA       !total net longwave rad (w/m2)      [+ to atm]
    real :: FSH        !total sensible heat flux (w/m2)    [+ to atm]
    real :: FGEV       !ground evaporation heat flux (w/m2)[+ to atm]
    real :: SSOIL      !soil heat flux (w/m2)             [+ to soil]
    real :: FIRE       !emitted ir (w/m2)
    real :: TRAD       !radiative temperature (k)
    real :: TAH        !"surface" temperature at height z0h+zpd (k)

    real :: CW         !water vapor exchange coefficient
    real :: FV         !friction velocity (m/s)
    real :: WSTAR      !friction velocity n vertical direction (m/s) (only for SFCDIF2)
    real :: Z0H        !roughness length, sensible heat, ground (m)
    real :: RB         !bulk leaf boundary layer resistance (s/m)
    real :: RAMB       !aerodynamic resistance for momentum (s/m)
    real :: RAHB       !aerodynamic resistance for sensible heat (s/m)
    real :: RAWB       !aerodynamic resistance for water vapor (s/m)
    real :: MOL        !Monin-Obukhov length (m)
    real :: DTG        !change in tg, last iteration (k)

    real :: CIR        !coefficients for ir as function of ts**4
    real :: CSH        !coefficients for sh as function of ts
    real :: CEV        !coefficients for ev as function of esat[ts]
    real :: CGH        !coefficients for st as function of ts

    !jref:start
    real :: RAHB2      !aerodynamic resistance for sensible heat 2m (s/m)
    real :: RAWB2      !aerodynamic resistance for water vapor 2m (s/m)
    real, intent(out) :: EHB2      !sensible heat conductance for diagnostics
    real :: CH2B       !exchange coefficient for 2m temp.
    real :: CQ2B       !exchange coefficient for 2m temp.
    real :: THVAIR     !virtual potential air temp
    real :: THGH       !potential ground temp
    real :: EMB        !momentum conductance
    real :: QFX        !moisture flux
    real :: ESTG2      !saturation vapor pressure at 2m (pa)
    real :: E1
    !jref:end

    real :: ESTG       !saturation vapor pressure at tg (pa)
    real :: DESTG      !d(es)/dt at tg (pa/K)
    real :: ESATW      !es for water
    real :: ESATI      !es for ice
    real :: DSATW      !d(es)/dt at tg (pa/K) for water
    real :: DSATI      !d(es)/dt at tg (pa/K) for ice

    real :: A          !temporary calculation
    real :: B          !temporary calculation
    real :: H          !temporary sensible heat flux (w/m2)
    real :: MOZ        !Monin-Obukhov stability parameter
    real :: MOZOLD     !Monin-Obukhov stability parameter from prior iteration
    real :: FM         !momentum stability correction, weighted by prior iters
    real :: FH         !sen heat stability correction, weighted by prior iters
    integer :: MOZSGN  !number of times MOZ changes sign
    real :: FM2          !Monin-Obukhov momentum adjustment at 2m
    real :: FH2          !Monin-Obukhov heat adjustment at 2m
    real :: CH2          !Surface exchange at 2m

    integer :: iter    !iteration index
    integer, save :: NITERB = 5  !number of iterations for surface temperature
    real    :: MPE     !prevents overflow error if division by zero
    !jref:start

    real :: T, TDC     !Kelvin to degree Celsius with limit -50 to +50
    TDC(T)   = min( 50., max(-50.,(T-TFRZ)) )

    ! initialization variables that do not depend on stability iteration
    MPE = 1.0E-6
    DTG = 0.0
    MOZSGN = 0
    MOZOLD = 0.0
    H      = 0.0
    QFX    = 0.0
    FV     = 0.1

    CIR = EMG * SB
    CGH = 2.0 * DF(ISNOW+1) / DZSNSO(ISNOW+1)

    loop3: do iter = 1, NITERB  ! begin stability iteration

       if (iter == 1) then
          Z0H = Z0M
       else
          Z0H = Z0M !* EXP(-CZIL*0.4*258.2*SQRT(FV*Z0M))
       end if

       if (opt_sfc == 1) then
          call SFCDIF1(iter   ,SFCTMP ,RHOAIR ,H      ,QAIR   , & !in
               ZLVL   ,ZPD    ,Z0M    ,Z0H    ,UR     , & !in
               MPE    ,ILOC   ,JLOC   ,                 & !in
               MOZ    ,MOZSGN ,FM     ,FH     ,FM2,FH2, & !inout
               CM     ,CH     ,FV     ,CH2     )          !out
       end if

       if (opt_sfc == 2) then
          call SFCDIF2(iter   ,Z0M    ,TGB    ,THAIR  ,UR     , & !in
               KK_CZIL   ,ZLVL   ,ILOC   ,JLOC   ,         & !in
               CM     ,CH     ,MOZ    ,WSTAR  ,         & !in
               FV     )                                   !out
          ! Undo the multiplication by windspeed that SFCDIF2
          ! applies to exchange coefficients CH and CM:
          CH = CH / UR
          CM = CM / UR
          if (snowh > 0.0) then
             CM = min(0.01, CM)   ! CM & CH are too large, causing
             CH = min(0.01, CH)   ! computational instability
          end if

       end if

       RAMB = max(1.0, 1.0 / (CM * UR))
       RAHB = max(1.0, 1.0 / (CH * UR))
       RAWB = RAHB

       !jref - variables for diagnostics
       EMB = 1.0 / RAMB
       EHB = 1.0 / RAHB

       ! es and d(es)/dt evaluated at tg

       T = TDC(TGB)
       call ESAT(T, ESATW, ESATI, DSATW, DSATI)
       if (T > 0.0) then
          ESTG  = ESATW
          DESTG = DSATW
       else
          ESTG  = ESATI
          DESTG = DSATI
       end if

       CSH = RHOAIR * CPAIR / RAHB
       CEV = RHOAIR * CPAIR / GAMMA / (RSURF + RAWB)

       ! surface fluxes and dtg

       IRB   = CIR * TGB ** 4 - EMG * LWDN
       SHB   = CSH * (TGB        - SFCTMP)
       EVB   = CEV * (ESTG * RHSUR - EAIR)
       GHB   = CGH * (TGB        - STC(ISNOW+1))

       B     = SAG - IRB - SHB - EVB - GHB
       A     = 4.0 * CIR * TGB ** 3 + CSH + CEV * DESTG + CGH
       DTG   = B / A

       IRB = IRB + 4.0 * CIR * TGB ** 3 * DTG
       SHB = SHB + CSH * DTG
       EVB = EVB + CEV * DESTG * DTG
       GHB = GHB + CGH * DTG

       ! update ground surface temperature
       TGB = TGB + DTG

       ! for M-O length
       H = CSH * (TGB - SFCTMP)

       T = TDC(TGB)
       call ESAT(T, ESATW, ESATI, DSATW, DSATI)
       if (T > 0.0) then
          ESTG  = ESATW
       else
          ESTG  = ESATI
       end if
       QSFC = 0.622 * (ESTG * RHSUR) / (PSFC - 0.378 * (ESTG * RHSUR))

       QFX = (QSFC - QAIR) * CEV * GAMMA / CPAIR

    end do loop3 ! end stability iteration

    ! if snow on ground and TG > TFRZ: reset TG = TFRZ. reevaluate ground fluxes.
    if (opt_stc == 1) then
       if (snowh > 0.05 .and. TGB > TFRZ) then
          TGB = TFRZ
          IRB = CIR * TGB ** 4 - EMG * LWDN
          SHB = CSH * (TGB - SFCTMP)
          EVB = CEV * (ESTG * RHSUR - EAIR )          !ESTG reevaluate ?
          GHB = SAG - (IRB + SHB + EVB)
       end if
    end if

    ! wind stresses
    TAUXB = -RHOAIR * CM * UR * UU
    TAUYB = -RHOAIR * CM * UR * VV

    !jref:start; errors in original equation corrected.
    ! 2m air temperature
    if (opt_sfc == 1 .or. opt_sfc ==2) then
       EHB2  = FV * KARMAN / log((2.0 + Z0H) / Z0H)
       EHB2  = FV * KARMAN / (log((2.0 + Z0H) / Z0H) - FH2)
       CQ2B  = EHB2
       if (EHB2 < 1.0E-5 ) then
          T2MB  = TGB
          Q2B   = QSFC
       else
          T2MB  = TGB - SHB / (RHOAIR * CPAIR) * 1.0 / EHB2
          Q2B   = QSFC - EVB / (LATHEA * RHOAIR) * (1.0 / CQ2B + RSURF)
       end if
       if (lutyp == ISURBAN) Q2B = QSFC
    end if

    ! update CH
    CH = EHB
  end subroutine bare_flux


  subroutine ragrb(iter   ,VAI    ,RHOAIR ,HG     ,TAH    , & !in
       ZPD    ,Z0MG   ,Z0HG   ,HCAN   ,UC     , & !in
       Z0H    ,FV     ,CWP    ,lutyp ,MPE    , & !in
       TV     ,MOZG   ,FHG    ,ILOC   ,JLOC   , & !inout
       RAMG   ,RAHG   ,RAWG   ,RB     )           !out
    ! compute under-canopy aerodynamic resistance RAG and leaf boundary layer
    ! resistance RB
    use noahmp_veg_param, only: LK_DLEAF
    implicit none
    ! inputs
    integer,              intent(in) :: ILOC   !grid index
    integer,              intent(in) :: JLOC   !grid index
    integer,              intent(in) :: iter   !iteration index
    integer,              intent(in) :: lutyp !vegetation physiology type
    real,                 intent(in) :: VAI    !total LAI + stem area index, one sided
    real,                 intent(in) :: RHOAIR !density air (kg/m3)
    real,                 intent(in) :: HG     !ground sensible heat flux (w/m2)
    real,                 intent(in) :: TV     !vegetation temperature (k)
    real,                 intent(in) :: TAH    !air temperature at height z0h+zpd (k)
    real,                 intent(in) :: ZPD    !zero plane displacement (m)
    real,                 intent(in) :: Z0MG   !roughness length, momentum, ground (m)
    real,                 intent(in) :: HCAN   !canopy height (m) [note: hcan >= z0mg]
    real,                 intent(in) :: UC     !wind speed at top of canopy (m/s)
    real,                 intent(in) :: Z0H    !roughness length, sensible heat (m)
    real,                 intent(in) :: Z0HG   !roughness length, sensible heat, ground (m)
    real,                 intent(in) :: FV     !friction velocity (m/s)
    real,                 intent(in) :: CWP    !canopy wind parameter
    real,                 intent(in) :: MPE    !prevents overflow error if division by zero

    ! in & out
    real,              intent(inout) :: MOZG   !Monin-Obukhov stability parameter
    real,              intent(inout) :: FHG    !stability correction

    ! outputs
    real                             :: RAMG   !aerodynamic resistance for momentum (s/m)
    real                             :: RAHG   !aerodynamic resistance for sensible heat (s/m)
    real                             :: RAWG   !aerodynamic resistance for water vapor (s/m)
    real                             :: RB     !bulk leaf boundary layer resistance (s/m)


    real :: KH           !turbulent transfer coefficient, sensible heat, (m2/s)
    real :: TMP1         !temporary calculation
    real :: TMP2         !temporary calculation
    real :: TMPRAH2      !temporary calculation for aerodynamic resistances
    real :: TMPRB        !temporary calculation for rb
    real :: MOLG,FHGNEW,CWPC

    ! stability correction to below canopy resistance

    MOZG = 0.0
    MOLG = 0.0

    if (iter > 1) then
       TMP1 = KARMAN * (GRAV / TAH) * HG / (RHOAIR * CPAIR)
       if (abs(TMP1) <= MPE) TMP1 = MPE
       MOLG = -1. * FV ** 3 / TMP1
       MOZG = min((ZPD - Z0MG) / MOLG, 1.0)
    end if

    if (MOZG < 0.0) then
       FHGNEW  = (1.0 - 15.0 * MOZG) ** (-0.25)
    else
       FHGNEW  = 1.0 + 4.7 * MOZG
    end if

    if (iter == 1) then
       FHG = FHGNEW
    else
       FHG = 0.5 * (FHG + FHGNEW)
    end if

    CWPC = (CWP * VAI * HCAN * FHG) ** 0.5
    !       CWPC = (CWP*FHG)**0.5

    TMP1 = exp(-CWPC * Z0HG / HCAN)
    TMP2 = exp(-CWPC * (Z0H + ZPD) / HCAN )
    TMPRAH2 = HCAN * exp(CWPC) / CWPC * (TMP1 - TMP2)

    ! aerodynamic resistances raw and rah between heights zpd+z0h and z0hg.
    KH  = max(KARMAN * FV * (HCAN - ZPD), MPE)
    RAMG = 0.0
    RAHG = TMPRAH2 / KH
    RAWG = RAHG

    ! leaf boundary layer resistance

    TMPRB  = CWPC * 50.0 / (1.0 - exp(-CWPC / 2.0))
    RB     = TMPRB * sqrt(LK_DLEAF(lutyp) / UC)
    !       RB = 200

  end subroutine ragrb


  subroutine sfcdif1(iter   ,SFCTMP ,RHOAIR ,H      ,QAIR   , & !in
       &             ZLVL   ,ZPD    ,Z0M    ,Z0H    ,UR     , & !in
       &             MPE    ,ILOC   ,JLOC   ,                 & !in
       &             MOZ    ,MOZSGN ,FM     ,FH     ,FM2,FH2, & !inout
       &             CM     ,CH     ,FV     ,CH2     )          !out
    ! computing surface drag coefficient CM for momentum and CH for heat
    implicit none

    ! inputs
    integer,              intent(in) :: ILOC   !grid index
    integer,              intent(in) :: JLOC   !grid index
    integer,              intent(in) :: iter   !iteration index
    real,                 intent(in) :: SFCTMP !temperature at reference height (k)
    real,                 intent(in) :: RHOAIR !density air (kg/m**3)
    real,                 intent(in) :: H      !sensible heat flux (w/m2) [+ to atm]
    real,                 intent(in) :: QAIR   !specific humidity at reference height (kg/kg)
    real,                 intent(in) :: ZLVL   !reference height  (m)
    real,                 intent(in) :: ZPD    !zero plane displacement (m)
    real,                 intent(in) :: Z0H    !roughness length, sensible heat, ground (m)
    real,                 intent(in) :: Z0M    !roughness length, momentum, ground (m)
    real,                 intent(in) :: UR     !wind speed (m/s)
    real,                 intent(in) :: MPE    !prevents overflow error if division by zero

    ! in & out
    integer,           intent(inout) :: MOZSGN !number of times moz changes sign
    real,              intent(inout) :: MOZ    !Monin-Obukhov stability (z/L)
    real,              intent(inout) :: FM     !momentum stability correction, weighted by prior iters
    real,              intent(inout) :: FH     !sen heat stability correction, weighted by prior iters
    real,              intent(inout) :: FM2    !sen heat stability correction, weighted by prior iters
    real,              intent(inout) :: FH2    !sen heat stability correction, weighted by prior iters

    ! outputs
    real,                intent(out) :: CM     !drag coefficient for momentum
    real,                intent(out) :: CH     !drag coefficient for heat
    real,                intent(out) :: FV     !friction velocity (m/s)
    real,                intent(out) :: CH2    !drag coefficient for heat

    ! locals
    real    :: MOL                      !Monin-Obukhov length (m)
    real    :: TMPCM                    !temporary calculation for CM
    real    :: TMPCH                    !temporary calculation for CH
    real    :: FMNEW                    !stability correction factor, momentum, for current moz
    real    :: FHNEW                    !stability correction factor, sen heat, for current moz
    real    :: MOZOLD                   !Monin-Obukhov stability parameter from prior iteration
    real    :: TMP1,TMP2,TMP3,TMP4,TMP5 !temporary calculation
    real    :: TVIR                     !temporary virtual temperature (k)
    real    :: MOZ2                     !2/L
    real    :: TMPCM2                   !temporary calculation for CM2
    real    :: TMPCH2                   !temporary calculation for CH2
    real    :: FM2NEW                   !stability correction factor, momentum, for current moz
    real    :: FH2NEW                   !stability correction factor, sen heat, for current moz
    real    :: TMP12,TMP22,TMP32        !temporary calculation

    real    :: CMFM, CHFH, CM2FM2, CH2FH2

    ! Monin-Obukhov stability parameter moz for next iteration

    MOZOLD = MOZ

    if (ZLVL <= ZPD) then
       write(*,*) 'critical problem: ZLVL <= ZPD; model stops'
       call wrf_error_fatal("STOP in Noah-MP")
    end if

    TMPCM = log((ZLVL-ZPD) / Z0M)
    TMPCH = log((ZLVL-ZPD) / Z0H)
    TMPCM2 = log((2.0 + Z0M) / Z0M)
    TMPCH2 = log((2.0 + Z0H) / Z0H)

    if (iter == 1) then
       FV   = 0.0
       MOZ  = 0.0
       MOL  = 0.0
       MOZ2 = 0.0
    else
       TVIR = (1.0 + 0.61 * QAIR) * SFCTMP
       TMP1 = KARMAN * (GRAV / TVIR) * H / (RHOAIR * CPAIR)
       if (abs(TMP1) <= MPE) TMP1 = MPE
       MOL  = -1.0 * FV ** 3 / TMP1
       MOZ  = min((ZLVL-ZPD) / MOL, 1.0)
       MOZ2  = min((2.0 + Z0H) / MOL, 1.0)
    end if

    ! accumulate number of times moz changes sign.

    if (MOZOLD*MOZ < 0.0) MOZSGN = MOZSGN + 1
    if (MOZSGN >= 2) then
       MOZ = 0.0
       FM = 0.0
       FH = 0.0
       MOZ2 = 0.0
       FM2 = 0.0
       FH2 = 0.0
    end if

    ! evaluate stability-dependent variables using moz from prior iteration
    if (MOZ < 0.0) then
       TMP1 = (1.0 - 16.0 * MOZ) ** 0.25
       TMP2 = log((1.0 + TMP1 * TMP1) / 2.0)
       TMP3 = log((1.0 + TMP1) / 2.0)
       FMNEW = 2.0 * TMP3 + TMP2 - 2.0 * atan(TMP1) + 1.5707963
       FHNEW = 2.0 * TMP2

       ! 2-meter
       TMP12 = (1.0 - 16.0 * MOZ2) ** 0.25
       TMP22 = log((1.0 + TMP12 * TMP12) / 2.0)
       TMP32 = log((1.0 + TMP12) / 2.0)
       FM2NEW = 2.0 * TMP32 + TMP22 - 2.0 * atan(TMP12) + 1.5707963
       FH2NEW = 2.0 * TMP22
    else
       FMNEW = -5.0 * MOZ
       FHNEW = FMNEW
       FM2NEW = -5.0 * MOZ2
       FH2NEW = FM2NEW
    end if

    ! except for first iteration, weight stability factors for previous
    ! iteration to help avoid flip-flops from one iteration to the next

    if (iter == 1) then
       FM = FMNEW
       FH = FHNEW
       FM2 = FM2NEW
       FH2 = FH2NEW
    else
       FM = 0.5 * (FM + FMNEW)
       FH = 0.5 * (FH + FHNEW)
       FM2 = 0.5 * (FM2 + FM2NEW)
       FH2 = 0.5 * (FH2 + FH2NEW)
    end if

    ! exchange coefficients

    FH = min(FH, 0.9 * TMPCH)
    FM = min(FM, 0.9 * TMPCM)
    FH2 = min(FH2, 0.9 * TMPCH2)
    FM2 = min(FM2, 0.9 * TMPCM2)

    CMFM = TMPCM - FM
    CHFH = TMPCH - FH
    CM2FM2 = TMPCM2 - FM2
    CH2FH2 = TMPCH2 - FH2
    if (abs(CMFM) <= MPE) CMFM = MPE
    if (abs(CHFH) <= MPE) CHFH = MPE
    if (abs(CM2FM2) <= MPE) CM2FM2 = MPE
    if (abs(CH2FH2) <= MPE) CH2FH2 = MPE
    CM  = KARMAN * KARMAN / (CMFM * CMFM)
    CH  = KARMAN * KARMAN / (CMFM * CHFH)
    CH2  = KARMAN * KARMAN / (CM2FM2 * CH2FH2)

    ! friction velocity

    FV = UR * sqrt(CM)
    CH2  = KARMAN * FV / CH2FH2

  end subroutine sfcdif1


  subroutine sfcdif2(iter   ,Z0     ,THZ0   ,THLM   ,SFCSPD , & !in
       CZIL   ,ZLM    ,ILOC   ,JLOC   ,         & !in
       AKMS   ,AKHS   ,RLMO   ,WSTAR2 ,         & !in
       USTAR  )                                   !out
    ! SUBROUTINE SFCDIF (renamed SFCDIF_off to avoid clash with Eta PBL)
    ! CALCULATE SURFACE LAYER EXCHANGE COEFFICIENTS VIA ITERATIVE PROCESS.
    ! SEE CHEN ET AL (1997, BLM)
    implicit none
    integer, intent(in) :: ILOC
    integer, intent(in) :: JLOC
    integer, intent(in) :: iter
    real,    intent(in) :: ZLM, Z0, THZ0, THLM, SFCSPD, CZIL
    real, intent(inout) :: AKMS
    real, intent(inout) :: AKHS
    real, intent(inout) :: RLMO
    real, intent(inout) :: WSTAR2
    real,   intent(out) :: USTAR

    real     ZZ, PSLMU, PSLMS, PSLHU, PSLHS
    real     XX, PSPMU, YY, PSPMS, PSPHU, PSPHS
    real     ZILFC, ZU, ZT, RDZ, CXCH
    real     DTHV, DU2, BTGH, ZSLU, ZSLT, RLOGU, RLOGT
    real     ZETALT, ZETALU, ZETAU, ZETAT, XLU4, XLT4, XU4, XT4

    real     XLU, XLT, XU, XT, PSMZ, SIMM, PSHZ, SIMH, USTARK, RLMN,  &
         &         RLMA

    integer  ILECH, ITR

    integer, parameter :: ITRMX  = 5
    real,    parameter :: WWST   = 1.2
    real,    parameter :: WWST2  = WWST * WWST
    real,    parameter :: VKRM   = 0.40
    real,    parameter :: EXCM   = 0.001
    real,    parameter :: BETA   = 1.0 / 270.0
    real,    parameter :: BTG    = BETA * GRAV
    real,    parameter :: ELFC   = VKRM * BTG
    real,    parameter :: WOLD   = 0.15
    real,    parameter :: WNEW   = 1.0 - WOLD
    real,    parameter :: PIHF   = 3.14159265 / 2.
    real,    parameter :: EPSU2  = 1.E-4
    real,    parameter :: EPSUST = 0.07
    real,    parameter :: EPSIT  = 1.E-4
    real,    parameter :: EPSA   = 1.E-8
    real,    parameter :: ZTMIN  = -5.0
    real,    parameter :: ZTMAX  = 1.0
    real,    parameter :: HPBL   = 1000.0
    real,    parameter :: SQVISC = 258.2
    real,    parameter :: RIC    = 0.183
    real,    parameter :: RRIC   = 1.0 / RIC
    real,    parameter :: FHNEU  = 0.8
    real,    parameter :: RFC    = 0.191
    real,    parameter :: RFAC   = RIC / ( FHNEU * RFC * RFC )

    ! NOTE: THE TWO CODE BLOCKS BELOW DEFINE FUNCTIONS
    ! LECH'S SURFACE FUNCTIONS
    PSLMU (ZZ)= -0.96 * log(1.0 - 4.5 * ZZ)
    PSLMS (ZZ)= ZZ * RRIC -2.076 * (1.0 - 1.0 / (ZZ + 1.0))
    PSLHU (ZZ)= -0.96 * log(1.0 - 4.5 * ZZ)
    PSLHS (ZZ)= ZZ * RFAC -2.076 * (1.0 - 1.0 / (ZZ + 1.0))
    ! PAULSON'S SURFACE FUNCTIONS
    PSPMU (XX)= -2.0 * log((XX +1.0) * 0.5) - log((XX * XX + 1.0) * 0.5)   &
         &        + 2.0 * atan(XX)                                            &
         &- PIHF
    PSPMS (YY)= 5.0 * YY
    PSPHU (XX)= -2.0 * log((XX * XX + 1.0) * 0.5)
    PSPHS (YY)= 5.0 * YY

    ! THIS ROUTINE SFCDIF CAN HANDLE BOTH OVER OPEN water (SEA, OCEAN) AND
    ! OVER SOLID SURFACE (LAND, SEA-ICE).
    !
    !     ZTFC: RATIO OF ZOH/ZOM  LESS OR EQUAL THAN 1
    !     C......ZTFC=0.1
    !     CZIL: CONSTANT C IN Zilitinkevich, S. S.1995,:NOTE ABOUT ZT
    !
    ILECH = 0

    !
    ZILFC = -CZIL * VKRM * SQVISC
    ZU = Z0
    RDZ = 1.0 / ZLM
    CXCH = EXCM * RDZ
    DTHV = THLM - THZ0

    ! BELJARS CORRECTION OF USTAR
    DU2 = max(SFCSPD * SFCSPD, EPSU2)
    BTGH = BTG * HPBL

    if (iter == 1) then
       if (BTGH * AKHS * DTHV /= 0.0) then
          WSTAR2 = WWST2  * abs(BTGH * AKHS * DTHV) ** (2.0 / 3.0)
       else
          WSTAR2 = 0.0
       end if
       USTAR = max(sqrt(AKMS * sqrt(DU2 + WSTAR2)), EPSUST)
       RLMO = ELFC * AKHS * DTHV / USTAR ** 3
    end if

    ! ZILITINKEVITCH APPROACH FOR ZT
    ZT = max(1.0E-6, exp(ZILFC * sqrt(USTAR * Z0)) * Z0)
    ZSLU = ZLM + ZU
    ZSLT = ZLM + ZT
    RLOGU = log(ZSLU / ZU)
    RLOGT = log(ZSLT / ZT)

    ! 1./MONIN-OBUKKHOV LENGTH-SCALE
    ZETALT = max(ZSLT * RLMO, ZTMIN)
    RLMO = ZETALT / ZSLT
    ZETALU = ZSLU * RLMO
    ZETAU = ZU * RLMO
    ZETAT = ZT * RLMO

    if (ILECH == 0) then
       if (RLMO < 0.0) then
          XLU4 = 1.0 - 16.0 * ZETALU
          XLT4 = 1.0 - 16.0 * ZETALT
          XU4  = 1.0 - 16.0 * ZETAU
          XT4  = 1.0 - 16.0 * ZETAT
          XLU  = sqrt(sqrt(XLU4))
          XLT  = sqrt(sqrt(XLT4))
          XU   = sqrt(sqrt(XU4))

          XT = sqrt (sqrt (XT4))
          PSMZ = PSPMU (XU)
          SIMM = PSPMU (XLU) - PSMZ + RLOGU
          PSHZ = PSPHU (XT)
          SIMH = PSPHU (XLT) - PSHZ + RLOGT
       else
          ZETALU = min (ZETALU,ZTMAX)
          ZETALT = min (ZETALT,ZTMAX)
          PSMZ = PSPMS (ZETAU)
          SIMM = PSPMS (ZETALU) - PSMZ + RLOGU
          PSHZ = PSPHS (ZETAT)
          SIMH = PSPHS (ZETALT) - PSHZ + RLOGT
       end if
       ! LECH'S FUNCTIONS
    else
       if (RLMO < 0.0) then
          PSMZ = PSLMU (ZETAU)
          SIMM = PSLMU (ZETALU) - PSMZ + RLOGU
          PSHZ = PSLHU (ZETAT)
          SIMH = PSLHU (ZETALT) - PSHZ + RLOGT
       else
          ZETALU = min (ZETALU,ZTMAX)
          ZETALT = min (ZETALT,ZTMAX)
          PSMZ = PSLMS (ZETAU)
          SIMM = PSLMS (ZETALU) - PSMZ + RLOGU
          PSHZ = PSLHS (ZETAT)
          SIMH = PSLHS (ZETALT) - PSHZ + RLOGT
       end if
    end if

    ! BELJAARS CORRECTION FOR USTAR
    USTAR = max (sqrt (AKMS * sqrt (DU2+ WSTAR2)),EPSUST)

    ! ZILITINKEVITCH FIX FOR ZT
    ZT = max(1.E-6,exp (ZILFC * sqrt (USTAR * Z0))* Z0)
    ZSLT = ZLM + ZT
    !
    RLOGT = log (ZSLT / ZT)
    USTARK = USTAR * VKRM
    AKMS = max (USTARK / SIMM,CXCH)
    !
    ! IF STATEMENTS TO AVOID TANGENT LINEAR PROBLEMS NEAR ZERO
    !
    AKHS = max (USTARK / SIMH,CXCH)

    if (BTGH * AKHS * DTHV /= 0.0) then
       WSTAR2 = WWST2* abs (BTGH * AKHS * DTHV)** (2.0 / 3.0)
    else
       WSTAR2 = 0.0
    end if
    !
    RLMN = ELFC * AKHS * DTHV / USTAR **3

    RLMA = RLMO * WOLD+ RLMN * WNEW
    RLMO = RLMA

  end subroutine sfcdif2


  subroutine esat(T, ESW, ESI, DESW, DESI)
    ! use polynomials to calculate saturation vapor pressure and derivative with
    ! respect to temperature: over water when t > 0 c and over ice when t <= 0 c
    implicit none
    ! in
    real, intent(in)  :: T              !temperature

    !out
    real, intent(out) :: ESW            !saturation vapor pressure over water (pa)
    real, intent(out) :: ESI            !saturation vapor pressure over ice (pa)
    real, intent(out) :: DESW           !d(esat)/dt over water (pa/K)
    real, intent(out) :: DESI           !d(esat)/dt over ice (pa/K)

    ! local
    real :: A0,A1,A2,A3,A4,A5,A6  !coefficients for esat over water
    real :: B0,B1,B2,B3,B4,B5,B6  !coefficients for esat over ice
    real :: C0,C1,C2,C3,C4,C5,C6  !coefficients for dsat over water
    real :: D0,D1,D2,D3,D4,D5,D6  !coefficients for dsat over ice

    parameter (A0=6.107799961    , A1=4.436518521E-01,  &
         A2=1.428945805E-02, A3=2.650648471E-04,  &
         A4=3.031240396E-06, A5=2.034080948E-08,  &
         A6=6.136820929E-11)

    parameter (B0=6.109177956    , B1=5.034698970E-01,  &
         B2=1.886013408E-02, B3=4.176223716E-04,  &
         B4=5.824720280E-06, B5=4.838803174E-08,  &
         B6=1.838826904E-10)

    parameter (C0= 4.438099984E-01, C1=2.857002636E-02,  &
         C2= 7.938054040E-04, C3=1.215215065E-05,  &
         C4= 1.036561403E-07, C5=3.532421810e-10,  &
         C6=-7.090244804E-13)

    parameter (D0=5.030305237E-01, D1=3.773255020E-02,  &
         D2=1.267995369E-03, D3=2.477563108E-05,  &
         D4=3.005693132E-07, D5=2.158542548E-09,  &
         D6=7.131097725E-12)

    ESW  = 100.*(A0+T*(A1+T*(A2+T*(A3+T*(A4+T*(A5+T*A6))))))
    ESI  = 100.*(B0+T*(B1+T*(B2+T*(B3+T*(B4+T*(B5+T*B6))))))
    DESW = 100.*(C0+T*(C1+T*(C2+T*(C3+T*(C4+T*(C5+T*C6))))))
    DESI = 100.*(D0+T*(D1+T*(D2+T*(D3+T*(D4+T*(D5+T*D6))))))

  end subroutine esat

  
  subroutine stomata(lutyp, igs, sfcprs, sfctmp, apar, tv, &
       &             ea, ei, o2, co2, foln, btran, &
       &             rb, rs, psn)
    use noahmp_const, only: MPE
    use noahmp_const, only: RGAS
    use noahmp_const, only: TFRZ
    use noahmp_veg_param, only: LK_C3C4
    use noahmp_veg_param, only: LK_KC25
    use noahmp_veg_param, only: LK_AKC
    use noahmp_veg_param, only: LK_KO25
    use noahmp_veg_param, only: LK_AKO
    use noahmp_veg_param, only: LK_VCMX25
    use noahmp_veg_param, only: LK_AVCMX
    use noahmp_veg_param, only: LK_BP
    use noahmp_veg_param, only: LK_MP
    use noahmp_veg_param, only: LK_QE25
    use noahmp_veg_param, only: LK_FOLNMX
    implicit none
    ! input
    integer,intent(in) :: lutyp   !vegetation physiology type

    real(r4), intent(in) :: igs   !growing season index (0=off, 1=on)

    real(r4), intent(in) :: sfcprs  !air pressure at reference height (pa)
    real(r4), intent(in) :: sfctmp  !air temperature at reference height (k)
    real(r4), intent(in) :: apar    !par absorbed per unit lai (w/m2)
    real(r4), intent(in) :: tv      !foliage temperature (k)
    real(r4), intent(in) :: ea      !vapor pressure of canopy air (pa)
    real(r4), intent(in) :: ei      !vapor pressure inside leaf (sat vapor press at tv) (pa)
    real(r4), intent(in) :: o2      !atmospheric o2 concentration (pa)
    real(r4), intent(in) :: co2     !atmospheric co2 concentration (pa)
    real(r4), intent(in) :: foln    !foliage nitrogen concentration (%)
    real(r4), intent(in) :: btran   !soil water transpiration factor (0 to 1)
    real(r4), intent(in) :: rb      !boundary layer resistance (s/m)

    ! output
    real(r4), intent(out) :: rs     !leaf stomatal resistance (s/m)
    real(r4), intent(out) :: psn    !foliage photosynthesis (umol co2 /m2/ s) [always +]

    ! locals
    real(r4) :: rlb    !boundary layer resistance (s m2 / umol)

    real(r4) :: ci                        !internal co2 (pa)
    real(r4), parameter :: CIERR = 5.0E-2 !threshold of terminating the bisection (Pa)
    real(r4) :: cihigh, cilow             !intermediate inner leaf CO2 pressure (Pa)
    real(r4) :: fcihigh, fcilow, fci      !intermediate inner leaf CO2 pressure (Pa)
    integer :: iter                       !iteration index
    integer, parameter :: NITER = 20      !number of iterations

    real(r4) :: tc          !foliage temperature (degree Celsius)
    real(r4) :: kc          !co2 Michaelis-Menten constant (pa)
    real(r4) :: ko          !o2 Michaelis-Menten constant (pa)
    real(r4) :: fnf         !foliage nitrogen adjustment factor (0 to 1)
    real(r4) :: ppf         !absorb photosynthetic photon flux (umol photons/m2/s)
    real(r4) :: cp          !co2 compensation point (pa)
    real(r4) :: awc         !intermediate calculation for wc
    real(r4) :: vcmx        !maximum rate of carbonylation (umol co2/m2/s)
    real(r4) :: j           !electron transport (umol co2/m2/s)
    real(r4) :: cf          !s m2/umol -> s/m

    ! initialize RS=RSMAX and PSN=0 because will only do calculations
    ! for APAR > 0, in which case RS <= RSMAX and PSN >= 0
    cf = sfcprs / (RGAS * sfctmp) * 1.0e06
    rs = 1.0 / LK_BP(lutyp) * cf
    psn = 0.0
    ci = co2

    if (apar <= 0.0) return

    fnf = min(foln / max(MPE, LK_FOLNMX(lutyp)), 1.0 )
    tc  = tv - TFRZ
    ppf = 4.6 * apar
    j   = ppf * LK_QE25(lutyp)
    ! F1 = AB ** ((BC - 25.0) / 10.0)
    ! KC  = KC25(lutyp) * F1(AKC(lutyp),TC)
    ! KO  = KO25(lutyp) * F1(AKO(lutyp),TC)
    kc  = LK_KC25(lutyp) * LK_AKC(lutyp) ** ((tc - 25.0) / 10.0)
    ko  = LK_KO25(lutyp) * LK_AKO(lutyp) ** ((tc - 25.0) / 10.0)
    awc = kc * (1.0 + o2 / ko)
    cp  = 0.5 * kc / ko * o2 * 0.21
    ! VCMX = VCMX25(lutyp) / F2(TC) * FNF * BTRAN * F1(AVCMX(lutyp),TC)
    ! F2 = 1.0 + EXP((-2.2E05 + 710.0 * (AB + 273.16)) / (8.314 * (AB + 273.16)))
    vcmx = LK_VCMX25(lutyp) &
         & / (1.0 + exp((-2.2E05 + 710.0 * (tc + TFRZ)) / (8.314 * (tc + TFRZ)))) &
         & * fnf * btran &
         & * (LK_AVCMX(lutyp) ** ((tc - 25.0) / 10.0))

    ! rb: s/m -> s m**2 / umol
    rlb = rb / cf

    ! endpoints of the search intervals
    cihigh = 1.5 * co2
    cilow = 0.0
    do iter = 1, niter
       ci = 0.5 * (cihigh + cilow)
       call ci2ci(ci, fci, RS, PSN)
       if (((cihigh - cilow) <= CIERR) .or. abs(fci - ci) <= MPE) then
          exit
       elseif (fci > ci) then
          cilow = ci
       else
          cihigh = ci
       end if
    end do

    rs = rs * cf

  contains
    subroutine ci2ci(ci, fci, rs, psn)
      !function for serching the fixed point of CI, that is CI = fci(CI)
      implicit none
      real(r4), intent(in) :: ci
      real(r4), intent(out) :: fci
      real(r4), intent(out) :: rs
      real(r4), intent(out) :: psn
      real(r4) :: wc = nan4   !Rubisco limited photosynthesis (umol co2/m2/s)
      real(r4) :: wj = nan4   !light limited photosynthesis (umol co2/m2/s)
      real(r4) :: we = nan4   !export limited photosynthesis (umol co2/m2/s)
      real(r4) :: cs          !co2 concentration at leaf surface (pa)
      real(r4) :: a, b, c, q  !intermediate calculations for RS
      real(r4) :: r1, r2      !roots for RS

      if (LK_C3C4(lutyp) == 1) then ! C3
        wj = max(ci - cp, 0.0) * j / (ci + 2.0 * cp)
        wc = max(ci - cp, 0.0) * vcmx / (ci + awc)
        we = 0.5 * vcmx
      elseif (LK_C3C4(lutyp) == 2) then  ! C4
        wj = j
        wc = vcmx
        we = 4000.0 * vcmx * ci / sfcprs
      end if
      psn = min(wj, wc, we) * igs

      cs = max(co2 - 1.37 * rlb * sfcprs * psn, MPE )
      a = LK_MP(lutyp) * psn * sfcprs * ea / (cs * ei) + LK_BP(lutyp)
      b = (LK_MP(lutyp) * psn * sfcprs / cs + LK_BP(lutyp)) * rlb - 1.0
      c = -rlb
      if (b >= 0.0) then
         q = -0.5 * (b + sqrt(b * b - 4.0 * a * c))
      else
         q = -0.5 * (b - sqrt(b * b - 4.0 * a * c))
      end if
      r1 = q / a
      r2 = c / q
      rs = max(r1, r2)

      fci = max(cs - psn * sfcprs * 1.65 * rs, 0.0)
    end subroutine ci2ci
  end subroutine stomata


  subroutine canres(lutyp, sfcprs, tv, par, eah, btran, rs, psn)
    ! calculate canopy resistance which depends on incoming solar radiation,
    ! air temperature, atmospheric water vapor pressure deficit at the
    ! lowest model level, and soil moisture (preferably unfrozen soil
    ! moisture rather than total)
    !
    ! References:
    !     Jarvis, 1976, doi:10.1098/rstb.1976.0035
    !     Noilhan and Planton, 1989, MWR, doi:10.1175/1520-0493(1989)117<0536:ASPOLS>2.0.CO;2
    !     Jacquemin and Noilhan, 1990, BLM, doi:10.1007/BF00123180
    !     Chen et al., 1996, JGR, doi:10.1029/95JD02165
    !
    use noahmp_veg_param, only: LK_RGL
    use noahmp_veg_param, only: LK_RSMAX
    use noahmp_veg_param, only: LK_RSMIN
    use noahmp_veg_param, only: LK_HS
    use noahmp_veg_param, only: LK_TOPT
    implicit none
    ! inputs
    integer, intent(in) :: lutyp !vegetation type
    real(r4), intent(in) :: sfcprs  !surface pressure (pa)
    real(r4), intent(in) :: par     !par absorbed per unit sunlit lai (w/m2)
    real(r4), intent(in) :: tv      !foliage temperature (K) (FIXME: air or foliage?)
    real(r4), intent(in) :: eah     !water vapor pressure (pa)
    real(r4), intent(in) :: btran   !soil moisture stress factor

    !outputs
    real(r4), intent(out) :: rs     !canopy resistance per unit LAI
    real(r4), intent(out) :: psn    !foliage photosynthesis (umolco2/m2/s)

    !local
    real(r4) :: rcs
    real(r4) :: rcq
    real(r4) :: rct
    real(r4) :: ff
    real(r4) :: q2     !water vapor mixing ratio (kg/kg)
    real(r4) :: q2sat  !saturation Q2
    real(r4) :: dqsdt2 !d(Q2SAT)/d(T)

    ! initialize canopy conductance multiplier terms.
    rcs = 1.0
    rct = 1.0
    rcq = 1.0

    !  compute Q2 and Q2SAT
    q2 = 0.622 * eah / (sfcprs - 0.378 * eah)   !specific humidity [kg/kg]
    q2 = q2 / (1.0 + q2)                        !mixing ratio [kg/kg]
    call calhum(tv, sfcprs, q2sat, dqsdt2)

    ! contribution due to incoming solar radiation
    ff  = 2.0 * par / LK_RGL(lutyp)  !FIXME: LAI? Chen 1996
    rcs = (ff + LK_RSMIN(lutyp) / LK_RSMAX(lutyp)) / (1.0 + ff)
    rcs = min(max(rcs, 0.0001), 1.0)

    ! contribution due to air temperature
    rct = 1.0 - 0.0016 * ((LK_TOPT(lutyp) - tv) ** 2)
    rct = min(max(rct, 0.0001), 1.0)

    ! contribution due to vapor pressure deficit
    rcq = 1.0 / (1.0 + LK_HS(lutyp) * max(0.0, q2sat - q2))
    rcq = min(max(rcq, 0.01), 1.0)

    ! determine canopy resistance due to all factors
    rs = LK_RSMIN(lutyp) / (rcs * rct * rcq * btran)
    psn = nan4       ! PSN not applied for dynamic carbon
  end subroutine canres


  subroutine calhum(SFCTMP, SFCPRS, Q2SAT, DQSDT2)

    implicit none

    real, intent(in)       :: SFCTMP, SFCPRS
    real, intent(out)      :: Q2SAT, DQSDT2
    real, parameter        :: A2=17.67,A3=273.15,A4=29.65, ELWV=2.501E6,         &
         A23M4=A2*(A3-A4), E0=0.611, RV=461.0,             &
         EPSILON=0.622
    real                   :: ES, SFCPRSX

    ! Q2SAT: saturated mixing ratio
    ES = E0 * exp ( ELWV/RV*(1./A3 - 1./SFCTMP) )
    ! convert SFCPRS from Pa to KPa
    SFCPRSX = SFCPRS*1.E-3
    Q2SAT = EPSILON * ES / (SFCPRSX-ES)
    ! convert from  g/g to g/kg
    Q2SAT = Q2SAT * 1.E3
    ! Q2SAT is currently a 'mixing ratio'

    ! DQSDT2 is calculated assuming Q2SAT is a specific humidity
    DQSDT2=(Q2SAT/(1+Q2SAT))*A23M4/(SFCTMP-A4)**2

    ! DG Q2SAT needs to be in g/g when returned for SFLX
    Q2SAT = Q2SAT / 1.E3

  end subroutine calhum


  subroutine tsnosoi(ICE     ,NSOIL   ,NSNOW   ,ISNOW   ,IST     , & !in
       TBOT    ,ZSNSO   ,SSOIL   ,DF      ,HCPCT   , & !in
       ZBOT    ,SAG     ,DT      ,snowh   ,DZSNSO  , & !in
       TG      ,ILOC    ,JLOC    ,                   & !in
       STC     )                                       !inout
    ! Compute snow (up to 3L) and soil (4L) temperature. Note that snow temperatures
    ! during melting season may exceed melting point (TFRZ) but later in PHASECHANGE
    ! subroutine the snow temperatures are reset to TFRZ for melting snow.
    implicit none
    !input
    integer,                         intent(in)  :: ILOC
    integer,                         intent(in)  :: JLOC
    integer,                         intent(in)  :: ICE    !
    integer,                         intent(in)  :: NSOIL  !no of soil layers (4)
    integer,                         intent(in)  :: NSNOW  !maximum no of snow layers (3)
    integer,                         intent(in)  :: ISNOW  !actual no of snow layers
    integer,                         intent(in)  :: IST    !surface type

    real,                            intent(in)  :: DT     !time step (s)
    real,                            intent(in)  :: TBOT   !
    real,                            intent(in)  :: SSOIL  !ground heat flux (w/m2)
    real,                            intent(in)  :: SAG    !solar rad. absorbed by ground (w/m2)
    real,                            intent(in)  :: snowh  !snow depth (m)
    real,                            intent(in)  :: ZBOT   !from soil surface (m)
    real,                            intent(in)  :: TG     !ground temperature (k)
    real, dimension(-NSNOW+1:NSOIL), intent(in)  :: ZSNSO  !layer-bot. depth from snow surf.(m)
    real, dimension(-NSNOW+1:NSOIL), intent(in)  :: DZSNSO !snow/soil layer thickness (m)
    real, dimension(-NSNOW+1:NSOIL), intent(in)  :: DF     !thermal conductivity
    real, dimension(-NSNOW+1:NSOIL), intent(in)  :: HCPCT  !heat capacity (J/m3/k)

    !input and output
    real, dimension(-NSNOW+1:NSOIL), intent(inout) :: STC

    !local
    integer                                      :: IZ
    real                                         :: ZBOTSNO   !ZBOT from snow surface
    real, dimension(-NSNOW+1:NSOIL)              :: AI, BI, CI, RHSTS
    real                                         :: EFLXB !energy influx from soil bottom (w/m2)
    real, dimension(-NSNOW+1:NSOIL)              :: PHI   !light through water (w/m2)

    real, dimension(-NSNOW+1:NSOIL) :: TBEG
    real                            :: ERR_EST !heat storage error  (w/m2)
    real                            :: SSOIL2  !ground heat flux (w/m2) (for energy check)
    real                            :: EFLXB2  !heat flux from the bottom (w/m2) (for energy check)
    character(len=256)              :: message

    ! compute solar penetration through water, needs more work

    PHI(ISNOW+1:NSOIL) = 0.0

    ! adjust ZBOT from soil surface to ZBOTSNO from snow surface

    ZBOTSNO = ZBOT - snowh    !from snow surface

    ! snow/soil heat storage for energy balance check

    do IZ = ISNOW+1, NSOIL
       TBEG(IZ) = STC(IZ)
    end do

    ! compute soil temperatures

    call HRT   (NSNOW     ,NSOIL     ,ISNOW     ,ZSNSO     , &
         STC       ,TBOT      ,ZBOTSNO   ,DT        , &
         DF        ,HCPCT     ,SSOIL     ,PHI       , &
         AI        ,BI        ,CI        ,RHSTS     , &
         EFLXB     )

    call HSTEP (NSNOW     ,NSOIL     ,ISNOW     ,DT        , &
         AI        ,BI        ,CI        ,RHSTS     , &
         STC       )

    ! update ground heat flux just for energy check, but not for final output
    ! otherwise, it would break the surface energy balance

    if (opt_tbot == 1) then
       EFLXB2  = 0.
    else if (opt_tbot == 2) then
       EFLXB2  = DF(NSOIL)*(TBOT-STC(NSOIL)) / &
            (0.5*(ZSNSO(NSOIL-1)+ZSNSO(NSOIL)) - ZBOTSNO)
    end if

    ! Skip the energy balance check for now, until we can make it work
    ! right for small time steps.
    return

    ! energy balance check

    ERR_EST = 0.0
    do IZ = ISNOW+1, NSOIL
       ERR_EST = ERR_EST + (STC(IZ)-TBEG(IZ)) * DZSNSO(IZ) * HCPCT(IZ) / DT
    end do

    if (opt_stc == 1) then   ! semi-implicit
       ERR_EST = ERR_EST - (SSOIL +EFLXB)
    else                     ! full-implicit
       SSOIL2 = DF(ISNOW+1)*(TG-STC(ISNOW+1))/(0.5*DZSNSO(ISNOW+1))   !M. Barlage
       ERR_EST = ERR_EST - (SSOIL2+EFLXB2)
    end if

    if (abs(ERR_EST) > 1.) then    ! W/m2
       write(message,*) 'TSNOSOI is losing(-)/gaining(+) false energy',ERR_EST,' W/m2'
       call wrf_message(trim(message))
       write(message,'(i6,1x,i6,1x,i3,F18.13,5F20.12)') &
            ILOC, JLOC, IST,ERR_EST,SSOIL,snowh,TG,STC(ISNOW+1),EFLXB
       call wrf_message(trim(message))
       !niu      STOP
    end if

  end subroutine tsnosoi


  subroutine hrt(NSNOW     ,NSOIL     ,ISNOW     ,ZSNSO     , &
       STC       ,TBOT      ,ZBOT      ,DT        , &
       DF        ,HCPCT     ,SSOIL     ,PHI       , &
       AI        ,BI        ,CI        ,RHSTS     , &
       BOTFLX    )
    ! calculate the right hand side of the time tendency term of the soil
    ! thermal diffusion equation.  also to compute ( prepare ) the matrix
    ! coefficients for the tri-diagonal matrix of the implicit time scheme.
    implicit none
    ! input
    integer,                         intent(in)  :: NSOIL  !no of soil layers (4)
    integer,                         intent(in)  :: NSNOW  !maximum no of snow layers (3)
    integer,                         intent(in)  :: ISNOW  !actual no of snow layers
    real,                            intent(in)  :: TBOT   !bottom soil temp. at ZBOT (k)
    real,                            intent(in)  :: ZBOT   !depth of lower boundary condition (m)
    !from soil surface not snow surface
    real,                            intent(in)  :: DT     !time step (s)
    real,                            intent(in)  :: SSOIL  !ground heat flux (w/m2)
    real, dimension(-NSNOW+1:NSOIL), intent(in)  :: ZSNSO  !depth of layer-bottom of snow/soil (m)
    real, dimension(-NSNOW+1:NSOIL), intent(in)  :: STC    !snow/soil temperature (k)
    real, dimension(-NSNOW+1:NSOIL), intent(in)  :: DF     !thermal conductivity [w/m/k]
    real, dimension(-NSNOW+1:NSOIL), intent(in)  :: HCPCT  !heat capacity [j/m3/k]
    real, dimension(-NSNOW+1:NSOIL), intent(in)  :: PHI    !light through water (w/m2)

    ! output
    real, dimension(-NSNOW+1:NSOIL), intent(out) :: RHSTS  !right-hand side of the matrix
    real, dimension(-NSNOW+1:NSOIL), intent(out) :: AI     !left-hand side coefficient
    real, dimension(-NSNOW+1:NSOIL), intent(out) :: BI     !left-hand side coefficient
    real, dimension(-NSNOW+1:NSOIL), intent(out) :: CI     !left-hand side coefficient
    real,                            intent(out) :: BOTFLX !energy influx from soil bottom (w/m2)

    ! local
    integer                                      :: K
    real, dimension(-NSNOW+1:NSOIL)              :: DDZ
    real, dimension(-NSNOW+1:NSOIL)              :: DZ
    real, dimension(-NSNOW+1:NSOIL)              :: DENOM
    real, dimension(-NSNOW+1:NSOIL)              :: DTSDZ
    real, dimension(-NSNOW+1:NSOIL)              :: EFLUX
    real                                         :: TEMP1

    do K = ISNOW+1, NSOIL
       if (K == ISNOW+1) then
          DENOM(K)  = - ZSNSO(K) * HCPCT(K)
          TEMP1     = - ZSNSO(K+1)
          DDZ(K)    = 2.0 / TEMP1
          DTSDZ(K)  = 2.0 * (STC(K) - STC(K+1)) / TEMP1
          EFLUX(K)  = DF(K) * DTSDZ(K) - SSOIL - PHI(K)
       else if (K < NSOIL) then
          DENOM(K)  = (ZSNSO(K-1) - ZSNSO(K)) * HCPCT(K)
          TEMP1     = ZSNSO(K-1) - ZSNSO(K+1)
          DDZ(K)    = 2.0 / TEMP1
          DTSDZ(K)  = 2.0 * (STC(K) - STC(K+1)) / TEMP1
          EFLUX(K)  = (DF(K)*DTSDZ(K) - DF(K-1)*DTSDZ(K-1)) - PHI(K)
       else if (K == NSOIL) then
          DENOM(K)  = (ZSNSO(K-1) - ZSNSO(K)) * HCPCT(K)
          TEMP1     =  ZSNSO(K-1) - ZSNSO(K)
          if (opt_tbot == 1) then
             BOTFLX     = 0.
          end if
          if (opt_tbot == 2) then
             DTSDZ(K)  = (STC(K) - TBOT) / ( 0.5*(ZSNSO(K-1)+ZSNSO(K)) - ZBOT)
             BOTFLX    = -DF(K) * DTSDZ(K)
          end if
          EFLUX(K)  = (-BOTFLX - DF(K-1)*DTSDZ(K-1) ) - PHI(K)
       end if
    end do

    do K = ISNOW+1, NSOIL
       if (K == ISNOW+1) then
          AI(K)    =   0.0
          CI(K)    = - DF(K)   * DDZ(K) / DENOM(K)
          if (opt_stc == 1) then
             BI(K) = - CI(K)
          end if
          if (opt_stc == 2) then
             BI(K) = - CI(K) + DF(K)/(0.5*ZSNSO(K)*ZSNSO(K)*HCPCT(K))
          end if
       else if (K < NSOIL) then
          AI(K)    = - DF(K-1) * DDZ(K-1) / DENOM(K)
          CI(K)    = - DF(K  ) * DDZ(K  ) / DENOM(K)
          BI(K)    = - (AI(K) + CI (K))
       else if (K == NSOIL) then
          AI(K)    = - DF(K-1) * DDZ(K-1) / DENOM(K)
          CI(K)    = 0.0
          BI(K)    = - (AI(K) + CI(K))
       end if
       RHSTS(K)  = EFLUX(K)/ (-DENOM(K))
    end do

  end subroutine hrt

  subroutine hstep(NSNOW     ,NSOIL     ,ISNOW     ,DT        ,  &
       AI        ,BI        ,CI        ,RHSTS     ,  &
       STC       )
    ! CALCULATE/UPDATE THE SOIL TEMPERATURE FIELD.
    implicit none
    ! input
    integer,                         intent(in)    :: NSOIL
    integer,                         intent(in)    :: NSNOW
    integer,                         intent(in)    :: ISNOW
    real,                            intent(in)    :: DT

    ! output & input
    real, dimension(-NSNOW+1:NSOIL), intent(inout) :: RHSTS
    real, dimension(-NSNOW+1:NSOIL), intent(inout) :: AI
    real, dimension(-NSNOW+1:NSOIL), intent(inout) :: BI
    real, dimension(-NSNOW+1:NSOIL), intent(inout) :: CI
    real, dimension(-NSNOW+1:NSOIL), intent(inout) :: STC

    ! local
    integer                                        :: K
    real, dimension(-NSNOW+1:NSOIL)                :: RHSTSIN
    real, dimension(-NSNOW+1:NSOIL)                :: CIIN

    do K = ISNOW+1,NSOIL
       RHSTS(K) =   RHSTS(K) * DT
       AI(K)    =      AI(K) * DT
       BI(K)    = 1. + BI(K) * DT
       CI(K)    =      CI(K) * DT
    end do

    ! copy values for input variables before call to rosr12

    do K = ISNOW+1,NSOIL
       RHSTSIN(K) = RHSTS(K)
       CIIN(K)    = CI(K)
    end do

    ! solve the tri-diagonal matrix equation

    call rosr12(CI, AI, BI, CIIN, RHSTSIN, RHSTS, ISNOW+1, NSOIL, NSNOW)

    ! update snow & soil temperature

    do K = ISNOW+1,NSOIL
       STC(K) = STC(K) + CI(K)
    end do

  end subroutine hstep


  subroutine rosr12(P, A, B, C, D, DELTA, NTOP, NSOIL, NSNOW)
    ! SUBROUTINE rosr12
    ! INVERT (SOLVE) THE TRI-DIAGONAL MATRIX PROBLEM SHOWN BELOW:
    ! ----------------------------------------------------------------------
    ! ###                                            ### ###  ###   ###  ###
    ! #B(1), C(1),  0  ,  0  ,  0  ,   . . .  ,    0   # #      #   #      #
    ! #A(2), B(2), C(2),  0  ,  0  ,   . . .  ,    0   # #      #   #      #
    ! # 0  , A(3), B(3), C(3),  0  ,   . . .  ,    0   # #      #   # D(3) #
    ! # 0  ,  0  , A(4), B(4), C(4),   . . .  ,    0   # # P(4) #   # D(4) #
    ! # 0  ,  0  ,  0  , A(5), B(5),   . . .  ,    0   # # P(5) #   # D(5) #
    ! # .                                          .   # #  .   # = #   .  #
    ! # .                                          .   # #  .   #   #   .  #
    ! # .                                          .   # #  .   #   #   .  #
    ! # 0  , . . . , 0 , A(M-2), B(M-2), C(M-2),   0   # #P(M-2)#   #D(M-2)#
    ! # 0  , . . . , 0 ,   0   , A(M-1), B(M-1), C(M-1)# #P(M-1)#   #D(M-1)#
    ! # 0  , . . . , 0 ,   0   ,   0   ,  A(M) ,  B(M) # # P(M) #   # D(M) #
    ! ###                                            ### ###  ###   ###  ###
    ! ----------------------------------------------------------------------
    implicit none

    integer, intent(in)   :: NTOP
    integer, intent(in)   :: NSOIL, NSNOW
    integer               :: K, KK

    real, dimension(-NSNOW+1:NSOIL),intent(in):: A, B, D
    real, dimension(-NSNOW+1:NSOIL),intent(inout):: C, P, DELTA

    ! INITIALIZE EQN COEF C FOR THE LOWEST SOIL LAYER
    C(NSOIL) = 0.0
    P(NTOP) = - C(NTOP) / B(NTOP)

    ! SOLVE THE COEFS FOR THE 1ST SOIL LAYER
    DELTA(NTOP) = D(NTOP) / B(NTOP)

    ! SOLVE THE COEFS FOR SOIL LAYERS 2 THRU NSOIL
    do K = NTOP+1, NSOIL
       P(K) = - C(K) * (1.0 / (B(K) + A(K) * P(K - 1)))
       DELTA(K) = (D(K) - A (K) * DELTA(K - 1))* (1.0 / (B(K) + A(K) * P(K - 1)))
    end do

    ! SET P TO DELTA FOR LOWEST SOIL LAYER
    P(NSOIL) = DELTA(NSOIL)

    ! ADJUST P FOR SOIL LAYERS 2 THRU NSOIL
    do K = NTOP+1, NSOIL
       KK = NSOIL - K + (NTOP-1) + 1
       P(KK) = P(KK) * P(KK +1) + DELTA(KK)
    end do
  end subroutine rosr12


  subroutine phasechange(sltyp, NSNOW, NSOIL, ISNOW, DT, FACT, & !in
       & DZSNSO, HCPCT, IST, ILOC, JLOC, & !in
       & STC, SNICE, SNLIQ, sneqv, snowh, & !inout
       & SMC, soilwat,                            & !inout
       & QMELT, IMELT, PONDING)                     !out
    ! melting/freezing of snow water and soil water
    use noahmp_soil_param, only: LK_BEXP
    use noahmp_soil_param, only: LK_PSISAT
    use noahmp_soil_param, only: LK_SMCMAX
    implicit none
    ! inputs
    integer, intent(in) :: sltyp
    integer, intent(in)                             :: ILOC   !grid index
    integer, intent(in)                             :: JLOC   !grid index
    integer, intent(in)                             :: NSNOW  !maximum no. of snow layers [=3]
    integer, intent(in)                             :: NSOIL  !No. of soil layers [=4]
    integer, intent(in)                             :: ISNOW  !actual no. of snow layers [<=3]
    integer, intent(in)                             :: IST    !surface type: 1->soil; 2->lake
    real, intent(in)                                :: DT     !land model time step (sec)
    real, dimension(-NSNOW+1:NSOIL), intent(in)     :: FACT   !temporary
    real, dimension(-NSNOW+1:NSOIL), intent(in)     :: DZSNSO !snow/soil layer thickness [m]
    real, dimension(-NSNOW+1:NSOIL), intent(in)     :: HCPCT  !heat capacity (J/m3/k)

    ! outputs
    integer, dimension(-NSNOW+1:NSOIL), intent(out) :: IMELT  !phase change index
    real,                               intent(out) :: QMELT  !snowmelt rate [mm/s]
    real,                               intent(out) :: PONDING!snowmelt when snow has no layer [mm]

    ! inputs and outputs
    real, intent(inout) :: sneqv
    real, intent(inout) :: snowh
    real, dimension(-NSNOW+1:NSOIL), intent(inout)  :: STC    !snow/soil layer temperature [k]
    real, dimension(       1:NSOIL), intent(inout)  :: soilwat   !soil liquid water [m3/m3]
    real, dimension(       1:NSOIL), intent(inout)  :: SMC    !total soil water [m3/m3]
    real, dimension(-NSNOW+1:0)    , intent(inout)  :: SNICE  !snow layer ice [mm]
    real, dimension(-NSNOW+1:0)    , intent(inout)  :: SNLIQ  !snow layer liquid water [mm]

    ! local
    integer                         :: J         !do loop index
    real, dimension(-NSNOW+1:NSOIL) :: HM        !energy residual [w/m2]
    real, dimension(-NSNOW+1:NSOIL) :: XM        !melting or freezing water [kg/m2]
    real, dimension(-NSNOW+1:NSOIL) :: WMASS0
    real, dimension(-NSNOW+1:NSOIL) :: WICE0
    real, dimension(-NSNOW+1:NSOIL) :: WLIQ0
    real, dimension(-NSNOW+1:NSOIL) :: MICE      !soil/snow ice mass [mm]
    real, dimension(-NSNOW+1:NSOIL) :: MLIQ      !soil/snow liquid water mass [mm]
    real, dimension(-NSNOW+1:NSOIL) :: SUPERCOOL !supercooled water in soil (kg/m2)
    real                            :: HEATR     !energy residual or loss after melting/freezing
    real                            :: TEMP1     !temporary variables [kg/m2]
    real                            :: PROPOR
    real                            :: SMP       !frozen water potential (mm)
    real                            :: XMF       !total latent heat of phase change

    ! Initialization

    QMELT   = 0.0
    PONDING = 0.0
    XMF     = 0.0

    do J = -NSNOW+1, NSOIL
       SUPERCOOL(J) = 0.0
    end do

    do J = ISNOW+1, 0       ! all layers
       MICE(J) = SNICE(J)
       MLIQ(J) = SNLIQ(J)
    end do

    do J = 1, NSOIL               ! soil
       MLIQ(J) =  soilwat(J)            * DZSNSO(J) * 1000.0
       MICE(J) = (SMC(J) - soilwat(J))  * DZSNSO(J) * 1000.0
    end do

    do J = ISNOW+1,NSOIL       ! all layers
       IMELT(J)    = 0
       HM(J)       = 0.0
       XM(J)       = 0.0
       WICE0(J)    = MICE(J)
       WLIQ0(J)    = MLIQ(J)
       WMASS0(J)   = MICE(J) + MLIQ(J)
    end do

    if (ist == 1) then
       do J = 1, NSOIL
          if (opt_frz == 1) then
             if (STC(J) < TFRZ) then
                SMP = HFUS * (TFRZ-STC(J)) / (GRAV * STC(J))             !(m)
                SUPERCOOL(J) = LK_SMCMAX(sltyp) * (SMP / LK_PSISAT(sltyp)) ** (-1.0 / LK_BEXP(sltyp))
                SUPERCOOL(J) = SUPERCOOL(J) * DZSNSO(J) * 1000.0        !(mm)
             end if
          end if
          if (opt_frz == 2) then
             call frh2o(sltyp, SUPERCOOL(J), STC(J), SMC(J), soilwat(J))
             SUPERCOOL(J) = SUPERCOOL(J) * DZSNSO(J) * 1000.0        !(mm)
          end if
       end do
    end if

    do J = ISNOW+1,NSOIL
       if (MICE(J) > 0.0 .and. STC(J) >= TFRZ) then  !melting
          IMELT(J) = 1
       end if
       if (MLIQ(J) > SUPERCOOL(J) .and. STC(J) < TFRZ) then
          IMELT(J) = 2
       end if

       ! If snow exists, but its thickness is not enough to create a layer
       if (ISNOW == 0 .and. sneqv > 0.0 .and. J == 1) then
          if (STC(J) >= TFRZ) then
             IMELT(J) = 1
          end if
       end if
    end do

    ! Calculate the energy surplus and loss for melting and freezing
    do J = ISNOW+1,NSOIL
       if (IMELT(J) > 0) then
          HM(J) = (STC(J)-TFRZ)/FACT(J)
          STC(J) = TFRZ
       end if

       if (IMELT(J) == 1 .and. HM(J) < 0.0) then
          HM(J) = 0.0
          IMELT(J) = 0
       end if
       if (IMELT(J) == 2 .and. HM(J) > 0.0) then
          HM(J) = 0.0
          IMELT(J) = 0
       end if
       XM(J) = HM(J)*DT/HFUS
    end do

    ! The rate of melting and freezing for snow without a layer, needs more work.
    if (ISNOW == 0 .and. sneqv > 0.0 .and. XM(1) > 0.0) then
       TEMP1  = sneqv
       sneqv  = max(0.0, TEMP1-XM(1))
       PROPOR = sneqv / TEMP1
       snowh  = max(0.0, PROPOR * snowh)
       HEATR  = HM(1) - HFUS * (TEMP1 - sneqv) / DT
       if (HEATR > 0.0) then
          XM(1) = HEATR * DT / HFUS
          HM(1) = HEATR
       else
          XM(1) = 0.0
          HM(1) = 0.0
       end if
       QMELT   = max(0.0, (TEMP1 - sneqv)) / DT
       XMF     = HFUS * QMELT
       PONDING = TEMP1 - sneqv
    end if

    ! The rate of melting and freezing for snow and soil
    do J = ISNOW+1,NSOIL
       if (IMELT(J) > 0 .and. abs(HM(J)) > 0.0) then

          HEATR = 0.0
          if (XM(J) > 0.0) then
             MICE(J) = max(0.0, WICE0(J) - XM(J))
             HEATR = HM(J) - HFUS * (WICE0(J) - MICE(J)) / DT
          else if (XM(J) < 0.0) then
             if (J <= 0) then                             ! snow
                MICE(J) = min(WMASS0(J), WICE0(J) - XM(J))
             else                                         ! soil
                if (WMASS0(J) < SUPERCOOL(J)) then
                   MICE(J) = 0.0
                else
                   MICE(J) = min(WMASS0(J) - SUPERCOOL(J), WICE0(J) - XM(J))
                   MICE(J) = max(MICE(J), 0.0)
                end if
             end if
             HEATR = HM(J) - HFUS * (WICE0(J) - MICE(J)) / DT
          end if

          MLIQ(J) = max(0.0, WMASS0(J) - MICE(J))

          if (abs(HEATR) > 0.0) then
             STC(J) = STC(J) + FACT(J) * HEATR
             if (J <= 0) then                             ! snow
                if (MLIQ(J) * MICE(J) > 0.0) STC(J) = TFRZ
             end if
          end if

          XMF = XMF + HFUS * (WICE0(J) - MICE(J)) / DT

          if (J < 1) then
             QMELT = QMELT + max(0.0, (WICE0(J) - MICE(J))) / DT
          end if
       end if
    end do

    do J = ISNOW+1, 0             ! snow
       SNLIQ(J) = MLIQ(J)
       SNICE(J) = MICE(J)
    end do

    do J = 1, NSOIL              ! soil
       soilwat(J) =  MLIQ(J)            / (1000.0 * DZSNSO(J))
       SMC(J)  = (MLIQ(J) + MICE(J)) / (1000.0 * DZSNSO(J))
    end do

  end subroutine phasechange


  subroutine frh2o(sltyp, FREE, TKELV, SMC, soilwat)
    ! SUBROUTINE FRH2O
    !
    ! CALCULATE AMOUNT OF SUPERCOOLED LIQUID SOIL water CONTENT IF
    ! TEMPERATURE IS BELOW 273.15K (TFRZ).  REQUIRES NEWTON-TYPE ITERATION
    ! TO SOLVE THE NONLINEAR IMPLICIT EQUATION GIVEN IN EQN 17 OF KOREN ET AL
    ! (1999, JGR, VOL 104(D16), 19569-19585).
    !
    ! NEW VERSION (JUNE 2001): MUCH FASTER AND MORE ACCURATE NEWTON
    ! ITERATION ACHIEVED BY FIRST TAKING LOG OF EQN CITED ABOVE -- LESS THAN
    ! 4 (TYPICALLY 1 OR 2) ITERATIONS ACHIEVES CONVERGENCE.  ALSO, EXPLICIT
    ! 1-STEP SOLUTION OPTION FOR SPECIAL CASE OF PARAMETER CK=0, WHICH
    ! REDUCES THE ORIGINAL IMPLICIT EQUATION TO A SIMPLER EXPLICIT FORM,
    ! KNOWN AS THE "FLERCHINGER EQN". IMPROVED HANDLING OF SOLUTION IN THE
    ! LIMIT OF FREEZING POINT TEMPERATURE TFRZ.
    !
    ! INPUT:

    !   TKELV.........TEMPERATURE (Kelvin)
    !   SMC...........TOTAL SOIL MOISTURE CONTENT (VOLUMETRIC)
    !   soilwat..........LIQUID SOIL MOISTURE CONTENT (VOLUMETRIC)
    !   B.............SOIL TYPE "B" PARAMETER
    !   PSISAT........SATURATED SOIL MATRIC POTENTIAL

    ! OUTPUT:
    !   FREE..........SUPERCOOLED LIQUID water CONTENT [m3/m3]
    use noahmp_soil_param, only: LK_BEXP
    use noahmp_soil_param, only: LK_PSISAT
    use noahmp_soil_param, only: LK_SMCMAX
    implicit none
    integer, intent(in) :: sltyp
    real, intent(in)    :: soilwat
    real, intent(in)    :: SMC
    real, intent(in)    :: TKELV
    real, intent(out)   :: FREE
    real                :: BX,DENOM,DF,DSWL,FK,SWL,SWLK
    integer             :: NLOG,KCOUNT
    !      PARAMETER(CK = 0.0)
    real, parameter     :: CK = 8.0, BLIM = 5.5, ERROR = 0.005,       &
         DICE = 920.0
    character(80)       :: message

    ! LIMITS ON PARAMETER B: B < 5.5  (use parameter BLIM)
    ! SIMULATIONS SHOWED IF B > 5.5 UNFROZEN water CONTENT IS
    ! NON-REALISTICALLY HIGH AT VERY LOW TEMPERATURES.
    BX = LK_BEXP(sltyp)

    ! INITIALIZING ITERATIONS COUNTER AND ITERATIVE SOLUTION FLAG.
    if (LK_BEXP(sltyp) > BLIM) BX = BLIM
    NLOG = 0

    ! IF TEMPERATURE NOT SIGNIFICANTLY BELOW FREEZING (TFRZ), soilwat = SMC
    KCOUNT = 0
    if (TKELV > (TFRZ - 1.0E-3)) then
       FREE = SMC
    else
       ! OPTION 1: ITERATED SOLUTION IN KOREN ET AL, JGR, 1999, EQN 17
       ! INITIAL GUESS FOR SWL (frozen content)
       if (CK /= 0.0) then
          SWL = SMC - soilwat

          ! KEEP WITHIN BOUNDS.
          if (SWL > (SMC - 0.02)) SWL = SMC -0.02

          ! START OF ITERATIONS
          if (SWL < 0.0) SWL = 0.0
          do while ((NLOG < 10) .and. (KCOUNT == 0))
             NLOG = NLOG + 1
             DF = ALOG((LK_PSISAT(sltyp) * GRAV / hfus ) * ((1.0 + CK * SWL ) ** 2.0) * &
                  & (LK_SMCMAX(sltyp) / (SMC - SWL)) ** BX) - ALOG(-(TKELV - TFRZ) / TKELV)
             DENOM = 2.0 * CK / (1.0 + CK * SWL) + BX / (SMC - SWL)
             SWLK = SWL - DF / DENOM

             ! BOUNDS USEFUL FOR MATHEMATICAL SOLUTION.
             if (SWLK > (SMC - 0.02)) SWLK = SMC - 0.02
             if (SWLK < 0.0) SWLK = 0.0

             ! MATHEMATICAL SOLUTION BOUNDS APPLIED.
             DSWL = abs(SWLK - SWL)
             ! IF MORE THAN 10 ITERATIONS, USE EXPLICIT METHOD (CK=0 APPROX.)
             ! WHEN DSWL LESS OR EQ. ERROR, NO MORE ITERATIONS REQUIRED.
             SWL = SWLK
             if (DSWL <= ERROR) then
                KCOUNT = KCOUNT + 1
             end if
             ! END OF ITERATIONS
             ! BOUNDS APPLIED WITHIN DO-BLOCK ARE VALID FOR PHYSICAL SOLUTION.
          end do
          FREE = SMC - SWL
       end if ! END OPTION 1

       ! OPTION 2: EXPLICIT SOLUTION FOR FLERCHINGER EQ. i.e. CK=0
       ! IN KOREN ET AL., JGR, 1999, EQN 17
       ! APPLY PHYSICAL BOUNDS TO FLERCHINGER SOLUTION
       if (KCOUNT == 0) then
          write(message, '("Flerchinger used in NEW version. Iterations=", I6)') NLOG
          call wrf_message(trim(message))
          FK = (((hfus / (GRAV * (-LK_PSISAT(sltyp)))) * &
               & ((TKELV - TFRZ)/ TKELV)) ** (-1.0 / BX)) * LK_SMCMAX(sltyp)
          if (FK < 0.02) FK = 0.02
          FREE = min(FK, SMC)
       end if ! END OF OPTION 2
    end if

  end subroutine frh2o


  subroutine water(slptyp, sltyp, lutyp, dt, nsnow, nsoil, IMELT,UU, & !in
       VV     ,FCEV   ,FCTR   ,QPRECC ,QPRECL ,elai   , & !in
       esai   ,SFCTMP ,QVAP   ,QDEW   ,ZSOIL  ,BTRANI , & !in
       FICEOLD,PONDING,TG     ,IST    ,fveg   ,ILOC   ,JLOC, & !in
       LATHEAV , LATHEAG , frozen_canopy,frozen_ground,                       & !in  MB
       ISNOW  ,CANLIQ ,CANICE ,TV     ,snowh  ,sneqv  , & !inout
       SNICE  ,SNLIQ  ,STC    ,ZSNSO  ,soilwat   ,SMC    , & !inout
       soilice   ,ZWT    ,WA     ,WT     ,DZSNSO ,WSLAKE , & !inout
       CMC    ,ECAN   ,ETRAN  ,FWET   ,runsrf ,runsub , & !out
       QIN    ,QDIS   ,QSNOW  ,PONDING1       ,PONDING2,&
       QSNBOT,FPICE)  !out
    !
    ! Code history:
    ! Initial code: Guo-Yue Niu, Oct. 2007
    !
    use noahmp_veg_param, only: ISURBAN
    use noahmp_veg_param, only: LK_NROOT
    implicit none
    !
    ! input
    integer,                         intent(in)    :: ILOC    !grid index
    integer,                         intent(in)    :: JLOC    !grid index
    integer, intent(in) :: slptyp
    integer, intent(in) :: sltyp
    integer,                         intent(in)    :: lutyp  !vegetation type
    integer,                         intent(in)    :: NSNOW   !maximum no. of snow layers
    integer,                         intent(in)    :: IST     !surface type 1-soil; 2-lake
    integer,                         intent(in)    :: NSOIL   !no. of soil layers
    integer, dimension(-NSNOW+1:0) , intent(in)    :: IMELT   !melting state index [1-melt; 2-freeze]
    real,                            intent(in)    :: DT      !main time step (s)
    real,                            intent(in)    :: UU      !u-direction wind speed [m/s]
    real,                            intent(in)    :: VV      !v-direction wind speed [m/s]
    real,                            intent(in)    :: FCEV    !canopy evaporation (w/m2) [+ to atm ]
    real,                            intent(in)    :: FCTR    !transpiration (w/m2) [+ to atm]
    real,                            intent(in)    :: QPRECC  !convective precipitation (mm/s)
    real,                            intent(in)    :: QPRECL  !large-scale precipitation (mm/s)
    real,                            intent(in)    :: elai    !leaf area index, after burying by snow
    real,                            intent(in)    :: esai    !stem area index, after burying by snow
    real,                            intent(in)    :: SFCTMP  !surface air temperature [k]
    real,                            intent(in)    :: QVAP    !soil surface evaporation rate[mm/s]
    real,                            intent(in)    :: QDEW    !soil surface dew rate[mm/s]
    real, dimension(       1:NSOIL), intent(in)    :: ZSOIL   !depth of layer-bottom from soil surface
    real, dimension(       1:NSOIL), intent(in)    :: BTRANI  !soil water stress factor (0 to 1)
    real, dimension(-NSNOW+1:    0), intent(in)    :: FICEOLD !ice fraction at last timestep
    !  REAL                           , INTENT(in)    :: PONDING ![mm]
    real                           , intent(in)    :: TG      !ground temperature (k)
    real                           , intent(in)    :: fveg    !greeness vegetation fraction (-)

    ! input/output
    integer,                         intent(inout) :: ISNOW   !actual no. of snow layers
    real,                            intent(inout) :: CANLIQ  !intercepted liquid water (mm)
    real,                            intent(inout) :: CANICE  !intercepted ice mass (mm)
    real,                            intent(inout) :: TV      !vegetation temperature (k)
    real,                            intent(inout) :: snowh   !snow height [m]
    real,                            intent(inout) :: sneqv   !snow water eqv. [mm]
    real, dimension(-NSNOW+1:    0), intent(inout) :: SNICE   !snow layer ice [mm]
    real, dimension(-NSNOW+1:    0), intent(inout) :: SNLIQ   !snow layer liquid water [mm]
    real, dimension(-NSNOW+1:NSOIL), intent(inout) :: STC     !snow/soil layer temperature [k]
    real, dimension(-NSNOW+1:NSOIL), intent(inout) :: ZSNSO   !depth of snow/soil layer-bottom
    real, dimension(-NSNOW+1:NSOIL), intent(inout) :: DZSNSO  !snow/soil layer thickness [m]
    real, dimension(       1:NSOIL), intent(inout) :: soilwat    !soil liquid water content [m3/m3]
    real, dimension(       1:NSOIL), intent(inout) :: soilice    !soil ice content [m3/m3]
    real, dimension(       1:NSOIL), intent(inout) :: SMC     !total soil water content [m3/m3]
    real,                            intent(inout) :: ZWT     !the depth to water table [m]
    real,                            intent(inout) :: WA      !water storage in aquifer [mm]
    real,                            intent(inout) :: WT      !water storage in aquifer
    !+ stuarated soil [mm]
    real,                            intent(inout) :: WSLAKE  !water storage in lake (can be -) (mm)
    real,                            intent(inout) :: PONDING ![mm]

    ! output
    real,                            intent(out)   :: CMC     !intercepted water per ground area (mm)
    real,                            intent(out)   :: ECAN    !evap of intercepted water (mm/s) [+]
    real,                            intent(out)   :: ETRAN   !transpiration rate (mm/s) [+]
    real,                            intent(out)   :: FWET    !wetted/snowed fraction of canopy (-)
    real,                            intent(out)   :: runsrf  !surface runoff [mm/s]
    real,                            intent(out)   :: runsub  !baseflow (sturation excess) [mm/s]
    real,                            intent(out)   :: QIN     !groundwater recharge [mm/s]
    real,                            intent(out)   :: QDIS    !groundwater discharge [mm/s]
    real,                            intent(out)   :: QSNOW   !snow at ground srf (mm/s) [+]
    real,                            intent(out)   :: PONDING1
    real,                            intent(out)   :: PONDING2
    real,                            intent(out)   :: QSNBOT  !melting water out of snow bottom [mm/s]
    real,                            intent(out)   :: FPICE   !snow fraction in precipitation
    real,                             intent(in)   :: LATHEAV !latent heat vap./sublimation (j/kg)
    real,                             intent(in)   :: LATHEAG !latent heat vap./sublimation (j/kg)
    logical,                          intent(in)   :: FROZEN_GROUND ! used to define latent heat pathway
    logical,                          intent(in)   :: FROZEN_CANOPY ! used to define latent heat pathway

    ! local
    integer                                        :: IZ
    real                                           :: qinsrf  !water input on soil surface [m/s]
    real                                           :: QRAIN   !rain at ground srf (mm) [+]
    real                                           :: QSEVA   !soil surface evap rate [mm/s]
    real                                           :: QSDEW   !soil surface dew rate [mm/s]
    real                                           :: QSNFRO  !snow surface frost rate[mm/s]
    real                                           :: QSNSUB  !snow surface sublimation rate [mm/s]
    real                                           :: snowhin !snow depth increasing rate (m/s)
    real, dimension(       1:NSOIL)                :: ETRANI  !transpiration rate (mm/s) [+]
    real, dimension(       1:NSOIL)                :: WCND    !hydraulic conductivity (m/s)
    real                                           :: QDRAIN  !soil-bottom free drainage [mm/s]
    real                                           :: SNOFLOW !glacier flow [mm/s]
    real                                           :: FCRMAX !maximum of FCR (-)

    real, parameter ::  WSLMAX = 5000.      !maximum lake water storage (mm)

    ! initialize

    ETRANI(1:NSOIL) = 0.0
    SNOFLOW         = 0.0
    runsub          = 0.0
    qinsrf          = 0.0

    ! canopy-intercepted snowfall/rainfall, drips, and throughfall
    call canwater(lutyp, dt, SFCTMP, UU, VV, & !in
         FCEV, FCTR, QPRECC, QPRECL, elai, & !in
         esai, IST, TG, fveg, ILOC, JLOC, & !in
         FROZEN_CANOPY, & !in
         CANLIQ, CANICE, TV, & !inout
         CMC, ECAN, ETRAN, QRAIN, QSNOW, & !out
         snowhin, FWET, FPICE) !out

    ! sublimation, frost, evaporation, and dew

    QSNSUB = 0.0
    if (sneqv > 0.0) then
       QSNSUB = min(QVAP, sneqv / dt)
    end if
    QSEVA = QVAP - QSNSUB

    QSNFRO = 0.0
    if (sneqv > 0.0) then
       QSNFRO = QDEW
    end if
    QSDEW = QDEW - QSNFRO

    call snowwater(dt, nsnow, nsoil, zsoil, IMELT, & !in
         &          SFCTMP ,snowhin,QSNOW  ,QSNFRO ,QSNSUB , & !in
         &          QRAIN  ,FICEOLD,ILOC   ,JLOC   ,         & !in
         &          ISNOW  ,snowh  ,sneqv  ,SNICE  ,SNLIQ  , & !inout
         &          soilwat   ,soilice   ,STC    ,ZSNSO  ,DZSNSO , & !inout
         &          QSNBOT ,SNOFLOW,PONDING1       ,PONDING2)  !out

    if (FROZEN_GROUND) then
       soilice(1) =  soilice(1) + (QSDEW - QSEVA) * DT / (DZSNSO(1) * 1000.0)
       QSDEW = 0.0
       QSEVA = 0.0
       if (soilice(1) < 0.0) then
          soilwat(1) = soilwat(1) + soilice(1)
          soilice(1) = 0.0
       end if
    end if

    ! convert units (mm/s -> m/s)

    !PONDING: melting water from snow when there is no layer
    qinsrf = (PONDING + PONDING1 + PONDING2) / DT * 0.001
    !    qinsrf = PONDING/DT * 0.001

    if (ISNOW == 0) then
       qinsrf = qinsrf + (QSNBOT + QSDEW + QRAIN) * 0.001
    else
       qinsrf = qinsrf + (QSNBOT + QSDEW) * 0.001
    end if

    QSEVA  = QSEVA * 0.001

    do IZ = 1, LK_NROOT(lutyp)
       ETRANI(IZ) = ETRAN * BTRANI(IZ) * 0.001
    end do

    ! lake/soil water balances

    if (IST == 2) then                                        ! lake
       runsrf = 0.
       if (WSLAKE >= WSLMAX) runsrf = qinsrf * 1000.0             !mm/s
       WSLAKE = WSLAKE + (qinsrf - QSEVA) * 1000.0 * DT - runsrf * DT !mm
    else                                                      ! soil
       call soilh2o(sltyp, lutyp, &
            &  dt, nsoil, nsnow, zsoil, dzsnso, slptyp, & !in
            &  qinsrf, QSEVA, ETRANI, soilice, & !in
            &  soilwat, SMC, ZWT, & !inout
            &  runsrf ,QDRAIN ,runsub ,WCND   ,FCRMAX )   !out

       if (opt_run == 1) then
          call groundwater(sltyp, dt, nsnow, nsoil, zsoil, soilice, & !in
               STC    ,WCND   ,FCRMAX ,ILOC   ,JLOC   , & !in
               soilwat   ,ZWT    ,WA     ,WT     ,         & !inout
               QIN    ,QDIS   )                           !out
          runsub       = QDIS          !mm/s
       end if

       if (opt_run == 3 .or. opt_run == 4) then
          runsub = runsub + QDRAIN        !mm/s
       end if

       do IZ = 1,NSOIL
          SMC(IZ) = soilwat(IZ) + soilice(IZ)
       end do
    end if

    runsub = runsub + SNOFLOW         !mm/s

  end subroutine water


  subroutine canwater(lutyp, dt, SFCTMP, UU, VV, & !in
       FCEV,FCTR, QPRECC, QPRECL, elai, & !in
       esai, IST, TG, fveg, ILOC, JLOC, & !in
       FROZEN_CANOPY, & !in
       CANLIQ, CANICE, TV, & !inout
       CMC, ECAN, ETRAN, QRAIN, QSNOW, & !out
       snowhin, FWET, FPICE)                           !out
    ! canopy hydrology
    use noahmp_veg_param, only: LK_CANWMXP
    implicit none
    ! input
    integer,intent(in)  :: ILOC    !grid index
    integer,intent(in)  :: JLOC    !grid index
    integer,intent(in)  :: lutyp  !vegetation type
    real,   intent(in)  :: dt      !main time step (s)
    real,   intent(in)  :: SFCTMP  !air temperature (k)
    real,   intent(in)  :: UU      !u-direction wind speed [m/s]
    real,   intent(in)  :: VV      !v-direction wind speed [m/s]
    real,   intent(in)  :: FCEV    !canopy evaporation (w/m2) [+ = to atm]
    real,   intent(in)  :: FCTR    !transpiration (w/m2) [+ = to atm]
    real,   intent(in)  :: QPRECC  !convective precipitation (mm/s)
    real,   intent(in)  :: QPRECL  !large-scale precipitation (mm/s)
    real,   intent(in)  :: elai    !leaf area index, after burying by snow
    real,   intent(in)  :: esai    !stem area index, after burying by snow
    integer,intent(in)  :: IST     !surface type 1-soil; 2-lake
    real,   intent(in)  :: TG      !ground temperature (k)
    real,   intent(in)  :: fveg    !greeness vegetation fraction (-)
    logical                           , intent(in)   :: FROZEN_CANOPY ! used to define latent heat pathway

    ! input & output
    real, intent(inout) :: CANLIQ  !intercepted liquid water (mm)
    real, intent(inout) :: CANICE  !intercepted ice mass (mm)
    real, intent(inout) :: TV      !vegetation temperature (k)

    ! output
    real, intent(out)   :: CMC     !intercepted water (mm)
    real, intent(out)   :: ECAN    !evaporation of intercepted water (mm/s) [+]
    real, intent(out)   :: ETRAN   !transpiration rate (mm/s) [+]
    real, intent(out)   :: QRAIN   !rain at ground srf (mm/s) [+]
    real, intent(out)   :: QSNOW   !snow at ground srf (mm/s) [+]
    real, intent(out)   :: snowhin !snow depth increasing rate (m/s)
    real, intent(out)   :: FWET    !wetted or snowed fraction of the canopy (-)
    real, intent(out)   :: FPICE   !snow fraction in precipitation

    ! locals
    real                :: MAXSNO  !canopy capacity for snow interception (mm)
    real                :: MAXLIQ  !canopy capacity for rain interception (mm)
    real                :: FP      !fraction of the gridcell that receives precipitation
    real                :: BDFALL  !bulk density of snowfall (kg/m3)
    real                :: QINTR   !interception rate for rain (mm/s)
    real                :: QDRIPR  !drip rate for rain (mm/s)
    real                :: QTHROR  !throughfall for rain (mm/s)
    real                :: QINTS   !interception (loading) rate for snowfall (mm/s)
    real                :: QDRIPS  !drip (unloading) rate for intercepted snow (mm/s)
    real                :: QTHROS  !throughfall of snowfall (mm/s)
    real                :: QEVAC   !evaporation rate (mm/s)
    real                :: QDEWC   !dew rate (mm/s)
    real                :: QFROC   !frost rate (mm/s)
    real                :: QSUBC   !sublimation rate (mm/s)
    real                :: FT      !temperature factor for unloading rate
    real                :: FV      !wind factor for unloading rate
    real                :: QMELTC  !melting rate of canopy snow (mm/s)
    real                :: QFRZC   !refreezing rate of canopy liquid water (mm/s)
    real                :: RAIN    !rainfall (mm/s)
    real                :: SNOW    !snowfall (mm/s)
    real                :: CANMAS  !total canopy mass (kg/m2)

    ! initialization
    FP      = 0.0
    RAIN    = 0.0
    SNOW    = 0.0
    QINTR   = 0.
    QDRIPR  = 0.
    QTHROR  = 0.
    QINTR   = 0.
    QINTS   = 0.
    QDRIPS  = 0.0
    QTHROS  = 0.
    QRAIN   = 0.0
    QSNOW   = 0.0
    snowhin = 0.0
    ECAN    = 0.0

    ! partition precipitation into rain and snow.

    ! Jordan (1991)
    if (opt_snf == 1) then
       if (SFCTMP > TFRZ + 2.5) then
          FPICE = 0.0
       else
          if (SFCTMP <= TFRZ + 0.5) then
             FPICE = 1.0
          else if (SFCTMP <= TFRZ + 2.0) then
             FPICE = 1.0 - (-54.632 + 0.2 * SFCTMP)
          else
             FPICE = 0.6
          end if
       end if
    end if

    if (opt_snf == 2) then
       if (SFCTMP >= TFRZ + 2.2) then
          FPICE = 0.
       else
          FPICE = 1.0
       end if
    end if

    if (opt_snf == 3) then
       if (SFCTMP >= TFRZ) then
          FPICE = 0.0
       else
          FPICE = 1.0
       end if
    end if

    ! Hedstrom NR and JW Pomeroy (1998), Hydrol. Processes, 12, 1611-1625
    ! fresh snow density

    BDFALL = min(120.0, 67.92 + 51.25 * exp((SFCTMP - TFRZ) / 2.59))   ! Barlage: change to MIN in v3.6

    RAIN   = (QPRECC + QPRECL) * (1.0 - FPICE)
    SNOW   = (QPRECC + QPRECL) * FPICE

    ! fractional area that receives precipitation (see, Niu et al. 2005)

    if (QPRECC + QPRECL > 0.0) &
         FP = (QPRECC + QPRECL) / (10.0 * QPRECC + QPRECL)

    ! liquid water

    ! maximum canopy water
    MAXLIQ =  LK_CANWMXP(lutyp) * (elai + esai)

    ! average interception and throughfall

    if ((elai + esai) > 0.0) then
       QINTR  = fveg * RAIN * FP  ! interception capability
       QINTR  = min(QINTR, (MAXLIQ - CANLIQ) / DT * (1.0 - exp(-RAIN * DT / MAXLIQ)))
       QINTR  = max(QINTR, 0.0)
       QDRIPR = fveg * RAIN - QINTR
       QTHROR = (1.0 - fveg) * RAIN
    else
       QINTR  = 0.0
       QDRIPR = 0.0
       QTHROR = RAIN
    end if

    ! evaporation, transpiration, and dew
    if (.not. FROZEN_CANOPY) then             ! Barlage: change to frozen_canopy
       ETRAN = max(FCTR/HVAP, 0.0)
       QEVAC = max(FCEV/HVAP, 0.0)
       QDEWC = abs(min(FCEV/HVAP, 0.0))
       QSUBC = 0.0
       QFROC = 0.0
    else
       ETRAN = max( FCTR/HSUB, 0.0)
       QEVAC = 0.0
       QDEWC = 0.0
       QSUBC = max(FCEV/HSUB, 0.0)
       QFROC = abs(min(FCEV/HSUB, 0.0))
    end if

    ! canopy water balance. for convenience allow dew to bring CANLIQ above
    ! maxh2o or else would have to re-adjust drip

    QEVAC = min(CANLIQ / DT, QEVAC)
    CANLIQ = max(0.0, CANLIQ + (QINTR + QDEWC - QEVAC) * DT)
    if (CANLIQ <= 1.0E-6) CANLIQ = 0.0

    ! canopy ice
    MAXSNO = 6.6 * (0.27 + 46.0 / BDFALL) * (elai + esai)

    if ((elai + esai) > 0.0) then
       QINTS = fveg * SNOW * FP
       QINTS = min(QINTS, (MAXSNO - CANICE) / DT * (1.0 - exp(-SNOW * DT / MAXSNO)))
       QINTS = max(QINTS, 0.0)
       FT = max(0.0, (TV - 270.15) / 1.87E5)
       FV = sqrt(UU * UU + VV * VV) / 1.56E5
       QDRIPS = max(0.0, CANICE) * (FV + FT)
       QTHROS = (1.0 - fveg) * SNOW + (fveg * SNOW - QINTS)
    else
       QINTS  = 0.0
       QDRIPS = 0.0
       QTHROS = SNOW
    end if

    QSUBC = min(CANICE / DT, QSUBC)
    CANICE= max(0.0, CANICE + (QINTS - QDRIPS) * DT + (QFROC - QSUBC) * DT)
    if (CANICE <= 1.0E-6) CANICE = 0.0

    ! wetted fraction of canopy

    if (CANICE > 0.0) then
       FWET = max(0.0, CANICE) / max(MAXSNO, 1.0E-06)
    else
       FWET = max(0.0, CANLIQ) / max(MAXLIQ, 1.0E-06)
    end if
    FWET = min(FWET, 1.0) ** 0.667

    ! phase change

    QMELTC = 0.0
    QFRZC = 0.0

    if (CANICE > 1.0E-6 .and. TV > TFRZ) then
       QMELTC = min(CANICE / DT, (TV - TFRZ) * CICE * CANICE / DENICE / (DT * HFUS))
       CANICE = max(0.0, CANICE - QMELTC * DT)
       CANLIQ = max(0.0, CANLIQ + QMELTC * DT)
       TV     = FWET * TFRZ + (1.0 - FWET) * TV
    end if

    if (CANLIQ > 1.0E-6 .and. TV < TFRZ) then
       QFRZC  = min(CANLIQ / DT, (TFRZ - TV) * CWAT * CANLIQ / DENWAT / (DT * HFUS))
       CANLIQ = max(0.0, CANLIQ - QFRZC * DT)
       CANICE = max(0.0, CANICE + QFRZC * DT)
       TV     = FWET * TFRZ + (1.0 - FWET) * TV
    end if

    ! total canopy water

    CMC = CANLIQ + CANICE

    ! total canopy evaporation

    ECAN = QEVAC + QSUBC - QDEWC - QFROC

    ! rain or snow on the ground

    QRAIN   = QDRIPR + QTHROR
    QSNOW   = QDRIPS + QTHROS
    snowhin = QSNOW / BDFALL


    if (IST == 2 .and. TG > TFRZ) then
       QSNOW   = 0.0
       snowhin = 0.0
    end if

  end subroutine canwater


  subroutine snowwater(dt, nsnow, nsoil, zsoil, IMELT, & !in
       SFCTMP ,snowhin,QSNOW  ,QSNFRO ,QSNSUB , & !in
       QRAIN  ,FICEOLD,ILOC   ,JLOC   ,         & !in
       ISNOW  ,snowh  ,sneqv  ,SNICE  ,SNLIQ  , & !inout
       soilwat   ,soilice   ,STC    ,ZSNSO  ,DZSNSO , & !inout
       QSNBOT ,SNOFLOW,PONDING1       ,PONDING2)  !out
    implicit none
    ! input
    integer,                         intent(in)    :: ILOC   !grid index
    integer,                         intent(in)    :: JLOC   !grid index
    integer,                         intent(in)    :: NSNOW  !maximum no. of snow layers
    integer,                         intent(in)    :: NSOIL  !no. of soil layers
    integer, dimension(-NSNOW+1:0) , intent(in)    :: IMELT  !melting state index [0-no melt;1-melt]
    real,                            intent(in)    :: DT     !time step (s)
    real, dimension(       1:NSOIL), intent(in)    :: ZSOIL  !depth of layer-bottom from soil surface
    real,                            intent(in)    :: SFCTMP !surface air temperature [k]
    real,                            intent(in)    :: snowhin!snow depth increasing rate (m/s)
    real,                            intent(in)    :: QSNOW  !snow at ground srf (mm/s) [+]
    real,                            intent(in)    :: QSNFRO !snow surface frost rate[mm/s]
    real,                            intent(in)    :: QSNSUB !snow surface sublimation rate[mm/s]
    real,                            intent(in)    :: QRAIN  !snow surface rain rate[mm/s]
    real, dimension(-NSNOW+1:0)    , intent(in)    :: FICEOLD!ice fraction at last timestep

    ! input & output
    integer,                         intent(inout) :: ISNOW  !actual no. of snow layers
    real,                            intent(inout) :: snowh  !snow height [m]
    real,                            intent(inout) :: sneqv  !snow water eqv. [mm]
    real, dimension(-NSNOW+1:    0), intent(inout) :: SNICE  !snow layer ice [mm]
    real, dimension(-NSNOW+1:    0), intent(inout) :: SNLIQ  !snow layer liquid water [mm]
    real, dimension(       1:NSOIL), intent(inout) :: soilwat   !soil liquid moisture (m3/m3)
    real, dimension(       1:NSOIL), intent(inout) :: soilice   !soil ice moisture (m3/m3)
    real, dimension(-NSNOW+1:NSOIL), intent(inout) :: STC    !snow layer temperature [k]
    real, dimension(-NSNOW+1:NSOIL), intent(inout) :: ZSNSO  !depth of snow/soil layer-bottom
    real, dimension(-NSNOW+1:NSOIL), intent(inout) :: DZSNSO !snow/soil layer thickness [m]

    ! output
    real,                              intent(out) :: QSNBOT !melting water out of snow bottom [mm/s]
    real,                              intent(out) :: SNOFLOW!glacier flow [mm]
    real,                              intent(out) :: PONDING1
    real,                              intent(out) :: PONDING2

    ! local
    integer :: IZ,i
    real    :: BDSNOW  !bulk density of snow (kg/m3)
    SNOFLOW = 0.0
    PONDING1 = 0.0
    PONDING2 = 0.0

    call SNOWFALL(NSOIL  ,NSNOW  ,DT     ,QSNOW  ,snowhin, & !in
         SFCTMP ,ILOC   ,JLOC   ,                 & !in
         ISNOW  ,snowh  ,DZSNSO ,STC    ,SNICE  , & !inout
         SNLIQ  ,sneqv  )                           !inout

    ! MB: do each if block separately

    if (ISNOW < 0) &        ! when multi-layer
         call  COMPACT(NSNOW  ,NSOIL  ,DT     ,STC    ,SNICE  , & !in
         SNLIQ  ,ZSOIL  ,IMELT  ,FICEOLD,ILOC   , JLOC ,& !in
         ISNOW  ,DZSNSO ,ZSNSO  )                   !inout

    if (ISNOW < 0) &        !when multi-layer
         call  COMBINE(NSNOW  ,NSOIL  ,ILOC   ,JLOC   ,         & !in
         ISNOW  ,soilwat   ,STC    ,SNICE  ,SNLIQ  , & !inout
         DZSNSO ,soilice   ,snowh  ,sneqv  ,         & !inout
         PONDING1       ,PONDING2)                  !out

    if (ISNOW < 0) &        !when multi-layer
         call   DIVIDE(NSNOW  ,NSOIL  ,                         & !in
         ISNOW  ,STC    ,SNICE  ,SNLIQ  ,DZSNSO )   !inout

    call  snowh2o(NSNOW  ,NSOIL  ,DT     ,QSNFRO ,QSNSUB , & !in
         QRAIN  ,ILOC   ,JLOC   ,                 & !in
         ISNOW  ,DZSNSO ,snowh  ,sneqv  ,SNICE  , & !inout
         SNLIQ  ,soilwat   ,soilice   ,STC    ,         & !inout
         QSNBOT ,PONDING1       ,PONDING2)           !out

    !set empty snow layers to zero

    do iz = -nsnow+1, isnow
       snice(iz) = 0.0
       snliq(iz) = 0.0
       stc(iz)   = 0.0
       dzsnso(iz)= 0.0
       zsnso(iz) = 0.0
    end do

    !to obtain equilibrium state of snow in glacier region

    if (sneqv > 2000.0) then   ! 2000 mm -> maximum water depth
       BDSNOW      = SNICE(0) / DZSNSO(0)
       SNOFLOW     = (sneqv - 2000.0)
       SNICE(0)    = SNICE(0)  - SNOFLOW
       DZSNSO(0)   = DZSNSO(0) - SNOFLOW / BDSNOW
       SNOFLOW     = SNOFLOW / DT
    end if

    ! sum up snow mass for layered snow

    if (ISNOW < 0) then  ! MB: only do for multi-layer
       sneqv = 0.0
       do IZ = ISNOW+1, 0
          sneqv = sneqv + SNICE(IZ) + SNLIQ(IZ)
       end do
    end if

    ! Reset ZSNSO and layer thinkness DZSNSO

    do IZ = ISNOW+1, 0
       DZSNSO(IZ) = -DZSNSO(IZ)
    end do

    DZSNSO(1) = ZSOIL(1)
    do IZ = 2, NSOIL
       DZSNSO(IZ) = (ZSOIL(IZ) - ZSOIL(IZ-1))
    end do

    ZSNSO(ISNOW+1) = DZSNSO(ISNOW+1)
    do IZ = ISNOW+2, NSOIL
       ZSNSO(IZ) = ZSNSO(IZ-1) + DZSNSO(IZ)
    end do

    do IZ = ISNOW+1, NSOIL
       DZSNSO(IZ) = -DZSNSO(IZ)
    end do

  end subroutine snowwater


  subroutine snowfall(NSOIL  ,NSNOW  ,DT     ,QSNOW  ,snowhin , & !in
       SFCTMP ,ILOC   ,JLOC   ,                  & !in
       ISNOW  ,snowh  ,DZSNSO ,STC    ,SNICE   , & !inout
       SNLIQ  ,sneqv  )                            !inout
    ! snow depth and density to account for the new snowfall.
    ! new values of snow depth & density returned.
    implicit none
    ! input
    integer,                            intent(in) :: ILOC   !grid index
    integer,                            intent(in) :: JLOC   !grid index
    integer,                            intent(in) :: NSOIL  !no. of soil layers
    integer,                            intent(in) :: NSNOW  !maximum no. of snow layers
    real,                               intent(in) :: DT     !main time step (s)
    real,                               intent(in) :: QSNOW  !snow at ground srf (mm/s) [+]
    real,                               intent(in) :: snowhin!snow depth increasing rate (m/s)
    real,                               intent(in) :: SFCTMP !surface air temperature [k]

    ! input and output
    integer,                         intent(inout) :: ISNOW  !actual no. of snow layers
    real,                            intent(inout) :: snowh  !snow depth [m]
    real,                            intent(inout) :: sneqv  !swow water equivalent [m]
    real, dimension(-NSNOW+1:NSOIL), intent(inout) :: DZSNSO !thickness of snow/soil layers (m)
    real, dimension(-NSNOW+1:NSOIL), intent(inout) :: STC    !snow layer temperature [k]
    real, dimension(-NSNOW+1:    0), intent(inout) :: SNICE  !snow layer ice [mm]
    real, dimension(-NSNOW+1:    0), intent(inout) :: SNLIQ  !snow layer liquid water [mm]

    ! local
    integer :: NEWNODE            ! 0-no new layers, 1-creating new layers
    NEWNODE  = 0

    ! shallow snow / no layer

    if (ISNOW == 0 .and. QSNOW > 0.0)  then
       snowh = snowh + snowhin * DT
       sneqv = sneqv + QSNOW * DT
    end if

    ! creating a new layer

    if (ISNOW == 0  .and. QSNOW > 0.0 .and. snowh >= 0.025) then !MB: change limit
       !    if (ISNOW == 0  .AND. QSNOW>0. .AND. snowh >= 0.05) THEN
       ISNOW    = -1
       NEWNODE  =  1
       DZSNSO(0)= snowh
       snowh    = 0.0
       STC(0)   = min(273.16, SFCTMP)   ! temporary setup
       SNICE(0) = sneqv
       SNLIQ(0) = 0.0
    end if

    ! snow with layers

    if (ISNOW <  0 .and. NEWNODE == 0 .and. QSNOW > 0.0) then
       SNICE(ISNOW+1)  = SNICE(ISNOW+1)   + QSNOW   * DT
       DZSNSO(ISNOW+1) = DZSNSO(ISNOW+1)  + snowhin * DT
    end if
  end subroutine snowfall


  subroutine combine(NSNOW  ,NSOIL  ,ILOC   ,JLOC   ,         & !in
       ISNOW  ,soilwat   ,STC    ,SNICE  ,SNLIQ  , & !inout
       DZSNSO ,soilice   ,SNOWH  ,sneqv  ,         & !inout
       PONDING1       ,PONDING2)                  !out
    implicit none
    ! input

    integer, intent(in)     :: ILOC
    integer, intent(in)     :: JLOC
    integer, intent(in)     :: NSNOW                        !maximum no. of snow layers
    integer, intent(in)     :: NSOIL                        !no. of soil layers

    ! input and output

    integer,                         intent(inout) :: ISNOW !actual no. of snow layers
    real, dimension(       1:NSOIL), intent(inout) :: soilwat  !soil liquid moisture (m3/m3)
    real, dimension(       1:NSOIL), intent(inout) :: soilice  !soil ice moisture (m3/m3)
    real, dimension(-NSNOW+1:NSOIL), intent(inout) :: STC   !snow layer temperature [k]
    real, dimension(-NSNOW+1:    0), intent(inout) :: SNICE !snow layer ice [mm]
    real, dimension(-NSNOW+1:    0), intent(inout) :: SNLIQ !snow layer liquid water [mm]
    real, dimension(-NSNOW+1:NSOIL), intent(inout) :: DZSNSO!snow layer depth [m]
    real,                            intent(inout) :: sneqv !snow water equivalent [m]
    real,                            intent(inout) :: snowh !snow depth [m]
    real,                            intent(out) :: PONDING1
    real,                            intent(out) :: PONDING2

    ! local variables:

    integer :: I,J,K,L               ! node indices
    integer :: ISNOW_OLD             ! number of top snow layer
    integer :: MSSI                  ! node index
    integer :: NEIBOR                ! adjacent node selected for combination
    real    :: ZWICE                 ! total ice mass in snow
    real    :: ZWLIQ                 ! total liquid water in snow

    ! real(r4) :: DZMIN(3) = (/0.045, 0.05, 0.2/) !minimum of top snow layer
    real(r4) :: DZMIN(3) = (/0.025, 0.025, 0.1/) ! MB: change limit

    ISNOW_OLD = ISNOW

    do J = ISNOW_OLD+1, 0
       if (SNICE(J) <= 0.1) then
          if (J /= 0) then
             SNLIQ(J+1) = SNLIQ(J+1) + SNLIQ(J)
             SNICE(J+1) = SNICE(J+1) + SNICE(J)
          else
             if (ISNOW_OLD < -1) then    ! MB/KM: change to ISNOW
                SNLIQ(J-1) = SNLIQ(J-1) + SNLIQ(J)
                SNICE(J-1) = SNICE(J-1) + SNICE(J)
             else
                if (SNICE(J) >= 0.0) then
                   PONDING1 = SNLIQ(J)    ! ISNOW WILL GET SET TO ZERO BELOW; PONDING1 WILL GET
                   sneqv = SNICE(J)       ! ADDED TO PONDING FROM PHASECHANGE PONDING SHOULD BE
                   snowh = DZSNSO(J)      ! ZERO HERE BECAUSE IT WAS CALCULATED FOR THIN SNOW
                else   ! SNICE OVER-SUBLIMATED EARLIER
                   PONDING1 = SNLIQ(J) + SNICE(J)
                   if (PONDING1 < 0.0) then  ! IF SNICE AND SNLIQ SUBLIMATES REMOVE FROM SOIL
                      soilice(1) = max(0.0, soilice(1) + PONDING1 / (DZSNSO(1) * 1000.0))
                      PONDING1 = 0.0
                   end if
                   sneqv = 0.0
                   snowh = 0.0
                end if
                SNLIQ(J) = 0.0
                SNICE(J) = 0.0
                DZSNSO(J) = 0.0
             end if
             !                soilwat(1) = soilwat(1)+SNLIQ(J)/(DZSNSO(1)*1000.)
             !                soilice(1) = soilice(1)+SNICE(J)/(DZSNSO(1)*1000.)
          end if

          ! shift all elements above this down by one.
          if (J > ISNOW + 1 .and. ISNOW < -1) then
             do I = J, ISNOW+2, -1
                STC(I) = STC(I-1)
                SNLIQ(I) = SNLIQ(I-1)
                SNICE(I) = SNICE(I-1)
                DZSNSO(I)= DZSNSO(I-1)
             end do
          end if
          ISNOW = ISNOW + 1
       end if
    end do

    ! to conserve water in case of too large surface sublimation

    if (soilice(1) < 0.0) then
       soilwat(1) = soilwat(1) + soilice(1)
       soilice(1) = 0.0
    end if

    if (ISNOW == 0) return   ! MB: get out if no longer multi-layer

    sneqv  = 0.0
    snowh  = 0.0
    ZWICE  = 0.0
    ZWLIQ  = 0.0

    do J = ISNOW+1, 0
       sneqv = sneqv + SNICE(J) + SNLIQ(J)
       snowh = snowh + DZSNSO(J)
       ZWICE = ZWICE + SNICE(J)
       ZWLIQ = ZWLIQ + SNLIQ(J)
    end do

    ! check the snow depth - all snow gone
    ! the liquid water assumes ponding on soil surface.

    if (snowh < 0.025 .and. ISNOW < 0 ) then ! MB: change limit
       !       IF (SNOWH < 0.05 .AND. ISNOW < 0 ) THEN
       ISNOW  = 0
       sneqv = ZWICE
       PONDING2 = ZWLIQ           ! LIMIT OF ISNOW < 0 MEANS INPUT PONDING
       if (sneqv <= 0.0) snowh = 0.0 ! SHOULD BE ZERO; SEE ABOVE
    end if

    !       IF (snowh < 0.05 ) THEN
    !          ISNOW  = 0
    !          sneqv = ZWICE
    !          soilwat(1) = soilwat(1) + ZWLIQ / (DZSNSO(1) * 1000.)
    !          if (sneqv <= 0.) SNOWH = 0.
    !       END IF

    ! check the snow depth - snow layers combined

    if (ISNOW < -1) then

       ISNOW_OLD = ISNOW
       MSSI     = 1

       do I = ISNOW_OLD+1, 0
          if (DZSNSO(I) < DZMIN(MSSI)) then

             if (I == ISNOW + 1) then
                NEIBOR = I + 1
             else if (I == 0) then
                NEIBOR = I - 1
             else
                NEIBOR = I + 1
                if ((DZSNSO(I-1) + DZSNSO(I)) < (DZSNSO(I+1) + DZSNSO(I))) NEIBOR = I - 1
             end if

             ! Node l and j are combined and stored as node j.
             if (NEIBOR > I) then
                J = NEIBOR
                L = I
             else
                J = I
                L = NEIBOR
             end if

             call COMBO(DZSNSO(J), SNLIQ(J), SNICE(J), &
                  STC(J), DZSNSO(L), SNLIQ(L), SNICE(L), STC(L) )

             ! Now shift all elements above this down one.
             if (J-1 > ISNOW+1) then
                do K = J-1, ISNOW+2, -1
                   STC(K)   = STC(K-1)
                   SNICE(K) = SNICE(K-1)
                   SNLIQ(K) = SNLIQ(K-1)
                   DZSNSO(K) = DZSNSO(K-1)
                end do
             end if

             ! Decrease the number of snow layers
             ISNOW = ISNOW + 1
             if (ISNOW >= -1) exit
          else

             ! The layer thickness is greater than the prescribed minimum value
             MSSI = MSSI + 1

          end if
       end do

    end if

  end subroutine combine


  subroutine divide(NSNOW  ,NSOIL  ,                         & !in
       ISNOW  ,STC    ,SNICE  ,SNLIQ  ,DZSNSO  )  !inout
    implicit none
    ! input
    integer, intent(in)                            :: NSNOW !maximum no. of snow layers [ =3]
    integer, intent(in)                            :: NSOIL !no. of soil layers [ =4]

    ! input and output
    integer                        , intent(inout) :: ISNOW !actual no. of snow layers
    real, dimension(-NSNOW+1:NSOIL), intent(inout) :: STC   !snow layer temperature [k]
    real, dimension(-NSNOW+1:    0), intent(inout) :: SNICE !snow layer ice [mm]
    real, dimension(-NSNOW+1:    0), intent(inout) :: SNLIQ !snow layer liquid water [mm]
    real, dimension(-NSNOW+1:NSOIL), intent(inout) :: DZSNSO!snow layer depth [m]

    ! local variables:
    integer                                        :: J     !indices
    integer                                        :: MSNO  !number of layer (top) to MSNO (bot)
    real                                           :: DRR   !thickness of the combined [m]
    real, dimension(       1:NSNOW)                :: DZ    !snow layer thickness [m]
    real, dimension(       1:NSNOW)                :: SWICE !partial volume of ice [m3/m3]
    real, dimension(       1:NSNOW)                :: SWLIQ !partial volume of liquid water [m3/m3]
    real, dimension(       1:NSNOW)                :: TSNO  !node temperature [k]
    real                                           :: ZWICE !temporary
    real                                           :: ZWLIQ !temporary
    real                                           :: PROPOR!temporary
    real                                           :: DTDZ  !temporary

    do J = 1, NSNOW
       if (J <= abs(ISNOW)) then
          DZ(J)    = DZSNSO(J+ISNOW)
          SWICE(J) = SNICE(J+ISNOW)
          SWLIQ(J) = SNLIQ(J+ISNOW)
          TSNO(J)  = STC(J+ISNOW)
       end if
    end do

    MSNO = abs(ISNOW)

    if (MSNO == 1) then
       ! Specify a new snow layer
       if (DZ(1) > 0.05) then
          MSNO = 2
          DZ(1)    = DZ(1) / 2.0
          SWICE(1) = SWICE(1) / 2.0
          SWLIQ(1) = SWLIQ(1) / 2.0
          DZ(2)    = DZ(1)
          SWICE(2) = SWICE(1)
          SWLIQ(2) = SWLIQ(1)
          TSNO(2)  = TSNO(1)
       end if
    end if

    if (MSNO > 1) then
       if (DZ(1) > 0.05) then
          DRR      = DZ(1) - 0.05
          PROPOR   = DRR / DZ(1)
          ZWICE    = PROPOR * SWICE(1)
          ZWLIQ    = PROPOR * SWLIQ(1)
          PROPOR   = 0.05 / DZ(1)
          SWICE(1) = PROPOR * SWICE(1)
          SWLIQ(1) = PROPOR * SWLIQ(1)
          DZ(1)    = 0.05

          call COMBO(DZ(2), SWLIQ(2), SWICE(2), TSNO(2), DRR, &
               ZWLIQ, ZWICE, TSNO(1))

          ! subdivide a new layer
          if (MSNO <= 2 .and. DZ(2) > 0.20) then  ! MB: change limit
             !             IF (MSNO <= 2 .AND. DZ(2) > 0.10) THEN
             MSNO = 3
             DTDZ = (TSNO(1) - TSNO(2)) / ((DZ(1) + DZ(2))/2.)
             DZ(2)    = DZ(2) / 2.0
             SWICE(2) = SWICE(2) / 2.0
             SWLIQ(2) = SWLIQ(2) / 2.0
             DZ(3)    = DZ(2)
             SWICE(3) = SWICE(2)
             SWLIQ(3) = SWLIQ(2)
             TSNO(3) = TSNO(2) - DTDZ * DZ(2) / 2.0
             if (TSNO(3) >= TFRZ) then
                TSNO(3)  = TSNO(2)
             else
                TSNO(2) = TSNO(2) + DTDZ * DZ(2) / 2.0
             end if

          end if
       end if
    end if

    if (MSNO > 2) then
       if (DZ(2) > 0.2) then
          DRR = DZ(2) - 0.2
          PROPOR   = DRR / DZ(2)
          ZWICE    = PROPOR * SWICE(2)
          ZWLIQ    = PROPOR * SWLIQ(2)
          PROPOR   = 0.2 / DZ(2)
          SWICE(2) = PROPOR * SWICE(2)
          SWLIQ(2) = PROPOR * SWLIQ(2)
          DZ(2)    = 0.2
          call COMBO(DZ(3), SWLIQ(3), SWICE(3), TSNO(3), DRR, &
               ZWLIQ, ZWICE, TSNO(2))
       end if
    end if

    ISNOW = -MSNO

    do J = ISNOW+1, 0
       DZSNSO(J) = DZ(J-ISNOW)
       SNICE(J) = SWICE(J-ISNOW)
       SNLIQ(J) = SWLIQ(J-ISNOW)
       STC(J)   = TSNO(J-ISNOW)
    end do


    !    DO J = ISNOW+1,NSOIL
    !    WRITE(*,'(I5,7F10.3)') J, DZSNSO(J), SNICE(J), SNLIQ(J),STC(J)
    !    END DO

  end subroutine divide


  subroutine combo(DZ,  WLIQ,  WICE, T, DZ2, WLIQ2, WICE2, T2)
    implicit none
    ! input
    real, intent(in)    :: DZ2   !nodal thickness of 2 elements being combined [m]
    real, intent(in)    :: WLIQ2 !liquid water of element 2 [kg/m2]
    real, intent(in)    :: WICE2 !ice of element 2 [kg/m2]
    real, intent(in)    :: T2    !nodal temperature of element 2 [k]
    real, intent(inout) :: DZ    !nodal thickness of 1 elements being combined [m]
    real, intent(inout) :: WLIQ  !liquid water of element 1
    real, intent(inout) :: WICE  !ice of element 1 [kg/m2]
    real, intent(inout) :: T     !node temperature of element 1 [k]

    ! local
    real                :: DZC   !total thickness of nodes 1 and 2 (DZC=DZ+DZ2).
    real                :: WLIQC !combined liquid water [kg/m2]
    real                :: WICEC !combined ice [kg/m2]
    real                :: TC    !combined node temperature [k]
    real                :: H     !enthalpy of element 1 [J/m2]
    real                :: H2    !enthalpy of element 2 [J/m2]
    real                :: HC    !temporary

    DZC = DZ + DZ2
    WICEC = WICE + WICE2
    WLIQC = WLIQ + WLIQ2
    H = (CICE * WICE + CWAT * WLIQ) * (T - TFRZ) + HFUS * WLIQ
    H2= (CICE * WICE2 + CWAT * WLIQ2) * (T2 - TFRZ) + HFUS * WLIQ2

    HC = H + H2
    if (HC < 0.0) then
       TC = TFRZ + HC / (CICE * WICEC + CWAT * WLIQC)
    else if (HC <= HFUS * WLIQC) then
       TC = TFRZ
    else
       TC = TFRZ + (HC - HFUS * WLIQC) / (CICE * WICEC + CWAT * WLIQC)
    end if

    DZ = DZC
    WICE = WICEC
    WLIQ = WLIQC
    T = TC

  end subroutine combo


  subroutine compact(NSNOW, NSOIL, DT, STC, SNICE, & !in
       & SNLIQ, ZSOIL, IMELT, FICEOLD, ILOC, JLOC, & !in
       & ISNOW, DZSNSO, ZSNSO )                    !inout
    implicit none
    ! input
    integer,                         intent(in)    :: ILOC   !grid index
    integer,                         intent(in)    :: JLOC   !grid index
    integer,                         intent(in)    :: NSOIL  !no. of soil layers [ =4]
    integer,                         intent(in)    :: NSNOW  !maximum no. of snow layers [ =3]
    integer, dimension(-NSNOW+1:0) , intent(in)    :: IMELT  !melting state index [0-no melt;1-melt]
    real,                            intent(in)    :: DT     !time step (sec)
    real, dimension(-NSNOW+1:NSOIL), intent(in)    :: STC    !snow layer temperature [k]
    real, dimension(-NSNOW+1:    0), intent(in)    :: SNICE  !snow layer ice [mm]
    real, dimension(-NSNOW+1:    0), intent(in)    :: SNLIQ  !snow layer liquid water [mm]
    real, dimension(       1:NSOIL), intent(in)    :: ZSOIL  !depth of layer-bottom from soil srf
    real, dimension(-NSNOW+1:    0), intent(in)    :: FICEOLD!ice fraction at last timestep

    ! input and output
    integer,                         intent(inout) :: ISNOW  ! actual no. of snow layers
    real, dimension(-NSNOW+1:NSOIL), intent(inout) :: DZSNSO ! snow layer thickness [m]
    real, dimension(-NSNOW+1:NSOIL), intent(inout) :: ZSNSO  ! depth of snow/soil layer-bottom

    ! local
    real, parameter     :: C2 = 21.e-3   ![m3/kg] ! default 21.e-3
    real, parameter     :: C3 = 2.5e-6   ![1/s]
    real, parameter     :: C4 = 0.04     ![1/k]
    real, parameter     :: C5 = 2.0      !
    real, parameter     :: DM = 100.0    !upper Limit on destructive metamorphism compaction [kg/m3]
    real, parameter     :: ETA0 = 0.8e+6 !viscosity coefficient [kg-s/m2]
    !according to Anderson, it is between 0.52e6~1.38e6
    real :: BURDEN !pressure of overlying snow [kg/m2]
    real :: DDZ1   !rate of settling of snow pack due to destructive metamorphism.
    real :: DDZ2   !rate of compaction of snow pack due to overburden.
    real :: DDZ3   !rate of compaction of snow pack due to melt [1/s]
    real :: DEXPF  !EXPF=exp(-c4*(273.15-STC)).
    real :: TD     !STC - TFRZ [K]
    real :: PDZDTC !nodal rate of change in fractional-thickness due to compaction [fraction/s]
    real :: VOID   !void (1 - SNICE - SNLIQ)
    real :: WX     !water mass (ice + liquid) [kg/m2]
    real :: BI     !partial density of ice [kg/m3]
    real, dimension(-NSNOW+1:0) :: FICE   !fraction of ice at current time step

    integer  :: J

    BURDEN = 0.0

    do J = ISNOW+1, 0

       WX      = SNICE(J) + SNLIQ(J)
       FICE(J) = SNICE(J) / WX
       VOID    = 1. - (SNICE(J) / DENICE + SNLIQ(J) / DENWAT) / DZSNSO(J)

       ! Allow compaction only for non-saturated node and higher ice lens node.
       if (VOID > 0.001 .and. SNICE(J) > 0.1) then
          BI = SNICE(J) / DZSNSO(J)
          TD = max(0.0, TFRZ - STC(J))
          DEXPF = exp(-C4 * TD)

          ! Settling as a result of destructive metamorphism

          DDZ1 = -C3 * DEXPF

          if (BI > DM) DDZ1 = DDZ1 * exp(-46.0E-3 * (BI - DM))

          ! Liquid water term

          if (SNLIQ(J) > 0.01 * DZSNSO(J)) DDZ1 = DDZ1 * C5

          ! Compaction due to overburden

          DDZ2 = -(BURDEN + 0.5 * WX) * exp(-0.08 * TD - C2 * BI) / ETA0 ! 0.5*WX -> self-burden

          ! Compaction occurring during melt

          if (IMELT(J) == 1) then
             DDZ3 = max(0.0, (FICEOLD(J) - FICE(J)) / max(1.E-6, FICEOLD(J)))
             DDZ3 = - DDZ3 / DT           ! sometimes too large
          else
             DDZ3 = 0.0
          end if

          ! Time rate of fractional change in DZ (units of s-1)

          PDZDTC = (DDZ1 + DDZ2 + DDZ3) * DT
          PDZDTC = max(-0.5, PDZDTC)

          ! The change in DZ due to compaction

          DZSNSO(J) = DZSNSO(J) * (1.0 + PDZDTC)
       end if

       ! Pressure of overlying snow

       BURDEN = BURDEN + WX

    end do

  end subroutine compact


  subroutine snowh2o(NSNOW  ,NSOIL  ,DT     ,QSNFRO ,QSNSUB , & !in
       QRAIN  ,ILOC   ,JLOC   ,                 & !in
       ISNOW  ,DZSNSO ,snowh  ,sneqv  ,SNICE  , & !inout
       SNLIQ  ,soilwat   ,soilice   ,STC    ,         & !inout
       QSNBOT ,PONDING1       ,PONDING2)          !out
    ! Renew the mass of ice lens (SNICE) and liquid (SNLIQ) of the
    ! surface snow layer resulting from sublimation (frost) / evaporation (dew)
    use noahmp_gen_param, only: KK_SSI
    implicit none
    ! input
    integer,                         intent(in)    :: ILOC   !grid index
    integer,                         intent(in)    :: JLOC   !grid index
    integer,                         intent(in)    :: NSNOW  !maximum no. of snow layers[=3]
    integer,                         intent(in)    :: NSOIL  !No. of soil layers[=4]
    real,                            intent(in)    :: DT     !time step
    real,                            intent(in)    :: QSNFRO !snow surface frost rate[mm/s]
    real,                            intent(in)    :: QSNSUB !snow surface sublimation rate[mm/s]
    real,                            intent(in)    :: QRAIN  !snow surface rain rate[mm/s]

    ! output
    real,                            intent(out)   :: QSNBOT !melting water out of snow bottom [mm/s]

    ! input and output
    integer,                         intent(inout) :: ISNOW  !actual no. of snow layers
    real, dimension(-NSNOW+1:NSOIL), intent(inout) :: DZSNSO ! snow layer depth [m]
    real,                            intent(inout) :: snowh  !snow height [m]
    real,                            intent(inout) :: sneqv  !snow water eqv. [mm]
    real, dimension(-NSNOW+1:0),     intent(inout) :: SNICE  !snow layer ice [mm]
    real, dimension(-NSNOW+1:0),     intent(inout) :: SNLIQ  !snow layer liquid water [mm]
    real, dimension(       1:NSOIL), intent(inout) :: soilwat   !soil liquid moisture (m3/m3)
    real, dimension(       1:NSOIL), intent(inout) :: soilice   !soil ice moisture (m3/m3)
    real, dimension(-NSNOW+1:NSOIL), intent(inout) :: STC    !snow layer temperature [k]

    ! locals
    integer                     :: J         !do loop/array indices
    real                        :: QIN       !water flow into the element (mm/s)
    real                        :: QOUT      !water flow out of the element (mm/s)
    real                        :: WGDIF     !ice mass after minus sublimation
    real, dimension(-NSNOW+1:0) :: VOL_LIQ   !partial volume of liquid water in layer
    real, dimension(-NSNOW+1:0) :: VOL_ICE   !partial volume of ice lens in layer
    real, dimension(-NSNOW+1:0) :: epore     !effective porosity = porosity - VOL_ICE
    real :: PROPOR, TEMP
    real :: PONDING1, PONDING2

    !for the case when sneqv becomes '0' after 'COMBINE'

    if (sneqv == 0.0) then
       soilice(1) =  soilice(1) + (QSNFRO - QSNSUB) * DT / (DZSNSO(1) * 1000.0)  ! Barlage: soilwat->soilice v3.6
       if (soilice(1) < 0.0) then
          soilwat(1) = soilwat(1) + soilice(1)
          soilice(1) = 0.0
       end if
    end if

    ! for shallow snow without a layer
    ! snow surface sublimation may be larger than existing snow mass. To conserve water,
    ! excessive sublimation is used to reduce soil water. Smaller time steps would tend
    ! to aviod this problem.

    if (ISNOW == 0 .and. sneqv > 0.0) then
       TEMP   = sneqv
       sneqv  = sneqv - QSNSUB * DT + QSNFRO * DT
       PROPOR = sneqv / TEMP
       snowh  = max(0.0, PROPOR * snowh)

       if (sneqv < 0.0) then
          soilice(1) = soilice(1) + sneqv / (DZSNSO(1) * 1000.0)
          sneqv   = 0.0
          snowh   = 0.0
       end if
       if (soilice(1) < 0.0) then
          soilwat(1) = soilwat(1) + soilice(1)
          soilice(1) = 0.0
       end if
    end if

    if (snowh <= 1.0E-8 .or. sneqv <= 1.0E-6) then
       snowh = 0.0
       sneqv = 0.0
    end if

    ! for deep snow

    if (ISNOW < 0) then !KWM added this IF statement to prevent out-of-bounds array references
       WGDIF = SNICE(ISNOW+1) - QSNSUB * DT + QSNFRO * DT
       SNICE(ISNOW+1) = WGDIF
       if (WGDIF < 1.0E-6 .and. ISNOW < 0) then
          call  COMBINE(NSNOW  ,NSOIL  ,ILOC, JLOC   , & !in
               ISNOW  ,soilwat   ,STC    ,SNICE  ,SNLIQ  , & !inout
               DZSNSO ,soilice   ,snowh  ,sneqv  ,         & !inout
               PONDING1, PONDING2)                       !out
       end if
       !KWM:  Subroutine COMBINE can change ISNOW to make it 0 again?
       if (ISNOW < 0) then !KWM added this IF statement to prevent out-of-bounds array references
          SNLIQ(ISNOW+1) = SNLIQ(ISNOW+1) + QRAIN * DT
          SNLIQ(ISNOW+1) = max(0.0, SNLIQ(ISNOW+1))
       end if

    end if !KWM  -- Can the ENDIF be moved toward the end of the subroutine (Just set QSNBOT=0)?

    ! Porosity and partial volume

    !KWM Looks to me like loop index / IF test can be simplified.

    do J = -NSNOW+1, 0
       if (J >= ISNOW + 1) then
          VOL_ICE(J)      = min(1.0, SNICE(J) / (DZSNSO(J) * DENICE))
          epore(J)        = 1.0 - VOL_ICE(J)
          VOL_LIQ(J)      = min(epore(J), SNLIQ(J) / (DZSNSO(J) * DENWAT))
       end if
    end do

    QIN = 0.0
    QOUT = 0.0

    !KWM Looks to me like loop index / IF test can be simplified.

    do J = -NSNOW+1, 0
       if (J >= ISNOW + 1) then
          SNLIQ(J) = SNLIQ(J) + QIN
          if (J <= -1) then
             if (epore(J) < 0.05 .or. epore(J+1) < 0.05) then
                QOUT = 0.0
             else
                QOUT = max(0.0, (VOL_LIQ(J) - KK_SSI * epore(J)) * DZSNSO(J))
                QOUT = min(QOUT, (1.0 - VOL_ICE(J+1) - VOL_LIQ(J+1)) * DZSNSO(J+1))
             end if
          else
             QOUT = max(0.0, (VOL_LIQ(J) - KK_SSI * epore(J)) * DZSNSO(J))
          end if
          QOUT = QOUT * 1000.0
          SNLIQ(J) = SNLIQ(J) - QOUT
          QIN = QOUT
       end if
    end do

    ! Liquid water from snow bottom to soil

    QSNBOT = QOUT / DT           ! mm/s
  end subroutine snowh2o


  subroutine soilh2o(sltyp, lutyp, dt, nsoil, nsnow, zsoil, dzsnso, slptyp, & !in
       & qinsrf, QSEVA  ,ETRANI ,soilice   , & !in
       & soilwat, SMC, ZWT, & !inout
       & runsrf, QDRAIN, runsub, WCND, FCRMAX)   !out
    ! calculate surface runoff and soil moisture.
    use noahmp_gen_param, only: KK_TIMEAN
    use noahmp_gen_param, only: KK_FSATMAX
    use noahmp_soil_param, only: LK_SMCMAX
    use noahmp_veg_param, only: ISURBAN
    implicit none
    integer, intent(in) :: sltyp
    integer, intent(in) :: lutyp
    real,                        intent(in) :: dt     !time step (sec)
    integer,                     intent(in) :: nsoil  !no. of soil layers
    integer,                     intent(in) :: nsnow  !maximum no. of snow layers
    real, dimension(1:NSOIL),    intent(in) :: zsoil  !depth of soil layer-bottom [m]
    real, dimension(-NSNOW+1:NSOIL), intent(in) :: dzsnso !snow/soil layer depth [m]
    integer, intent(in) :: slptyp
    real, intent(in)                        :: qinsrf !water input on soil surface [mm/s]
    real, intent(in)                        :: QSEVA  !evap from soil surface [mm/s]
    real, dimension(1:NSOIL),    intent(in) :: ETRANI !evapotranspiration from soil layers [mm/s]
    real, dimension(1:NSOIL), intent(in)    :: soilice !soil ice content [m3/m3]

    ! input & output
    real, dimension(1:NSOIL), intent(inout) :: soilwat !soil liquid water content [m3/m3]
    real, dimension(1:NSOIL), intent(inout) :: SMC    !total soil water content [m3/m3]
    real, intent(inout)                     :: ZWT    !water table depth [m]

    ! output
    real, intent(out)                       :: QDRAIN !soil-bottom free drainage [mm/s]
    real, intent(out)                       :: runsrf !surface runoff [mm/s]
    real, intent(out)                       :: runsub !subsurface runoff [mm/s]
    real, intent(out)                       :: FCRMAX !maximum of FCR (-)
    real, dimension(1:NSOIL), intent(out)   :: WCND   !hydraulic conductivity (m/s)

    ! local
    integer                                 :: iz     !do-loop index
    integer                                 :: iter   !iteration index
    real                                    :: dtfine !fine time step (s)
    real, dimension(1:NSOIL)                :: RHSTT  !right-hand side term of the matrix
    real, dimension(1:NSOIL)                :: AI     !left-hand side term
    real, dimension(1:NSOIL)                :: BI     !left-hand side term
    real, dimension(1:NSOIL)                :: CI     !left-hand side term

    real                                    :: FFF    !runoff decay factor (m-1)
    real                                    :: RSBMX  !baseflow coefficient [mm/s]
    real                                    :: qinfil  !infiltration rate at surface (m/s)
    real                                    :: FICE   !ice fraction in frozen soil
    real                                    :: WPLUS  !saturation excess of the total soil [m]
    real                                    :: RSAT   !accumulation of WPLUS (saturation excess) [m]
    real                                    :: SICEMAX!maximum soil ice content (m3/m3)
    real                                    :: SH2OMIN!minimum soil liquid water content (m3/m3)
    real                                    :: WTSUB  !sum of WCND(K)*DZSNSO(K)
    real                                    :: MH2O   !water mass removal (mm)
    real                                    :: FSAT   !fractional saturated area (-)
    real, dimension(1:NSOIL)                :: MLIQ   !
    real                                    :: XS     !
    real                                    :: WATMIN !
    real                                    :: QDRAIN_SAVE !
    real                                    :: epore  !effective porosity [m3/m3]
    real, dimension(1:NSOIL)                :: FCR    !impermeable fraction due to frozen soil
    integer                                 :: niter  !iteration times soil moisture (-)
    real                                    :: SMCTOT !2-m averaged soil moisture (m3/m3)
    real                                    :: DZTOT  !2-m soil depth (m)
    real, parameter :: A = 4.0

    runsrf = 0.0
    qinfil  = 0.0
    RSAT   = 0.0

    ! for the case when snowmelt water is too large
    do iz = 1, nsoil
       epore   = max(1.0E-4 , (LK_SMCMAX(sltyp) - soilice(iz)))
       RSAT    = RSAT + max(0.0, soilwat(iz) - epore) * DZSNSO(iz)
       soilwat(iz) = min(epore, soilwat(iz))
    end do

    !impermeable fraction due to frozen soil
    do iz = 1, nsoil
       FICE = min(1.0, soilice(iz) / LK_SMCMAX(sltyp))
       FCR(iz) = max(0.0, exp(-A * (1.0 - FICE)) - exp(-A)) /  &
            (1.0 - exp(-A))
    end do

    ! maximum soil ice content and minimum liquid water of all layers
    SICEMAX = 0.0
    FCRMAX  = 0.0
    SH2OMIN = LK_SMCMAX(sltyp)
    do iz = 1, NSOIL
       if (soilice(iz) > SICEMAX) SICEMAX = soilice(iz)
       if (FCR(iz)  > FCRMAX)  FCRMAX  = FCR(iz)
       if (soilwat(iz) < SH2OMIN) SH2OMIN = soilwat(iz)
    end do

    !subsurface runoff for runoff scheme option 2
    if (opt_run == 2) then
       FFF   = 2.0
       RSBMX = 4.0
       call zwteq(sltyp, nsoil, nsnow, ZSOIL, DZSNSO, soilwat, ZWT)
       runsub = (1.0 - FCRMAX) * RSBMX * exp(-KK_TIMEAN) * exp(-FFF * ZWT)   ! mm/s
    end if

    !surface runoff and infiltration rate using different schemes

    !jref impermable surface at urban
    if (lutyp == ISURBAN) FCR(1) = 0.95

    if (opt_run == 1) then
       FFF = 6.0
       FSAT = KK_FSATMAX * exp(-0.5 * FFF * (ZWT - 2.0))
       if (qinsrf > 0.0) then
          runsrf = qinsrf * ((1.0 - FCR(1)) * FSAT + FCR(1))
          qinfil  = qinsrf - runsrf                          ! m/s
       end if
    end if

    if (opt_run == 2) then
       FFF = 2.0
       FSAT = KK_FSATMAX * exp(-0.5 * FFF * ZWT)
       if (qinsrf > 0.0) then
          runsrf = qinsrf * ((1.0 - FCR(1)) * FSAT + FCR(1))
          qinfil  = qinsrf - runsrf                          ! m/s
       end if
    end if

    if (opt_run == 3) then
       call infil(sltyp, dt, nsoil, zsoil, soilwat   ,soilice   , & !in
            SICEMAX,qinsrf ,                         & !in
            qinfil  ,runsrf )                           !out
    end if

    if (opt_run == 4) then
       SMCTOT = 0.0
       DZTOT  = 0.0
       do iz = 1, NSOIL
          DZTOT   = DZTOT  + DZSNSO(iz)
          SMCTOT  = SMCTOT + SMC(iz) * DZSNSO(iz)
          if (DZTOT >= 2.0) exit
       end do
       SMCTOT = SMCTOT / DZTOT
       FSAT   = max(0.01, SMCTOT / LK_SMCMAX(sltyp)) ** 4.0        !BATS

       if (qinsrf > 0.0) then
          runsrf = qinsrf * ((1.0 - FCR(1)) * FSAT + FCR(1))
          qinfil  = qinsrf - runsrf                       ! m/s
       end if
    end if

    ! determine iteration times and finer time step
    niter = 1
    if (opt_inf == 1) then    !opt_inf =2 may cause water imbalance
       niter = 3
       if (qinfil * dt > dzsnso(1) * LK_SMCMAX(sltyp)) then
          niter = niter * 2
       end if
    end if
    dtfine  = dt / niter

    ! solve soil moisture
    QDRAIN_SAVE = 0.0
    do iter = 1, niter
       call srt(sltyp, dtfine, nsoil, zsoil, slptyp, qinfil  ,ETRANI , & !in
            QSEVA  ,soilwat   ,SMC    ,ZWT    ,FCR    , & !in
            SICEMAX,FCRMAX, & !in
            RHSTT  ,AI     ,BI     ,CI     ,QDRAIN , & !out
            WCND   )                                   !out

       call sstep(sltyp, dtfine, nsoil, nsnow, zsoil, dzsnso, & !in
            soilice   ,ZWT            ,                 & !in
            soilwat   ,SMC    ,AI     ,BI     ,CI     , & !inout
            RHSTT, QDRAIN,  & !inout
            WPLUS)                                     !out
       RSAT =  RSAT + WPLUS
       QDRAIN_SAVE = QDRAIN_SAVE + QDRAIN
    end do

    QDRAIN = QDRAIN_SAVE / niter

    runsrf = runsrf * 1000.0 + RSAT * 1000.0 / dt  ! m/s -> mm/s
    QDRAIN = QDRAIN * 1000.0

    ! removal of soil water due to groundwater flow (option 2)
    if (opt_run == 2) then
       WTSUB = 0.0
       do iz = 1, nsoil
          WTSUB = WTSUB + WCND(iz) * DZSNSO(iz)
       end do

       do iz = 1, nsoil
          MH2O    = runsub * DT * (WCND(iz) * DZSNSO(iz)) / WTSUB       ! mm
          soilwat(iz) = soilwat(iz) - MH2O / (DZSNSO(iz) * 1000.0)
       end do
    end if

    ! Limit MLIQ to be greater than or equal to watmin.
    ! Get water needed to bring MLIQ equal WATMIN from lower layer.
    if (opt_run /= 1) then
       do iz = 1, nsoil
          MLIQ(iz) = soilwat(iz) * DZSNSO(iz) * 1000.0
       end do

       WATMIN = 0.01           ! mm
       do iz = 1, nsoil-1
          if (MLIQ(iz) < 0.0) then
             XS = WATMIN - MLIQ(iz)
          else
             XS = 0.0
          end if
          MLIQ(iz) = MLIQ(iz) + XS
          MLIQ(iz+1) = MLIQ(iz+1) - XS
       end do

       iz = NSOIL
       if (MLIQ(iz) < WATMIN) then
          XS = WATMIN - MLIQ(iz)
       else
          XS = 0.0
       end if
       MLIQ(iz) = MLIQ(iz) + XS
       runsub   = runsub - XS / DT

       do iz = 1, NSOIL
          soilwat(iz) = MLIQ(iz) / (DZSNSO(iz) * 1000.0)
       end do
    end if

  end subroutine soilh2o


  subroutine zwteq(sltyp, nsoil, nsnow, ZSOIL, DZSNSO, soilwat, ZWT)
    ! calculate equilibrium water table depth (Niu et al., 2005)
    use noahmp_soil_param, only: LK_BEXP
    use noahmp_soil_param, only: LK_PSISAT
    use noahmp_soil_param, only: LK_SMCMAX
    implicit none
    ! input
    integer, intent(in) :: sltyp
    integer,                         intent(in) :: NSOIL  !no. of soil layers
    integer,                         intent(in) :: NSNOW  !maximum no. of snow layers
    real, dimension(1:NSOIL),        intent(in) :: ZSOIL  !depth of soil layer-bottom [m]
    real, dimension(-NSNOW+1:NSOIL), intent(in) :: DZSNSO !snow/soil layer depth [m]
    real, dimension(1:NSOIL),        intent(in) :: soilwat   !soil liquid water content [m3/m3]

    ! output

    real,                           intent(out) :: ZWT    !water table depth [m]

    ! locals

    integer :: K                      !do-loop index
    integer, parameter :: NFINE = 100 !no. of fine soil layers of 6m soil
    real    :: WD1                    !water deficit from coarse (4-L) soil moisture profile
    real    :: WD2                    !water deficit from fine (100-L) soil moisture profile
    real    :: DZFINE                 !layer thickness of the 100-L soil layers to 6.0 m
    real    :: TEMP                   !temporary variable
    real, dimension(1:NFINE) :: ZFINE !layer-bottom depth of the 100-L soil layers to 6.0 m

    WD1 = 0.0
    do K = 1,NSOIL
       WD1 = WD1 + (LK_SMCMAX(sltyp) - soilwat(K)) * DZSNSO(K) ! [m]
    end do

    DZFINE = 3.0 * (-ZSOIL(NSOIL)) / NFINE
    do K =1, NFINE
       ZFINE(K) = FLOAT(K) * DZFINE
    end do

    ZWT = -3.0 * ZSOIL(NSOIL) - 0.001   ! initial value [m]

    WD2 = 0.0
    do K = 1, NFINE
       TEMP  = 1.0 + (ZWT - ZFINE(K)) / LK_PSISAT(sltyp)
       WD2   = WD2 + LK_SMCMAX(sltyp) * (1.0 - TEMP ** (-1.0 / LK_BEXP(sltyp))) * DZFINE
       if (abs(WD2 - WD1) <= 0.01) then
          ZWT = ZFINE(K)
          exit
       end if
    end do
  end subroutine zwteq


  subroutine infil(sltyp, dt, nsoil, zsoil, soilwat, soilice, & !in
       & SICEMAX, qinsrf ,                         & !in
       & qinfil, runsrf )                           !out
    ! compute inflitration rate at soil surface and surface runoff
    use noahmp_soil_param, only: LK_KDT
    use noahmp_soil_param, only: LK_FRZX
    use noahmp_soil_param, only: LK_SMCMAX
    use noahmp_soil_param, only: LK_SMCWLT
    implicit none
    ! inputs
    integer, intent(in) :: sltyp
    real,                     intent(in) :: dt     !time step (sec)
    integer,                  intent(in) :: nsoil  !no. of soil layers
    real, dimension(1:NSOIL), intent(in) :: zsoil  !depth of soil layer-bottom [m]
    real, dimension(1:NSOIL), intent(in) :: soilwat   !soil liquid water content [m3/m3]
    real, dimension(1:NSOIL), intent(in) :: soilice   !soil ice content [m3/m3]
    real,                     intent(in) :: qinsrf !water input on soil surface [mm/s]
    real,                     intent(in) :: SICEMAX!maximum soil ice content (m3/m3)

    ! outputs
    real,                    intent(out) :: runsrf !surface runoff [mm/s]
    real,                    intent(out) :: qinfil !infiltration rate at surface

    ! locals
    integer :: IALP1, J, JJ,  K
    real                     :: VAL
    real                     :: DDT
    real                     :: PX
    real                     :: DT1, DD, DICE
    real                     :: FCR
    real                     :: SUM
    real                     :: ACRT
    real                     :: WDF
    real                     :: WCND
    real                     :: SMCAV
    real                     :: INFMAX
    real, dimension(1:NSOIL) :: DMAX
    integer, parameter       :: CVFRZ = 3

    if (qinsrf >  0.0) then
       DT1 = DT / 86400.0
       SMCAV = LK_SMCMAX(sltyp) - LK_SMCWLT(sltyp)

       ! maximum infiltration rate
       DMAX(1)= -ZSOIL(1) * SMCAV
       DICE   = -ZSOIL(1) * soilice(1)
       DMAX(1)= DMAX(1)* (1.0 - (soilwat(1) + soilice(1) - LK_SMCWLT(sltyp)) / SMCAV)

       DD = DMAX(1)

       do K = 2, NSOIL
          DICE    = DICE + (ZSOIL(K-1) - ZSOIL(K) ) * soilice(K)
          DMAX(K) = (ZSOIL(K-1) - ZSOIL(K)) * SMCAV
          DMAX(K) = DMAX(K) * (1.0-(soilwat(K) + soilice(K) - LK_SMCWLT(sltyp)) / SMCAV)
          DD      = DD + DMAX(K)
       end do

       VAL = (1.0 - exp (-LK_KDT(sltyp) * dt1))
       DDT = DD * VAL
       PX  = max(0.0, qinsrf * DT)
       INFMAX = (PX * (DDT / (PX + DDT)))/ DT

       ! impermeable fraction due to frozen soil

       FCR = 1.0
       if (DICE >  1.0E-2) then
          ACRT = CVFRZ * LK_FRZX(sltyp) / DICE
          SUM = 1.0
          IALP1 = CVFRZ - 1
          do J = 1, IALP1
             K = 1
             do JJ = J+1, IALP1
                K = K * JJ
             end do
             SUM = SUM + (ACRT ** (CVFRZ - J)) / FLOAT(K)
          end do
          FCR = 1.0 - exp(-ACRT) * SUM
       end if

       ! correction of infiltration limitation

       INFMAX = INFMAX * FCR

       ! jref for urban areas
       !       IF (lutyp == ISURBAN ) INFMAX == INFMAX * 0.05

       call wdfcnd2(sltyp, WDF, WCND, soilwat(1), SICEMAX)
       INFMAX = max(INFMAX, WCND)
       INFMAX = min(INFMAX, PX)

       runsrf= max(0.0, qinsrf - INFMAX)
       qinfil = qinsrf - runsrf
    end if
  end subroutine infil


  subroutine srt(sltyp, dt, nsoil, zsoil, slptyp, &
       &         qinfil, ETRANI, QSEVA, soilwat, SMC, &
       &         ZWT, FCR, SICEMAX, FCRMAX, & !in
       &         RHSTT, AI, BI, CI, QDRAIN, & !out
       &         WCND)
    ! calculate the right hand side of the time tendency term of the soil
    ! water diffusion equation.  also to compute ( prepare ) the matrix
    ! coefficients for the tri-diagonal matrix of the implicit time scheme.
    use noahmp_gen_param, only: LK_SLOPE
    implicit none
    !input
    integer, intent(in) :: sltyp
    real,                     intent(in)  :: dt
    integer,                  intent(in)  :: nsoil
    real, dimension(1:NSOIL), intent(in)  :: zsoil
    integer, intent(in) :: slptyp
    real,                     intent(in)  :: qinfil
    real,                     intent(in)  :: QSEVA
    real, dimension(1:NSOIL), intent(in)  :: ETRANI
    real, dimension(1:NSOIL), intent(in)  :: soilwat
    real, dimension(1:NSOIL), intent(in)  :: SMC
    real,                     intent(in)  :: ZWT    ! water table depth [m]
    real, dimension(1:NSOIL), intent(in)  :: FCR
    real, intent(in)                      :: FCRMAX !maximum of FCR (-)
    real,                     intent(in)  :: SICEMAX!maximum soil ice content (m3/m3)

    ! output

    real, dimension(1:NSOIL), intent(out) :: RHSTT
    real, dimension(1:NSOIL), intent(out) :: AI
    real, dimension(1:NSOIL), intent(out) :: BI
    real, dimension(1:NSOIL), intent(out) :: CI
    real, dimension(1:NSOIL), intent(out) :: WCND    !hydraulic conductivity (m/s)
    real,                     intent(out) :: QDRAIN  !bottom drainage (m/s)

    ! local
    integer                               :: K
    real, dimension(1:NSOIL)              :: DDZ
    real, dimension(1:NSOIL)              :: DENOM
    real, dimension(1:NSOIL)              :: DSMDZ
    real, dimension(1:NSOIL)              :: WFLUX
    real, dimension(1:NSOIL)              :: WDF
    real, dimension(1:NSOIL)              :: SMX
    real                                  :: TEMP1

    ! Niu and Yang (2006), J. of Hydrometeorology

    if (opt_inf == 1) then
       do K = 1, nsoil
          call wdfcnd1(sltyp, WDF(K), WCND(K), SMC(K), FCR(K))
          SMX(K) = SMC(K)
       end do
    end if

    if (opt_inf == 2) then
       do K = 1, NSOIL
          call wdfcnd2(sltyp, WDF(K), WCND(K), soilwat(K), SICEMAX)
          SMX(K) = soilwat(K)
       end do
    end if

    do K = 1, nsoil
       if (K == 1) then
          DENOM(K) = -ZSOIL(K)
          TEMP1    = -ZSOIL(K+1)
          DDZ(K)   = 2.0 / TEMP1
          DSMDZ(K) = 2.0 * (SMX(K) - SMX(K+1)) / TEMP1
          WFLUX(K) = WDF(K) * DSMDZ(K) + WCND(K) - qinfil + ETRANI(K) + QSEVA
       else if (K < NSOIL) then
          DENOM(k) = (ZSOIL(K-1) - ZSOIL(K))
          TEMP1    = (ZSOIL(K-1) - ZSOIL(K+1))
          DDZ(K)   = 2.0 / TEMP1
          DSMDZ(K) = 2.0 * (SMX(K) - SMX(K+1)) / TEMP1
          WFLUX(K) = WDF(K) * DSMDZ(K) + WCND(K)         &
               - WDF(K-1) * DSMDZ(K-1) - WCND(K-1) + ETRANI(K)
       else
          DENOM(K) = (ZSOIL(K-1) - ZSOIL(K))
          if (opt_run == 1 .or. opt_run == 2) then
             QDRAIN   = 0.0
          end if
          if (opt_run == 3) then
             QDRAIN   = LK_SLOPE(slptyp) * WCND(K)
          end if
          if (opt_run == 4) then
             QDRAIN   = (1.0 - FCRMAX) * WCND(K)
          end if
          WFLUX(K) = -(WDF(K-1) * DSMDZ(K-1)) - WCND(K-1) + ETRANI(K) + QDRAIN
       end if
    end do

    do K = 1, nsoil
       if (K == 1) then
          AI(K)    =  0.0
          BI(K)    =  WDF(K) * DDZ(K) / DENOM(K)
          CI(K)    = -BI(K)
       else if (K < nsoil) then
          AI(K)    = -WDF(K-1) * DDZ(K-1) / DENOM(K)
          CI(K)    = -WDF(K) * DDZ(K) / DENOM(K)
          BI(K)    = -(AI(K) + CI(K))
       else
          AI(K)    = -WDF(K-1) * DDZ(K-1) / DENOM(K)
          CI(K)    = 0.0
          BI(K)    = -(AI(K) + CI(K))
       end if
       RHSTT(K) = WFLUX(K) / (-DENOM(K))
    end do
  end subroutine srt


  subroutine sstep(sltyp, dt, nsoil, nsnow, zsoil, dzsnso, & !in
       soilice   ,ZWT            ,                 & !in
       soilwat   ,SMC    ,AI     ,BI     ,CI     , & !inout
       RHSTT, QDRAIN, & !inout
       WPLUS  )                                   !out
    ! calculate/update soil moisture content values
    use noahmp_soil_param, only: LK_SMCMAX
    implicit none
    integer, intent(in) :: sltyp
    real, intent(in)                            :: dt
    integer,                         intent(in) :: nsoil  !
    integer,                         intent(in) :: nsnow  !
    real, dimension(       1:NSOIL), intent(in) :: zsoil
    real, dimension(-NSNOW+1:NSOIL), intent(in) :: dzsnso ! snow/soil layer thickness [m]
    real, intent(in)                            :: ZWT
    real, dimension(       1:NSOIL), intent(in) :: soilice

    !input and output
    real, dimension(1:NSOIL), intent(inout) :: soilwat
    real, dimension(1:NSOIL), intent(inout) :: SMC
    real, dimension(1:NSOIL), intent(inout) :: AI
    real, dimension(1:NSOIL), intent(inout) :: BI
    real, dimension(1:NSOIL), intent(inout) :: CI
    real, dimension(1:NSOIL), intent(inout) :: RHSTT
    real                    , intent(inout) :: QDRAIN

    !output
    real, intent(out)                       :: WPLUS     !saturation excess water (m)

    !local
    integer                                 :: K
    real, dimension(1:NSOIL)                :: RHSTTIN
    real, dimension(1:NSOIL)                :: CIIN
    real                                    :: STOT
    real                                    :: epore
    real                                    :: WMINUS

    WPLUS = 0.0

    do K = 1,NSOIL
       RHSTT (K) =   RHSTT(K) * DT
       AI (K)    =      AI(K) * DT
       BI (K)    = 1.0 + BI(K) * DT
       CI (K)    =      CI(K) * DT
    end do

    ! copy values for input variables before calling rosr12

    do K = 1,NSOIL
       RHSTTIN(k) = RHSTT(K)
       CIIN(k)    = CI(K)
    end do

    ! call rosr12 to solve the tri-diagonal matrix

    call rosr12(CI, AI, BI, CIIN, RHSTTIN, RHSTT, 1, NSOIL, 0)

    do K = 1, NSOIL
       soilwat(K) = soilwat(K) + CI(K)
    end do

    !  excessive water above saturation in a layer is moved to
    !  its unsaturated layer like in a bucket

    do K = NSOIL,2,-1
       epore = max(1.0E-4, LK_SMCMAX(sltyp) - soilice(K))
       WPLUS = max(soilwat(K) - epore, 0.0) * DZSNSO(K)
       soilwat(K) = min(epore, soilwat(K))
       soilwat(K-1) = soilwat(K-1) + WPLUS/DZSNSO(K-1)
    end do

    epore = max(1.0E-4, LK_SMCMAX(sltyp) - soilice(1))
    WPLUS = max(soilwat(1) - epore, 0.0) * DZSNSO(1)
    soilwat(1) = min(epore, soilwat(1))
    SMC = soilwat + soilice
  end subroutine sstep


  subroutine wdfcnd1(sltyp, WDF, WCND, SMC, FCR)
    ! calculate soil water diffusivity and soil hydraulic conductivity.
    use noahmp_soil_param, only: LK_BEXP
    use noahmp_soil_param, only: LK_SMCMAX
    use noahmp_soil_param, only: LK_DKSAT
    use noahmp_soil_param, only: LK_DWSAT
    implicit none
    ! input
    integer, intent(in) :: sltyp
    real, intent(in)  :: SMC
    real, intent(in)  :: FCR

    ! output
    real, intent(out) :: WCND
    real, intent(out) :: WDF

    ! local
    real :: EXPON
    real :: FACTR
    real :: VKWGT

    ! soil water diffusivity
    FACTR = max(0.01, SMC / LK_SMCMAX(sltyp))
    EXPON = LK_BEXP(sltyp) + 2.0
    WDF   = LK_DWSAT(sltyp) * FACTR ** EXPON
    WDF   = WDF * (1.0 - FCR)

    ! hydraulic conductivity
    EXPON = 2.0 * LK_BEXP(sltyp) + 3.0
    WCND  = LK_DKSAT(sltyp) * FACTR ** EXPON
    WCND  = WCND * (1.0 - FCR)
  end subroutine wdfcnd1


  subroutine wdfcnd2(sltyp, WDF, WCND, SMC, soilice)
    ! calculate soil water diffusivity and soil hydraulic conductivity.
    use noahmp_soil_param, only: LK_BEXP
    use noahmp_soil_param, only: LK_SMCMAX
    use noahmp_soil_param, only: LK_DKSAT
    use noahmp_soil_param, only: LK_DWSAT
    implicit none
    ! input
    integer, intent(in) :: sltyp
    real,intent(in)  :: SMC
    real,intent(in)  :: soilice

    ! output
    real,intent(out) :: WCND
    real,intent(out) :: WDF

    ! local
    real :: EXPON
    real :: FACTR
    real :: VKWGT
    ! soil water diffusivity

    FACTR = max(0.01, SMC / LK_SMCMAX(sltyp))
    EXPON = LK_BEXP(sltyp) + 2.0
    WDF   = LK_DWSAT(sltyp) * FACTR ** EXPON

    if (soilice > 0.0) then
       VKWGT = 1.0 / (1.0 + (500.0 * soilice) ** 3.0)
       WDF   = VKWGT * WDF &
            & + (1.0 - VKWGT) * LK_DWSAT(sltyp) * (0.2 / LK_SMCMAX(sltyp)) ** EXPON
    end if

    ! hydraulic conductivity
    EXPON = 2.0 * LK_BEXP(sltyp) + 3.0
    WCND  = LK_DKSAT(sltyp) * FACTR ** EXPON
  end subroutine wdfcnd2


  subroutine groundwater(sltyp, dt, nsnow, nsoil, zsoil, soilice, & !in
       & STC, WCND, FCRMAX, ILOC, JLOC, & !in
       & soilwat, ZWT, WA, WT, & !inout
       & QIN, QDIS)                           !out
    use noahmp_gen_param, only: KK_TIMEAN
    use noahmp_soil_param, only: LK_BEXP
    use noahmp_soil_param, only: LK_PSISAT
    use noahmp_soil_param, only: LK_SMCMAX
    implicit none
    ! input
    integer, intent(in) :: sltyp
    integer,                         intent(in) :: ILOC  !grid index
    integer,                         intent(in) :: JLOC  !grid index
    integer,                         intent(in) :: NSNOW !maximum no. of snow layers
    integer,                         intent(in) :: NSOIL !no. of soil layers
    real,                            intent(in) :: DT    !timestep [sec]
    real,                            intent(in) :: FCRMAX!maximum FCR (-)
    real, dimension(       1:NSOIL), intent(in) :: soilice  !soil ice content [m3/m3]
    real, dimension(       1:NSOIL), intent(in) :: ZSOIL !depth of soil layer-bottom [m]
    real, dimension(       1:NSOIL), intent(in) :: WCND  !hydraulic conductivity (m/s)
    real, dimension(-NSNOW+1:NSOIL), intent(in) :: STC   !snow/soil temperature (k)

    ! input and output
    real, dimension(    1:NSOIL), intent(inout) :: soilwat  !liquid soil water [m3/m3]
    real,                         intent(inout) :: ZWT   !the depth to water table [m]
    real,                         intent(inout) :: WA    !water storage in aquifer [mm]
    real,                         intent(inout) :: WT    !water storage in aquifer
    !+ saturated soil [mm]
    ! output
    real,                           intent(out) :: QIN   !groundwater recharge [mm/s]
    real,                           intent(out) :: QDIS  !groundwater discharge [mm/s]

    ! local
    real                                        :: FFF   !runoff decay factor (m-1)
    real                                        :: RSBMX !baseflow coefficient [mm/s]
    integer                                     :: IZ    !do-loop index
    integer                                     :: IWT   !layer index above water table layer
    real,  dimension(    1:NSOIL)               :: DZMM  !layer thickness [mm]
    real,  dimension(    1:NSOIL)               :: ZNODE !node depth [m]
    real,  dimension(    1:NSOIL)               :: MLIQ  !liquid water mass [kg/m2 or mm]
    real,  dimension(    1:NSOIL)               :: epore !effective porosity [-]
    real,  dimension(    1:NSOIL)               :: HK    !hydraulic conductivity [mm/s]
    real,  dimension(    1:NSOIL)               :: SMC   !total soil water  content [m3/m3]
    real(KIND=8)                                :: S_NODE!degree of saturation of IWT layer
    real                                        :: DZSUM !cumulative depth above water table [m]
    real                                        :: SMPFZ !matric potential (frozen effects) [mm]
    real                                        :: KA    !aquifer hydraulic conductivity [mm/s]
    real                                        :: WH_ZWT!water head at water table [mm]
    real                                        :: WH    !water head at layer above ZWT [mm]
    real                                        :: WS    !water used to fill air pore [mm]
    real                                        :: WTSUB !sum of HK*DZMM
    real                                        :: WATMIN!minimum soil vol soil moisture [m3/m3]
    real                                        :: XS    !excessive water above saturation [mm]
    real, parameter                             :: ROUS = 0.2    !specific yield [-]
    real, parameter                             :: CMIC = 0.20   !microprore content (0.0-1.0)
    !0.0-close to free drainage

    QDIS      = 0.0
    QIN       = 0.0

    ! Derive layer-bottom depth in [mm]
    !KWM:  Derive layer thickness in mm

    DZMM(1) = -ZSOIL(1) * 1.0E3
    do IZ = 2, NSOIL
       DZMM(IZ)  = 1.0E3 * (ZSOIL(IZ-1) - ZSOIL(IZ))
    end do

    ! Derive node (middle) depth in [m]
    !KWM:  Positive number, depth below ground surface in m
    ZNODE(1) = -ZSOIL(1) / 2.0
    do IZ = 2, NSOIL
       ZNODE(IZ)  = -ZSOIL(IZ-1) + 0.5 * (ZSOIL(IZ-1) - ZSOIL(IZ))
    end do

    ! Convert volumetric soil moisture "soilwat" to mass

    do IZ = 1, NSOIL
       SMC(IZ)      = soilwat(IZ) + soilice(IZ)
       MLIQ(IZ)     = soilwat(IZ) * DZMM(IZ)
       epore(IZ)    = max(0.01, LK_SMCMAX(sltyp) - soilice(IZ))
       HK(IZ)       = 1.0E3 * WCND(IZ)
    end do

    ! The layer index of the first unsaturated layer,
    ! i.e., the layer right above the water table

    IWT = NSOIL
    do IZ = 2, NSOIL
       if (ZWT <= -ZSOIL(IZ) ) then
          IWT = IZ - 1
          exit
       end if
    end do

    ! Groundwater discharge [mm/s]
    FFF   = 6.0
    RSBMX = 5.0

    QDIS = (1.0 - FCRMAX) * RSBMX * exp(-KK_TIMEAN) * exp(-FFF * (ZWT-2.0))

    ! Matric potential at the layer above the water table
    S_NODE = min(1.0, SMC(IWT) / LK_SMCMAX(sltyp))
    S_NODE = max(S_NODE, real(0.01,KIND=8))
    SMPFZ  = -LK_PSISAT(sltyp) * 1000.0 * S_NODE ** (-LK_BEXP(sltyp))   ! m --> mm
    SMPFZ  = max(-120000.0, CMIC * SMPFZ)

    ! Recharge rate qin to groundwater
    KA  = HK(IWT)

    WH_ZWT  = -ZWT * 1.0E3                          !(mm)
    WH      = SMPFZ  - ZNODE(IWT) * 1.0E3              !(mm)
    QIN     = -KA * (WH_ZWT - WH) / ((ZWT - ZNODE(IWT)) * 1.0E3)
    QIN     = max(-10.0 / DT, min(10.0 / DT, QIN))

    ! water storage in the aquifer + saturated soil

    WT  = WT + (QIN - QDIS) * DT     !(mm)

    if (IWT==NSOIL) then
       WA          = WA + (QIN - QDIS) * DT     !(mm)
       WT          = WA
       ZWT         = (-ZSOIL(NSOIL) + 25.0) - WA / 1000.0 / ROUS      !(m)
       MLIQ(NSOIL) = MLIQ(NSOIL) - QIN * DT        ! [mm]

       MLIQ(NSOIL) = MLIQ(NSOIL) + max(0.0, WA - 5000.0)
       WA          = min(WA, 5000.0)
    else

       if (IWT == NSOIL - 1) then
          ZWT = -ZSOIL(NSOIL) &
               & - (WT - ROUS * 1000.0 * 25.0) / (epore(NSOIL)) / 1000.0
       else
          WS = 0.0   ! water used to fill soil air pores
          do IZ = IWT+2, NSOIL
             WS = WS + epore(IZ) * DZMM(IZ)
          end do
          ZWT = -ZSOIL(IWT+1)                  &
               - (WT-ROUS * 1000.0 * 25.0 - WS) / (epore(IWT+1)) / 1000.0
       end if

       WTSUB = 0.0
       do IZ = 1, NSOIL
          WTSUB = WTSUB + HK(IZ) * DZMM(IZ)
       end do

       do IZ = 1, NSOIL           ! Removing subsurface runoff
          MLIQ(IZ) = MLIQ(IZ) - QDIS * DT * HK(IZ) * DZMM(IZ) / WTSUB
       end do
    end if

    ZWT = max(1.5, ZWT)

    !
    ! Limit MLIQ to be greater than or equal to watmin.
    ! Get water needed to bring MLIQ equal WATMIN from lower layer.
    !
    WATMIN = 0.01
    do IZ = 1, NSOIL-1
       if (MLIQ(IZ) < 0.0) then
          XS = WATMIN - MLIQ(IZ)
       else
          XS = 0.0
       end if
       MLIQ(IZ) = MLIQ(IZ) + XS
       MLIQ(IZ+1) = MLIQ(IZ+1) - XS
    end do

    IZ = NSOIL
    if (MLIQ(IZ) < WATMIN) then
       XS = WATMIN - MLIQ(IZ)
    else
       XS = 0.0
    end if
    MLIQ(IZ) = MLIQ(IZ) + XS
    WA = WA - XS
    WT = WT - XS

    do IZ = 1, NSOIL
       soilwat(IZ) = MLIQ(IZ) / DZMM(IZ)
    end do
  end subroutine groundwater


  subroutine carbon(NSNOW, NSOIL, lutyp, DT, ZSOIL, & !in
       DZSNSO, STC,    SMC,    TV,     TG,     PSN,    & !in
       FOLN,   SMCMAX, BTRAN,  APAR,   fveg,   IGS,    & !in
       TROOT,  IST,    LAT,    ILOC,   JLOC, & !in
       LFMASS, RTMASS, STMASS, WOOD,   STBLCP, FASTCP , & !inout
       GPP,    NPP,    NEE,    AUTORS, HETERS, TOTSC  , & !out
       TOTLB,  XLAI,   XSAI)                   !out
    use noahmp_veg_param, only: ISURBAN
    use noahmp_veg_param, only: ISWATER
    use noahmp_veg_param, only: ISBARREN
    use noahmp_veg_param, only: ISICE
    use noahmp_veg_param, only: LK_NROOT
    use noahmp_veg_param, only: LK_SLA
    implicit none
    ! inputs (carbon)
    integer                        , intent(in) :: ILOC   !grid index
    integer                        , intent(in) :: JLOC   !grid index
    integer                        , intent(in) :: lutyp !vegetation type
    integer                        , intent(in) :: NSNOW  !number of snow layers
    integer                        , intent(in) :: NSOIL  !number of soil layers
    real                           , intent(in) :: LAT    !latitude (radians)
    real                           , intent(in) :: DT     !time step (s)
    real, dimension(       1:NSOIL), intent(in) :: ZSOIL  !depth of layer-bottom from soil surface
    real, dimension(-NSNOW+1:NSOIL), intent(in) :: DZSNSO !snow/soil layer thickness [m]
    real, dimension(-NSNOW+1:NSOIL), intent(in) :: STC    !snow/soil temperature [k]
    real, dimension(       1:NSOIL), intent(in) :: SMC    !soil moisture (ice + liq.) [m3/m3]
    real                           , intent(in) :: TV     !vegetation temperature (k)
    real                           , intent(in) :: TG     !ground temperature (k)
    real                           , intent(in) :: FOLN   !foliage nitrogen (%)
    real                           , intent(in) :: SMCMAX !soil porosity (m3/m3)
    real                           , intent(in) :: BTRAN  !soil water transpiration factor (0 to 1)
    real                           , intent(in) :: PSN    !total leaf photosyn (umolco2/m2/s) [+]
    real                           , intent(in) :: APAR   !PAR by canopy (w/m2)
    real                           , intent(in) :: IGS    !growing season index (0=off, 1=on)
    real                           , intent(in) :: fveg   !vegetation greenness fraction
    real                           , intent(in) :: TROOT  !root-zone averaged temperature (k)
    integer                        , intent(in) :: IST    !surface type 1->soil; 2->lake

    ! input & output (carbon)

    real                        , intent(inout) :: LFMASS !leaf mass [g/m2]
    real                        , intent(inout) :: RTMASS !mass of fine roots [g/m2]
    real                        , intent(inout) :: STMASS !stem mass [g/m2]
    real                        , intent(inout) :: WOOD   !mass of wood (incl. woody roots) [g/m2]
    real                        , intent(inout) :: STBLCP !stable carbon in deep soil [g/m2]
    real                        , intent(inout) :: FASTCP !short-lived carbon in shallow soil [g/m2]

    ! outputs: (carbon)

    real                          , intent(out) :: GPP    !net instantaneous assimilation [g/m2/s C]
    real                          , intent(out) :: NPP    !net primary productivity [g/m2/s C]
    real                          , intent(out) :: NEE    !net ecosystem exchange [g/m2/s CO2]
    real                          , intent(out) :: AUTORS !net ecosystem respiration [g/m2/s C]
    real                          , intent(out) :: HETERS !organic respiration [g/m2/s C]
    real                          , intent(out) :: TOTSC  !total soil carbon [g/m2 C]
    real                          , intent(out) :: TOTLB  !total living carbon ([g/m2 C]
    real                          , intent(out) :: XLAI   !leaf area index [-]
    real                          , intent(out) :: XSAI   !stem area index [-]
    !  REAL                          , INTENT(out) :: VOCFLX(5) ! voc fluxes [ug C m-2 h-1]

    ! local variables

    integer :: J         !do-loop index
    real    :: WROOT     !root zone soil water [-]
    real    :: WSTRES    !water stress coeficient [-]  (1. for wilting )
    real    :: LAPM      !leaf area per unit mass [m2/g]

    if ((lutyp == ISWATER) .or. (lutyp == ISBARREN) .or. (lutyp == ISICE) .or. (lutyp == ISURBAN)) then
       XLAI   = 0.0
       XSAI   = 0.0
       GPP    = 0.0
       NPP    = 0.0
       NEE    = 0.0
       AUTORS = 0.0
       HETERS = 0.0
       TOTSC  = 0.0
       TOTLB  = 0.0
       LFMASS = 0.0
       RTMASS = 0.0
       STMASS = 0.0
       WOOD   = 0.0
       STBLCP = 0.0
       FASTCP = 0.0

       return
    end if

    LAPM = LK_SLA(lutyp) / 1000.0   ! m2/kg -> m2/g

    ! water stress

    WSTRES  = 1.0 - BTRAN

    WROOT  = 0.0
    do J = 1, LK_NROOT(lutyp)
       WROOT = WROOT + SMC(J) / SMCMAX *  DZSNSO(J) / (-ZSOIL(LK_NROOT(lutyp)))
    end do

    call co2flux(NSNOW  ,NSOIL  ,lutyp ,IGS    ,DT     , & !in
         DZSNSO ,STC    ,PSN    ,TROOT  ,TV     , & !in
         WROOT  ,WSTRES ,FOLN   ,LAPM   ,         & !in
         LAT    ,ILOC   ,JLOC   ,fveg   ,         & !in
         XLAI   ,XSAI   ,LFMASS ,RTMASS ,STMASS , & !inout
         FASTCP ,STBLCP ,WOOD   ,                 & !inout
         GPP    ,NPP    ,NEE    ,AUTORS ,HETERS , & !out
         TOTSC  ,TOTLB  )                           !out

    !   CALL BVOC (VOCFLX,  lutyp,  VEGFAC,   APAR,   TV)
    !   CALL CH4
  end subroutine carbon


  subroutine co2flux(NSNOW  ,NSOIL  ,lutyp ,IGS    ,DT     , & !in
       DZSNSO ,STC    ,PSN    ,TROOT  ,TV     , & !in
       WROOT  ,WSTRES ,FOLN   ,LAPM   ,         & !in
       LAT    ,ILOC   ,JLOC   ,fveg   ,         & !in
       XLAI   ,XSAI   ,LFMASS ,RTMASS ,STMASS , & !inout
       FASTCP ,STBLCP ,WOOD   ,                 & !inout
       GPP    ,NPP    ,NEE    ,AUTORS ,HETERS , & !out
       TOTSC  ,TOTLB  )                           !out
    ! The original code is from RE Dickinson et al.(1998), modifed by Guo-Yue Niu, 2004
    use noahmp_veg_param, only: ISEGBLF
    use noahmp_veg_param, only: LK_DILEFC
    use noahmp_veg_param, only: LK_DILEFW
    use noahmp_veg_param, only: LK_FRAGR
    use noahmp_veg_param, only: LK_LTOVRC
    use noahmp_veg_param, only: LK_WRRAT
    use noahmp_veg_param, only: LK_WDPOOL
    use noahmp_veg_param, only: LK_TDLEF
    use noahmp_veg_param, only: LK_RMF25
    use noahmp_veg_param, only: LK_RMS25
    use noahmp_veg_param, only: LK_RMR25
    use noahmp_veg_param, only: LK_FOLNMX
    use noahmp_veg_param, only: LK_TMIN
    use noahmp_veg_param, only: LK_ARM
    use noahmp_veg_param, only: LK_MRP
    implicit none
    ! input
    integer                        , intent(in) :: ILOC   !grid index
    integer                        , intent(in) :: JLOC   !grid index
    integer                        , intent(in) :: lutyp !vegetation physiology type
    integer                        , intent(in) :: NSNOW  !number of snow layers
    integer                        , intent(in) :: NSOIL  !number of soil layers
    real                           , intent(in) :: DT     !time step (s)
    real                           , intent(in) :: LAT    !latitude (radians)
    real                           , intent(in) :: IGS    !growing season index (0=off, 1=on)
    real, dimension(-NSNOW+1:NSOIL), intent(in) :: DZSNSO !snow/soil layer thickness [m]
    real, dimension(-NSNOW+1:NSOIL), intent(in) :: STC    !snow/soil temperature [k]
    real                           , intent(in) :: PSN    !total leaf photosynthesis (umolco2/m2/s)
    real                           , intent(in) :: TROOT  !root-zone averaged temperature (k)
    real                           , intent(in) :: TV     !leaf temperature (k)
    real                           , intent(in) :: WROOT  !root zone soil water
    real                           , intent(in) :: WSTRES !soil water stress
    real                           , intent(in) :: FOLN   !foliage nitrogen (%)
    real                           , intent(in) :: LAPM   !leaf area per unit mass [m2/g]
    real                           , intent(in) :: fveg   !vegetation greenness fraction

    ! input and output

    real                        , intent(inout) :: XLAI   !leaf  area index from leaf carbon [-]
    real                        , intent(inout) :: XSAI   !stem area index from leaf carbon [-]
    real                        , intent(inout) :: LFMASS !leaf mass [g/m2]
    real                        , intent(inout) :: RTMASS !mass of fine roots [g/m2]
    real                        , intent(inout) :: STMASS !stem mass [g/m2]
    real                        , intent(inout) :: FASTCP !short lived carbon [g/m2]
    real                        , intent(inout) :: STBLCP !stable carbon pool [g/m2]
    real                        , intent(inout) :: WOOD   !mass of wood (incl. woody roots) [g/m2]

    ! output

    real                          , intent(out) :: GPP    !net instantaneous assimilation [g/m2/s]
    real                          , intent(out) :: NPP    !net primary productivity [g/m2]
    real                          , intent(out) :: NEE    !net ecosystem exchange (autors+heters-gpp)
    real                          , intent(out) :: AUTORS !net ecosystem resp. (maintance and growth)
    real                          , intent(out) :: HETERS !organic respiration
    real                          , intent(out) :: TOTSC  !total soil carbon (g/m2)
    real                          , intent(out) :: TOTLB  !total living carbon (g/m2)

    ! local

    real                   :: CFLUX    !carbon flux to atmosphere [g/m2/s]
    real                   :: LFMSMN   !minimum leaf mass [g/m2]
    real                   :: RSWOOD   !wood respiration [g/m2]
    real                   :: RSLEAF   !leaf maintenance respiration per timestep [g/m2]
    real                   :: RSROOT   !fine root respiration per time step [g/m2]
    real                   :: NPPL     !leaf net primary productivity [g/m2/s]
    real                   :: NPPR     !root net primary productivity [g/m2/s]
    real                   :: NPPW     !wood net primary productivity [g/m2/s]
    real                   :: NPPS     !wood net primary productivity [g/m2/s]
    real                   :: DIELF    !death of leaf mass per time step [g/m2]

    real                   :: ADDNPPLF !leaf assimil after resp. losses removed [g/m2]
    real                   :: ADDNPPST !stem assimil after resp. losses removed [g/m2]
    real                   :: CARBFX   !carbon assimilated per model step [g/m2]
    real                   :: GRLEAF   !growth respiration rate for leaf [g/m2/s]
    real                   :: GRROOT   !growth respiration rate for root [g/m2/s]
    real                   :: GRWOOD   !growth respiration rate for wood [g/m2/s]
    real                   :: GRSTEM   !growth respiration rate for stem [g/m2/s]
    real                   :: LEAFPT   !fraction of carbon allocated to leaves [-]
    real                   :: LFDEL    !maximum  leaf mass  available to change [g/m2/s]
    real                   :: LFTOVR   !stem turnover per time step [g/m2]
    real                   :: STTOVR   !stem turnover per time step [g/m2]
    real                   :: WDTOVR   !wood turnover per time step [g/m2]
    real                   :: RSSOIL   !soil respiration per time step [g/m2]
    real                   :: RTTOVR   !root carbon loss per time step by turnover [g/m2]
    real                   :: STABLC   !decay rate of fast carbon to slow carbon [g/m2/s]
    real                   :: WOODF    !calculated wood to root ratio [-]
    real                   :: NONLEF   !fraction of carbon to root and wood [-]
    real                   :: ROOTPT   !fraction of carbon flux to roots [-]
    real                   :: WOODPT   !fraction of carbon flux to wood [-]
    real                   :: STEMPT   !fraction of carbon flux to stem [-]
    real                   :: RESP     !leaf respiration [umol/m2/s]
    real                   :: RSSTEM   !stem respiration [g/m2/s]

    real                   :: FSW      !soil water factor for microbial respiration
    real                   :: FST      !soil temperature factor for microbial respiration
    real                   :: FNF      !foliage nitrogen adjustemt to respiration (<= 1)
    real                   :: TF       !temperature factor
    real                   :: RF       !respiration reduction factor (<= 1)
    real                   :: STDEL
    real                   :: STMSMN
    real                   :: SAPM     !stem area per unit mass (m2/g)
    real                   :: DIEST

    ! constants
    real                   :: BF       !parameter for present wood allocation [-]
    real                   :: RSWOODC  !wood respiration coeficient [1/s]
    real                   :: STOVRC   !stem turnover coefficient [1/s]
    real                   :: RSDRYC   !degree of drying that reduces soil respiration [-]
    real                   :: RTOVRC   !root turnover coefficient [1/s]
    real                   :: WSTRC    !water stress coeficient [-]
    real                   :: LAIMIN   !minimum leaf area index [m2/m2]
    real                   :: XSAMIN   !minimum leaf area index [m2/m2]
    real                   :: SC
    real                   :: SD
    real                   :: VEGFRAC

    ! Respiration as a function of temperature

    real :: r, x
    r(x) = exp(0.08 * (x - 298.16))

    ! constants
    RTOVRC  = 2.0E-8        !original was 2.0e-8
    RSDRYC  = 40.0          !original was 40.0
    RSWOODC = 3.0E-10       !
    BF      = 0.90          !original was 0.90   ! carbon to roots
    WSTRC   = 100.0
    LAIMIN  = 0.05
    XSAMIN  = 0.01

    SAPM    = 3.0 * 0.001      ! m2/kg -->m2/g
    LFMSMN  = laimin / lapm
    STMSMN  = xsamin / sapm

    ! respiration

    if (IGS == 0.0) then
       RF = 0.5
    else
       RF = 1.0
    end if

    FNF     = min(FOLN / max(1.0E-06, LK_FOLNMX(lutyp)), 1.0)
    TF      = LK_ARM(lutyp) ** ((TV - 298.16) / 10.0)
    RESP    = LK_RMF25(lutyp) * TF * FNF * XLAI * RF * (1.0 - WSTRES) ! umol/m2/s
    RSLEAF  = min(LFMASS / DT, RESP * 12.0E-6)                         ! g/m2/s

    RSROOT  = LK_RMR25(lutyp) * (RTMASS * 1.0E-3) * TF * RF * 12.0E-6         ! g/m2/s
    RSSTEM  = LK_RMS25(lutyp) * (STMASS * 1.0E-3) * TF * RF * 12.0E-6         ! g/m2/s
    RSWOOD  = RSWOODC * R(TV) * WOOD * LK_WDPOOL(lutyp)

    ! carbon assimilation
    ! 1 mole -> 12 g carbon or 44 g CO2; 1 umol -> 12.e-6 g carbon;

    CARBFX  = PSN * 12.0E-6              ! umol co2 /m2/ s -> g/m2/s carbon

    ! fraction of carbon into leaf versus nonleaf

    LEAFPT = exp(0.01 * (1.0 - exp(0.75 * XLAI)) * XLAI)
    if (lutyp == ISEGBLF) LEAFPT = exp(0.01 * (1.0 - exp(0.50 * XLAI)) * XLAI)

    NONLEF = 1.0 - LEAFPT
    STEMPT = XLAI / 10.0
    LEAFPT = LEAFPT - STEMPT

    !  fraction of carbon into wood versus root

    if (WOOD > 0.0) then
       WOODF = (1.0 - exp(-BF * (LK_WRRAT(lutyp) * RTMASS / WOOD)) / BF) * LK_WDPOOL(lutyp)
    else
       WOODF = 0.0
    end if

    ROOTPT = NONLEF * (1.0 - WOODF)
    WOODPT = NONLEF * WOODF

    ! leaf and root turnover per time step

    LFTOVR = LK_LTOVRC(lutyp) * 1.0E-6 * LFMASS
    STTOVR = LK_LTOVRC(lutyp) * 1.0E-6 * STMASS
    RTTOVR = RTOVRC * RTMASS
    WDTOVR = 9.5E-10 * WOOD

    ! seasonal leaf die rate dependent on temp and water stress
    ! water stress is set to 1 at permanent wilting point

    SC  = exp(-0.3 * max(0.0, TV - LK_TDLEF(lutyp))) * (LFMASS / 120.0)
    SD  = exp((WSTRES - 1.0) * WSTRC)
    DIELF = LFMASS * 1.0E-6 * (LK_DILEFW(lutyp) * SD + LK_DILEFC(lutyp) * SC)
    DIEST = STMASS * 1.0E-6 * (LK_DILEFW(lutyp) * SD + LK_DILEFC(lutyp) * SC)

    ! calculate growth respiration for leaf, rtmass and wood

    GRLEAF = max(0.0, LK_FRAGR(lutyp) * (LEAFPT * CARBFX - RSLEAF))
    GRSTEM = max(0.0, LK_FRAGR(lutyp) * (STEMPT * CARBFX - RSSTEM))
    GRROOT = max(0.0, LK_FRAGR(lutyp) * (ROOTPT * CARBFX - RSROOT))
    GRWOOD = max(0.0, LK_FRAGR(lutyp) * (WOODPT * CARBFX - RSWOOD))

    ! Impose lower T limit for photosynthesis

    ADDNPPLF = max(0.0, LEAFPT * CARBFX - GRLEAF - RSLEAF)
    ADDNPPST = max(0.0, STEMPT * CARBFX - GRSTEM - RSSTEM)
    if (TV < LK_TMIN(lutyp)) ADDNPPLF = 0.0
    if (TV < LK_TMIN(lutyp)) ADDNPPST = 0.0

    ! update leaf, root, and wood carbon
    ! avoid reducing leaf mass below its minimum value but conserve mass

    LFDEL = (LFMASS - LFMSMN) / DT
    STDEL = (STMASS - STMSMN) / DT
    DIELF = min(DIELF, LFDEL + ADDNPPLF - LFTOVR)
    DIEST = min(DIEST, STDEL + ADDNPPST - STTOVR)

    ! net primary productivities

    NPPL   = max(ADDNPPLF, -LFDEL)
    NPPS   = max(ADDNPPST, -STDEL)
    NPPR   = ROOTPT * CARBFX - RSROOT - GRROOT
    NPPW   = WOODPT * CARBFX - RSWOOD - GRWOOD

    ! masses of plant components

    LFMASS = LFMASS + (NPPL - LFTOVR - DIELF) * DT
    STMASS = STMASS + (NPPS - STTOVR - DIEST) * DT   ! g/m2
    RTMASS = RTMASS + (NPPR - RTTOVR) * DT

    if (RTMASS < 0.0) then
       RTTOVR = NPPR
       RTMASS = 0.0
    end if
    WOOD = (WOOD + (NPPW - WDTOVR) * DT) * LK_WDPOOL(lutyp)

    ! soil carbon budgets

    FASTCP = FASTCP + (RTTOVR + LFTOVR + STTOVR + WDTOVR + DIELF) * DT

    FST = 2.0 ** ((STC(1) - 283.16) / 10.0)
    FSW = WROOT / (0.20 + WROOT) * 0.23 / (0.23 + WROOT)
    RSSOIL = FSW * FST * LK_MRP(lutyp) * max(0.0, FASTCP * 1.0E-3) * 12.0E-6

    STABLC = 0.1 * RSSOIL
    FASTCP = FASTCP - (RSSOIL + STABLC) * DT
    STBLCP = STBLCP + STABLC * DT

    !  total carbon flux
    CFLUX  = -CARBFX + RSLEAF + RSROOT + RSWOOD + RSSTEM &
         + RSSOIL + GRLEAF + GRROOT + GRWOOD                ! g/m2/s

    ! for outputs
    GPP    = CARBFX                                             !g/m2/s C
    NPP    = NPPL + NPPW + NPPR                                 !g/m2/s C
    AUTORS = RSROOT + RSWOOD  + RSLEAF +  &                     !g/m2/s C
         &  GRLEAF + GRROOT + GRWOOD                           !g/m2/s C
    HETERS = RSSOIL                                             !g/m2/s C
    NEE    = (AUTORS + HETERS - GPP) * 44.0 / 12.0              !g/m2/s CO2
    TOTSC  = FASTCP + STBLCP                                    !g/m2   C
    TOTLB  = LFMASS + RTMASS + WOOD                             !g/m2   C

    ! leaf area index and stem area index
    XLAI = max(LFMASS * LAPM, LAIMIN)
    XSAI = max(STMASS * SAPM, XSAMIN)

  end subroutine co2flux


  subroutine bvocflux(VOCFLX,  lutyp,  VEGFRAC,  APAR,   TV )
    use noahmp_const, only: RGAS
    use noahmp_veg_param, only: LK_SLAREA
    use noahmp_veg_param, only: LK_EPS
    implicit none
    ! source file:       BVOC
    ! purpose:           BVOC emissions
    ! DESCRIPTION:
    ! Volatile organic compound emission
    ! This code simulates volatile organic compound emissions
    ! following the algorithm presented in Guenther, A., 1999: Modeling
    ! Biogenic Volatile Organic Compound Emissions to the Atmosphere. In
    ! Reactive Hydrocarbons in the Atmosphere, Ch. 3
    ! This model relies on the assumption that 90% of isoprene and monoterpene
    ! emissions originate from canopy foliage:
    !    E = epsilon * gamma * density * delta
    ! The factor delta (longterm activity factor) applies to isoprene emission
    ! from deciduous plants only. We neglect this factor at the present time.
    ! This factor is discussed in Guenther (1997).
    ! Subroutine written to operate at the patch level.
    ! IN FINAL IMPLEMENTATION, REMEMBER:
    ! 1. may wish to call this routine only as freq. as rad. calculations
    ! 2. may wish to place epsilon values directly in pft-physiology file
    ! input
    integer                     ,intent(in) :: lutyp  !vegetation type
    real                        ,intent(in) :: vegfrac !green vegetation fraction [0.0-1.0]
    real                        ,intent(in) :: apar    !photosynthesis active energy by canopy (w/m2)
    real                        ,intent(in) :: tv      !vegetation canopy temperature (k)

    ! output
    real                        ,intent(out) :: vocflx(5) ! voc fluxes [ug C m-2 h-1]

    ! Locals
    real, parameter :: alpha  = 0.0027   ! empirical coefficient
    real, parameter :: cl1    = 1.066    ! empirical coefficient
    real, parameter :: ct1    = 95000.0  ! empirical coefficient [J mol-1]
    real, parameter :: ct2    = 230000.0 ! empirical coefficient [J mol-1]
    real, parameter :: ct3    = 0.961    ! empirical coefficient
    real, parameter :: tm     = 314.0    ! empirical coefficient [K]
    real, parameter :: tstd   = 303.0    ! std temperature [K]
    real, parameter :: bet    = 0.09     ! beta empirical coefficient [K-1]

    integer :: ivoc        ! do-loop index
    integer :: ityp        ! do-loop index
    real :: epsilon(5)
    real :: gamma(5)
    real :: density
    real :: elai
    real :: par, cl, reciprod, ct

    ! epsilon :

    do ivoc = 1, 5
       epsilon(ivoc) = LK_EPS(ivoc,lutyp)
    end do

    ! gamma : Activity factor. Units [dimensionless]

    reciprod = 1.0 / (RGAS * tv * tstd)
    ct = exp(ct1 * (tv - tstd) * reciprod) / &
         (ct3 + exp(ct2 * (tv - tm) * reciprod))

    par = apar * 4.6 ! (multiply w/m2 by 4.6 to get umol/m2/s)
    cl  = alpha * cl1 * par * (1. + alpha * alpha * par * par)**(-0.5)

    gamma(1) = cl * ct ! for isoprenes

    do ivoc = 2, 5
       gamma(ivoc) = exp(bet * (tv - tstd))
    end do

    ! Foliage density

    ! transform vegfrac to lai

    elai    = max(0.0, -6.5 / 2.5 * alog((1.0 - vegfrac)))
    density = elai / (LK_SLAREA(lutyp) * 0.5)

    ! calculate the voc flux

    do ivoc = 1, 5
       vocflx(ivoc) = epsilon(ivoc) * gamma(ivoc) * density
    end do

  end subroutine bvocflux
end module noahmp_func