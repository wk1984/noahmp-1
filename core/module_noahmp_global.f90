! Noah-MP Global Variables
! These variables should be read-only.
module noahmp_global
  use noahmp_const, only: r4
  implicit none

  !%% Noah-MP General Settings
  ! number of solar radiation wave bands
  integer, parameter :: nband = 2 !1=vis, 2=nir
  ! number of soil layers
  integer, parameter :: nsoil = 4
  ! maximum number of snow layers
  integer, parameter :: msnow = 3

  !%% Noah-MP Options
  ! options for dynamic vegetation
  integer :: opt_dveg    != 4   !
  ! 1 -> off (use table LAI; use fveg = SHDFAC from input)
  ! 2 -> on (together with OPT_CRS = 1)
  ! 3 -> off (use table LAI; calculate fveg)
  ! 4 -> off (use table LAI; use maximum vegetation fraction)

  ! options for canopy stomatal resistance
  integer :: opt_crs != 1    !(must 1 when OPT_DVEG = 2)
  ! 1-> Ball-Berry; 2->Jarvis

  ! options for soil moisture factor for stomatal resistance
  integer :: opt_btr != 1    !(suggested 1)
  ! 1-> Noah (soil moisture)
  ! 2-> CLM  (matric potential)
  ! 3-> SSiB (matric potential)

  ! options for runoff and groundwater
  integer :: opt_run != 1    !(suggested 1)
  ! 1 -> TOPMODEL with groundwater (Niu et al. 2007 JGR) ;
  ! 2 -> TOPMODEL with an equilibrium water table (Niu et al. 2005 JGR) ;
  ! 3 -> original surface and subsurface runoff (free drainage)
  ! 4 -> BATS surface and subsurface runoff (free drainage)

  ! options for surface layer drag coeff (CH & CM)
  integer :: opt_sfc != 1    !(1 or 2)
  ! 1->M-O ; 2->original Noah (Chen97)

  ! options for supercooled liquid water (or ice fraction)
  integer :: opt_frz != 1    !(1 or 2)
  ! 1-> no iteration (Niu and Yang, 2006 JHM); 2: Koren's iteration

  ! options for frozen soil permeability
  integer :: opt_inf != 1    !(suggested 1)
  ! 1 -> linear effects, more permeable (Niu and Yang, 2006, JHM)
  ! 2 -> nonlinear effects, less permeable (old)

  ! options for radiation transfer
  integer :: opt_rad != 1    !(suggested 1)
  ! 1 -> modified two-stream (gap = F(solar angle, 3D structure ...)<1-fveg)
  ! 2 -> two-stream applied to grid-cell (gap = 0)
  ! 3 -> two-stream applied to vegetated fraction (gap=1-fveg)

  ! options for ground snow surface albedo
  integer :: opt_alb != 2    !(suggested 2)
  ! 1-> BATS; 2 -> CLASS

  ! options for partitioning  precipitation into rainfall & snowfall
  integer :: opt_snf != 1    !(suggested 1)
  ! 1 -> Jordan (1991); 2 -> BATS: when SFCTMP<TFRZ+2.2 ; 3-> SFCTMP<TFRZ

  ! options for lower boundary condition of soil temperature
  integer :: opt_tbot != 2   !(suggested 2)
  ! 1 -> zero heat flux from bottom (ZBOT and TBOT not used)
  ! 2 -> TBOT at ZBOT (8m) read from a file (original Noah)

  ! options for snow/soil temperature time scheme (only layer 1)
  integer :: opt_stc != 1    !(suggested 1)
  ! 1 -> semi-implicit; 2 -> full implicit (original Noah)

contains
  subroutine noahmp_set_options(iopt_dveg, iopt_crs, iopt_btr, &
       & iopt_run, iopt_sfc, iopt_frz, &
       & iopt_inf, iopt_rad, iopt_alb, &
       & iopt_snf, iopt_tbot, iopt_stc)
    implicit none

    integer,  intent(in) :: iopt_dveg !dynamic vegetation (1 -> off ; 2 -> on) with opt_crs = 1
    integer,  intent(in) :: iopt_crs  !canopy stomatal resistance (1-> Ball-Berry; 2->Jarvis)
    integer,  intent(in) :: iopt_btr  !soil moisture factor for stomatal resistance (1-> Noah; 2-> CLM; 3-> SSiB)
    integer,  intent(in) :: iopt_run  !runoff and groundwater (1->SIMGM; 2->SIMTOP; 3->Schaake96; 4->BATS)
    integer,  intent(in) :: iopt_sfc  !surface layer drag coeff (CH & CM) (1->M-O; 2->Chen97)
    integer,  intent(in) :: iopt_frz  !supercooled liquid water (1-> NY06; 2->Koren99)
    integer,  intent(in) :: iopt_inf  !frozen soil permeability (1-> NY06; 2->Koren99)
    integer,  intent(in) :: iopt_rad  !radiation transfer (1->gap=F(3D,cosz); 2->gap=0; 3->gap=1-fveg)
    integer,  intent(in) :: iopt_alb  !snow surface albedo (1->BATS; 2->CLASS)
    integer,  intent(in) :: iopt_snf  !rainfall & snowfall (1-Jordan91; 2->BATS; 3->Noah)
    integer,  intent(in) :: iopt_tbot !lower boundary of soil temperature (1->zero-flux; 2->Noah)

    integer,  intent(in) :: iopt_stc  !snow/soil temperature time scheme (only layer 1)
    ! 1 -> semi-implicit; 2 -> full implicit (original Noah)

    opt_dveg = iopt_dveg

    opt_crs  = iopt_crs
    opt_btr  = iopt_btr
    opt_run  = iopt_run
    opt_sfc  = iopt_sfc
    opt_frz  = iopt_frz
    opt_inf  = iopt_inf
    opt_rad  = iopt_rad
    opt_alb  = iopt_alb
    opt_snf  = iopt_snf
    opt_tbot = iopt_tbot
    opt_stc  = iopt_stc

  end subroutine noahmp_set_options
end module noahmp_global