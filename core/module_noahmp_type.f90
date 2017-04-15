module noahmp_type
  use noahmp_const, only: r4
  use noahmp_const, only: nan4
  use noahmp_global, only: nband
  use noahmp_global, only: nsoil
  use noahmp_global, only: msnow
  implicit none


  type noahmp_state_t
     sequence
     ! geolocation
     real(r4) :: lat  !latitude [radius]
     real(r4) :: lon  !longitude [radius]
     ! grid shape
     real(r4) :: zlvl         !reference height [m]
     real(r4) :: zsoil(nsoil) !soil-layer-bottom depth [m,+:up]
     real(r4) :: zsnow(msnow) !snow-layer-top height [m,+:up]
     real(r4) :: nsnow        !number of snow layers [1]
     ! plant and soil geography
     integer :: lutyp  !vegetation type
     integer :: sltyp !soil type
     ! phenology
     real(r4) :: lai  !leaf area index [m2 m-2]
     real(r4) :: sai  !stem area index [m2 m-2]
     real(r4) :: fveg !green vegetation fraction [m2 m-2,0~1]
     ! canopy temperature and storage
     real(r4) :: cantmp !canopy temperature [K]
     real(r4) :: canwat !canopy-intercepted liquid water [kg m-2]
     real(r4) :: cansno !canopy-intercepted solid water [kg m-2]
     ! snow temperature and storage
     real(r4) :: snowtmp(msnow) !snow temperature [K]
     real(r4) :: snowwat(msnow) !snow liquid water mass [kg m-3]
     real(r4) :: snowice(msnow) !snow solid water mass [kg m-3]
     ! soil temperature and storage
     real(r4) :: soiltmp(nsoil) !soil temperature [K]
     real(r4) :: soilwat(nsoil) !soil liquid water content [m3 m-3]
     real(r4) :: soilice(nsoil) !soil solid water content [m3 m-3]
     ! groundwater storage
     real(r4) :: grndwat  !groundwater storage [kg m-2]
     real(r4) :: zwt      !groundwater table depth [m,+:up]
  end type noahmp_state_t


  type noahmp_atm_t
     sequence

     real(r4) :: cosz
     real(r4) :: srdnd(nband)
     real(r4) :: srdni(nband)
     real(r4) :: lrdn
     real(r4) :: uatm
     real(r4) :: vatm
     real(r4) :: watm
     real(r4) :: tatm
     real(r4) :: ptatm
     real(r4) :: vpatm
     real(r4) :: rhatm
     real(r4) :: pres

     real(r4) :: fprcp
     real(r4) :: prcpwatcc
     real(r4) :: prcpsnocc
     real(r4) :: prcpwatcl
     real(r4) :: prcpsnocl
     real(r4) :: prcpsnobd
  end type noahmp_atm_t


  type noahmp_flux_t
     sequence

     real(r4) :: srabc(nband)
     real(r4) :: srabg(nband)
     real(r4) :: lrabc
     real(r4) :: lrabg

     real(r4) :: shc
     real(r4) :: shg
     real(r4) :: lhc
     real(r4) :: lhg

     real(r4) :: evapc
     real(r4) :: evapg
     real(r4) :: transun
     real(r4) :: transha
     real(r4) :: tranrt(nsoil)

     real(r4) :: qthruwat
     real(r4) :: qthrusno
     real(r4) :: qintcwat
     real(r4) :: qintcsno
     real(r4) :: qdripwat
     real(r4) :: qdripsno
     real(r4) :: qsrfwat
     real(r4) :: qsrfsno

     real(r4) :: runsrf
     real(r4) :: runsub

     real(r4) :: qdisch
     real(r4) :: qrech
  end type noahmp_flux_t


  type noahmp_prop_t
     sequence

     real(r4) :: albv(nband)
     real(r4) :: albg(nband)
     real(r4) :: emsv
     real(r4) :: emsg

     real(r4) :: canwatmax
     real(r4) :: cansnomax
  end type noahmp_prop_t

end module noahmp_type
