module noahmp_veg_param
  use noahmp_const, only: r4
  use noahmp_const, only: nan4
  use noahmp_global, only: nband
  implicit none

  ! maximum number of land use type
  integer, parameter :: MLUTYP = 27

  character(256) :: lutyp_scheme = ' '
  integer :: nlutyp = 0

  integer :: ISURBAN  = -1
  integer :: ISWATER  = -1
  integer :: ISBARREN = -1
  integer :: ISICE   = -1
  integer :: ISEGBLF  = -1

  real(r4) :: LK_XL(MLUTYP)          = nan4 !leaf/stem orientation index
  real(r4) :: LK_RHOL(nband,MLUTYP)  = nan4 !leaf reflectance: 1=vis, 2=nir
  real(r4) :: LK_RHOS(nband,MLUTYP)  = nan4 !stem reflectance: 1=vis, 2=nir
  real(r4) :: LK_TAUL(nband,MLUTYP)  = nan4 !leaf transmittance: 1=vis, 2=nir
  real(r4) :: LK_TAUS(nband,MLUTYP)  = nan4 !stem transmittance: 1=vis, 2=nir

  integer :: LK_NROOT(MLUTYP)    =   -1!rooting depth [as the number of layers]
  real(r4) :: LK_CANWMXP(MLUTYP) = nan4!maximum intercepted h2o per unit lai+sai (mm)
  real(r4) :: LK_DLEAF(MLUTYP)   = nan4!characteristic leaf dimension (m)
  real(r4) :: LK_Z0MVT(MLUTYP)   = nan4!momentum roughness length (m)
  real(r4) :: LK_HVT(MLUTYP)     = nan4!top of canopy (m)
  real(r4) :: LK_HVB(MLUTYP)     = nan4!bottom of canopy (m)
  real(r4) :: LK_DEN(MLUTYP)     = nan4!tree density (no. of trunks per m2) [???useless]
  real(r4) :: LK_RCROWN(MLUTYP)  = nan4!tree crown radius (m)
  real(r4) :: LK_CWPVT(MLUTYP)   = nan4!empirical canopy wind parameter

  real(r4) :: LK_SAI12M(12,MLUTYP) = nan4!monthly stem area index, one-sided
  real(r4) :: LK_LAI12M(12,MLUTYP) = nan4!monthly leaf area index, one-sided

  real(r4) :: LK_SLA(MLUTYP)     = nan4!single-side leaf area per Kg [m2/kg]

  real(r4) :: LK_DILEFC(MLUTYP)  = nan4!coeficient for leaf stress death [1/s]
  real(r4) :: LK_DILEFW(MLUTYP)  = nan4!coeficient for leaf stress death [1/s]
  real(r4) :: LK_FRAGR(MLUTYP)   = nan4!fraction of growth respiration  !original was 0.3
  real(r4) :: LK_LTOVRC(MLUTYP)  = nan4!leaf turnover [1/s]
  real(r4) :: LK_WRRAT(MLUTYP)   = nan4 !wood to non-wood ratio
  real(r4) :: LK_WDPOOL(MLUTYP)  = nan4 !wood pool (switch 1 or 0) depending on woody or not [-]
  real(r4) :: LK_TDLEF(MLUTYP)   = nan4 !characteristic T for leaf freezing [K]

  real(r4) :: LK_RGL(MLUTYP)     = nan4!parameter used in radiation stress function
  real(r4) :: LK_HS(MLUTYP)      = nan4!parameter used in vapor pressure deficit function
  real(r4) :: LK_RSMAX(MLUTYP)   = nan4!maximum stomatal resistance
  real(r4) :: LK_RSMIN(MLUTYP)   = nan4!minimum Canopy Resistance [s/m]
  real(r4) :: LK_TOPT(MLUTYP)    = nan4!optimum transpiration air temperature.

  integer  :: LK_C3C4(MLUTYP)    =    0!photosynthetic pathway: 1 = c3, 2 = c4
  real(r4) :: LK_KC25(MLUTYP)    = nan4!co2 michaelis-menten constant at 25c (pa)
  real(r4) :: LK_AKC(MLUTYP)     = nan4!q10 for kc25
  real(r4) :: LK_KO25(MLUTYP)    = nan4!o2 michaelis-menten constant at 25c (pa)
  real(r4) :: LK_AKO(MLUTYP)     = nan4!q10 for ko25
  real(r4) :: LK_VCMX25(MLUTYP)  = nan4!maximum rate of carboxylation at 25c (umol co2/m**2/s)
  real(r4) :: LK_AVCMX(MLUTYP)   = nan4!q10 for vcmx25
  real(r4) :: LK_BP(MLUTYP)      = nan4!minimum leaf conductance (umol/m**2/s)
  real(r4) :: LK_MP(MLUTYP)      = nan4!slope of conductance-to-photosynthesis relationship
  real(r4) :: LK_QE25(MLUTYP)    = nan4!quantum efficiency at 25c (umol co2 / umol photon)
  real(r4) :: LK_AQE(MLUTYP)     = nan4!q10 for qe25 [???useless]
  real(r4) :: LK_FOLNMX(MLUTYP)  = nan4!foliage nitrogen concentration when f(n)=1 (%)
  real(r4) :: LK_TMIN(MLUTYP)    = nan4!minimum temperature for photosynthesis (k)
  real(r4) :: LK_RMF25(MLUTYP)   = nan4!leaf maintenance respiration at 25c (umol co2/m**2/s)
  real(r4) :: LK_RMS25(MLUTYP)   = nan4!stem maintenance respiration at 25c (umol co2/kg bio/s)
  real(r4) :: LK_RMR25(MLUTYP)   = nan4!root maintenance respiration at 25c (umol co2/kg bio/s)
  real(r4) :: LK_ARM(MLUTYP)     = nan4!q10 for maintenance respiration
  real(r4) :: LK_MRP(MLUTYP)     = nan4!microbial respiration parameter (umol co2 /kg c/ s)

  real(r4) :: LK_SLAREA(MLUTYP)  = nan4
  real(r4) :: LK_EPS(5,MLUTYP)   = nan4

contains
  subroutine noahmp_veg_param_readptable(tag)
    use noahmp_utils, only: assert
    use noahmp_utils, only: noahmp_ptable_read_var
    use noahmp_utils, only: noahmp_ptable_read_tablestr
    implicit none
    integer, parameter :: LINELEN = 512
    character(256), parameter :: veg_table_filename = 'VEGPARMMP.TBL'
    character(*), intent(in) :: tag
    character(LINELEN) :: tblstr(MLUTYP)
    integer :: ii, jj, index
    integer :: ierr

    call noahmp_ptable_read_var(veg_table_filename, 'ISURBAN', tag, ISURBAN)
    call noahmp_ptable_read_var(veg_table_filename, 'ISWATER', tag, ISWATER)
    call noahmp_ptable_read_var(veg_table_filename, 'ISBARREN', tag, ISBARREN)
    call noahmp_ptable_read_var(veg_table_filename, 'ISICE', tag, ISICE)
    call noahmp_ptable_read_var(veg_table_filename, 'ISEGBLF', tag, ISEGBLF)

    ! radiation parameters
    call noahmp_ptable_read_tablestr(veg_table_filename, 'RAD', tag, nlutyp, tblstr)
    call assert(nlutyp <= MLUTYP, 'Table RAD#'//trim(tag)//' is too big to be read')
    do ii = 1, nlutyp
       read(tblstr(ii), *, iostat=ierr) index, &
            & LK_XL(ii), (LK_RHOL(jj,ii), jj=1,nband), (LK_RHOS(jj,ii), jj=1,nband), &
            & (LK_TAUL(jj,ii), jj=1,nband), (LK_TAUS(jj,ii), jj=1,nband)
       call assert(ierr == 0, 'Unable to read RAD#'//trim(tag)//' from file '//trim(veg_table_filename))
    end do

    ! 12-month LAI
    call noahmp_ptable_read_tablestr(veg_table_filename, 'LAI12M', tag, nlutyp, tblstr)
    call assert(nlutyp <= MLUTYP, 'Table LAI12M#'//trim(tag)//' is too big to be read')
    do ii = 1, nlutyp
       read(tblstr(ii), *, iostat=ierr) index, (LK_LAI12M(jj,ii), jj=1,12)
       call assert(ierr == 0, 'Unable to read LAI12M#'//trim(tag)//' from file '//trim(veg_table_filename))
    end do

    ! 12-month SAI
    call noahmp_ptable_read_tablestr(veg_table_filename, 'SAI12M', tag, nlutyp, tblstr)
    call assert(nlutyp <= MLUTYP, 'Table SAI12M#'//trim(tag)//' is too big to be read')
    do ii = 1, nlutyp
       read(tblstr(ii), *, iostat=ierr) index, (LK_SAI12M(jj,ii), jj=1,12)
       call assert(ierr == 0, 'Unable to read SAI12M#'//trim(tag)//' from file '//trim(veg_table_filename))
    end do

    ! dynamic vegetation parameters
    call noahmp_ptable_read_tablestr(veg_table_filename, 'DVEG', tag, nlutyp, tblstr)
    call assert(nlutyp <= MLUTYP, 'Table DVEG#'//trim(tag)//' is too big to be read')
    do ii = 1, nlutyp
       read(tblstr(ii), *, iostat=ierr) index, &
            & LK_SLA(ii), LK_DILEFC(ii), LK_DILEFW(ii), LK_FRAGR(ii), LK_LTOVRC(ii), &
            & LK_WRRAT(ii), LK_WDPOOL(ii), LK_TDLEF(ii)
       call assert(ierr == 0, 'Unable to read DVEG#'//trim(tag)//' from file '//trim(veg_table_filename))
    end do

    ! physiological parameters
    call noahmp_ptable_read_tablestr(veg_table_filename, 'PHYS', tag, nlutyp, tblstr)
    call assert(nlutyp <= MLUTYP, 'Table PHYS#'//trim(tag)//' is too big to be read')
    do ii = 1, nlutyp
       read(tblstr(ii), *, iostat=ierr) index, &
            & LK_NROOT(ii), LK_CANWMXP(ii), LK_DLEAF(ii), LK_Z0MVT(ii), LK_HVT(ii), &
            & LK_HVB(ii), LK_DEN(ii), LK_RCROWN(ii), LK_CWPVT(ii)
       call assert(ierr == 0, 'Unable to read PHYS#'//trim(tag)//' from file '//trim(veg_table_filename))
    end do

    ! photosynthesis parameters
    call noahmp_ptable_read_tablestr(veg_table_filename, 'PHOTO', tag, nlutyp, tblstr)
    call assert(nlutyp <= MLUTYP, 'Table PHOTO#'//trim(tag)//' is too big to be read')
    do ii = 1, nlutyp
       read(tblstr(ii), *, iostat=ierr) index, &
            & LK_C3C4(ii), LK_RGL(ii), LK_HS(ii), LK_KC25(ii), LK_AKC(ii), &
            & LK_KO25(ii), LK_AKO(ii), LK_VCMX25(ii), LK_AVCMX(ii), LK_BP(ii), &
            & LK_RSMAX(ii), LK_RSMIN(ii), LK_MP(ii), LK_QE25(ii), LK_AQE(ii), &
            & LK_RMF25(ii), LK_RMS25(ii), LK_RMR25(ii), LK_FOLNMX(ii), LK_TOPT(ii), &
            & LK_TMIN(ii), LK_ARM(ii), LK_MRP(ii)
       call assert(ierr == 0, 'Unable to read PHOTO#'//trim(tag)//' from file '//trim(veg_table_filename))
    end do

    ! voc parameters
    call noahmp_ptable_read_tablestr(veg_table_filename, 'VOC', tag, nlutyp, tblstr)
    call assert(nlutyp <= MLUTYP, 'Table VOC#'//trim(tag)//' is too big to be read')
    do ii = 1, nlutyp
       read(tblstr(ii), *, iostat=ierr) index, LK_SLAREA(ii), (LK_EPS(jj,ii),jj=1,5)
       call assert(ierr == 0, 'Unable to read VOC#'//trim(tag)//' from file '//trim(veg_table_filename))
    end do
  end subroutine noahmp_veg_param_readptable
end module noahmp_veg_param