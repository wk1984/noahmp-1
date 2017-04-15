module noahmp_soil_param
  use noahmp_const, only: r4
  use noahmp_const, only: nan4
  use noahmp_global, only: nband
  implicit none

  integer, parameter :: MSLTYP = 30 !max. number of soil types
  integer, parameter :: MSLCOL = 20 !max. number of soil colors

  character(256) :: sltyp_scheme = ' '
  integer :: nsltyp = 0   !number of soil types

  real(r4) :: LK_BEXP(MSLTYP) = nan4   !B parameter
  real(r4) :: LK_SMCMAX(MSLTYP) = nan4 !soil porosity (m3 m-3)
  real(r4) :: LK_SMCREF(MSLTYP) = nan4 !field capacity (m3 m-3)
  real(r4) :: LK_SMCWLT(MSLTYP) = nan4 !wilting point (m3 m-3)
  real(r4) :: LK_PSISAT(MSLTYP) = nan4 !saturation soil matric potential (m)
  real(r4) :: LK_DKSAT(MSLTYP)  = nan4 !saturation soil hydraulic conductivity
  real(r4) :: LK_DWSAT(MSLTYP)  = nan4 !???
  real(r4) :: LK_QUARTZ(MSLTYP) = nan4 !???

  real(r4) :: LK_KDT(MSLTYP) = nan4    !???
  real(r4) :: LK_FRZX(MSLTYP) = nan4   !???

  integer :: nsoilcol = 0

  real(r4) :: LK_ALBSAT(nband, MSLCOL) = nan4 !saturated soil albedos: 1=vis, 2=nir
  real(r4) :: LK_ALBDRY(nband, MSLCOL) = nan4 !dry soil albedos: 1=vis, 2=nir

contains
  subroutine noahmp_soil_param_readptable(tag)
    use noahmp_gen_param, only: KK_KDTREF
    use noahmp_gen_param, only: KK_DKREF
    use noahmp_gen_param, only: KK_FRZK
    use noahmp_utils, only: assert
    use noahmp_utils, only: noahmp_ptable_read_tablestr
    implicit none
    integer, parameter :: PATHLEN = 512
    integer, parameter :: LINELEN = 512
    character(512), parameter :: soil_table_filename = 'SOILPARMMP.TBL'
    character(*), intent(in) :: tag
    character(LINELEN) :: tblstr(max(MSLTYP,MSLCOL))
    integer :: ii, jj, index
    integer :: ierr

    ! soil thermal and hydraulic parameters
    call noahmp_ptable_read_tablestr(soil_table_filename, 'PARM', tag, nsltyp, tblstr)
    call assert(nsltyp <= MSLTYP, 'Table PARM.'//trim(tag)//' is too big to be read')
    do ii = 1, nsltyp
       read(tblstr(ii), *, iostat=ierr) index, LK_BEXP(ii), &
            & LK_SMCMAX(ii), LK_SMCREF(ii), LK_SMCWLT(ii), &
            & LK_PSISAT(ii), LK_DKSAT(ii), LK_DWSAT(ii), &
            & LK_QUARTZ(ii)
       call assert(ierr == 0, &
            & 'Unable to read PARM.'//trim(tag)//' from file '//trim(soil_table_filename))
    end do

    LK_KDT(:) = KK_KDTREF * LK_DKSAT(:) / KK_DKREF
    where (LK_SMCREF > 0.0)
       LK_FRZX = KK_FRZK * (LK_SMCMAX / LK_SMCREF) * (0.412 / 0468)
    end where

    ! soil radiation parameters
    call noahmp_ptable_read_tablestr(soil_table_filename, 'COLOR', '', nsoilcol, tblstr)
    call assert(nsoilcol <= MSLCOL, 'Table COLOR is too big to be read')
    do ii = 1, nsoilcol
       read(tblstr(ii), *, iostat=ierr) index, &
            & (LK_ALBSAT(jj,ii), jj=1,nband), &
            & (LK_ALBDRY(jj,ii), jj=1,nband)
       call assert(ierr == 0, 'Unable to read COLOR from file '//trim(soil_table_filename))
    end do
  end subroutine noahmp_soil_param_readptable
end module noahmp_soil_param