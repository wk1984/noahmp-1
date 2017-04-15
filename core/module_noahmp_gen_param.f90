module noahmp_gen_param
  use noahmp_const, only: r4
  use noahmp_const, only: nan4
  use noahmp_global, only: nband
  implicit none

  ! maximum number of slope type
  integer, parameter :: MSLOPETYP = 30

  integer :: nslptyp = 0

  real(r4) :: LK_SLOPE(MSLOPETYP) = nan4 !slope index (0 - 1)
  real(r4) :: KK_CSOIL  = nan4  !vol. soil heat capacity [j/m3/K]
  real(r4) :: KK_ZBOT   = nan4  !Depth (m) of lower boundary soil temperature
  real(r4) :: KK_CZIL   = nan4  !Calculate roughness length of heat
                                ![constant in Zilithinkevich, S. S., 1995]

  real(r4) :: KK_DKREF  = nan4  !???
  real(r4) :: KK_KDTREF = nan4  !???
  real(r4) :: KK_FRZK   = nan4  !used in compute maximum infiltration rate (in INFIL)

  ! runoff parameters used for SIMTOP and SIMGM:
  real(r4) :: KK_TIMEAN  = nan4 !gridcell mean topgraphic index (global mean)
  real(r4) :: KK_FSATMAX = nan4 !maximum surface saturated fraction (global mean)

  ! adjustable parameters for snow processes
  ! melting factor (-) [Niu and Yang, 2007, JGR]
  real(r4) :: KK_MLTFCT = nan4
  ! snow surface roughness length (m)
  real(r4) :: KK_Z0SNO  = nan4
  !liquid water holding capacity for snowpack (m3 m-3)
  real(r4) :: KK_SSI    = nan4
  ! new snow mass to fully cover old snow (mm)
  ! (assuming snow density = 100 kg m-3)
  real(r4) :: KK_SWEMAX = nan4

  ! albedos of land ice, 1=vis, 2=nir
  real(r4) :: KK_ALBICE(nband) = nan4
  ! albedos of lake, 1=vis, 2=nir
  real(r4) :: KK_ALBLAKE(nband) = nan4
  ! two-stream parameters for snow
  real(r4) :: KK_OMEGAS(nband) = nan4
  real(r4) :: KK_BETADS = nan4
  real(r4) :: KK_BETAIS = nan4
  ! soil emissivity
  real(r4) :: KK_EMSSOIL = nan4
  ! lake emissivity
  real(r4) :: KK_EMSLAKE = nan4
contains

  subroutine noahmp_gen_param_readptable()
    use noahmp_utils, only: assert
    use noahmp_utils, only: noahmp_ptable_read_var
    use noahmp_utils, only: noahmp_ptable_read_tablestr
    implicit none
    integer, parameter :: PATHLEN = 512
    integer, parameter :: LINELEN = 512
    character(PATHLEN), parameter :: gen_table_filename = 'GENPARMMP.TBL'
    character(LINELEN) :: tblstr(MSLOPETYP)
    integer :: ii, index
    integer :: ierr

    call noahmp_ptable_read_tablestr(gen_table_filename, 'SLOPE', '', nslptyp, tblstr)
    call assert(nslptyp <= MSLOPETYP, 'Table SLOPE is too big to be read')
    do ii = 1, nslptyp
       read(tblstr(ii), *, iostat=ierr) index, LK_SLOPE(ii)
       call assert(ierr == 0, 'Unable to parse SLOPE in file '//trim(gen_table_filename))
    end do

    call noahmp_ptable_read_var(gen_table_filename, 'CSOIL', '', KK_CSOIL)
    call noahmp_ptable_read_var(gen_table_filename, 'DKREF', '', KK_DKREF)
    call noahmp_ptable_read_var(gen_table_filename, 'KDTREF', '', KK_KDTREF)
    call noahmp_ptable_read_var(gen_table_filename, 'FRZK', '', KK_FRZK)
    call noahmp_ptable_read_var(gen_table_filename, 'ZBOT', '', KK_ZBOT)
    call noahmp_ptable_read_var(gen_table_filename, 'CZIL', '', KK_CZIL)
    call noahmp_ptable_read_var(gen_table_filename, 'TIMEAN', '', KK_TIMEAN)
    call noahmp_ptable_read_var(gen_table_filename, 'FSATMAX', '', KK_FSATMAX)
    call noahmp_ptable_read_var(gen_table_filename, 'MLTFCT', '', KK_MLTFCT)
    call noahmp_ptable_read_var(gen_table_filename, 'Z0SNO', '', KK_Z0SNO)
    call noahmp_ptable_read_var(gen_table_filename, 'SSI', '', KK_SSI)
    call noahmp_ptable_read_var(gen_table_filename, 'SWEMAX', '', KK_SWEMAX)
    call noahmp_ptable_read_var(gen_table_filename, 'ALBICE', '', KK_ALBICE)
    call noahmp_ptable_read_var(gen_table_filename, 'ALBLAKE', '', KK_ALBLAKE)
    call noahmp_ptable_read_var(gen_table_filename, 'OMEGAS', '', KK_OMEGAS)
    call noahmp_ptable_read_var(gen_table_filename, 'BETADS', '', KK_BETADS)
    call noahmp_ptable_read_var(gen_table_filename, 'BETAIS', '', KK_BETAIS)
    call noahmp_ptable_read_var(gen_table_filename, 'EMSSOIL', '', KK_EMSSOIL)
    call noahmp_ptable_read_var(gen_table_filename, 'EMSLAKE', '', KK_EMSLAKE)
  end subroutine noahmp_gen_param_readptable
end module noahmp_gen_param