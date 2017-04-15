module noahmp_const
  use iso_c_binding
  implicit none

  !%% Computer Constants
  integer, parameter :: i1 = C_INT8_T
  integer, parameter :: i2 = C_INT16_T
  integer, parameter :: i4 = C_INT32_T
  integer, parameter :: i8 = C_INT64_T
  integer, parameter :: r4 = C_FLOAT
  integer, parameter :: r8 = C_DOUBLE
  real(r8), parameter :: nan8 = transfer(-2251799813685248_i8, 1.0_r8)
  real(r4), parameter :: nan4 = transfer(-4194304_i4, 1.0_r4)
  real(r4), parameter :: MPE = 1.0E-6       !prevents division by zero erros

  !%% Physical Constants
  real(r4), parameter :: GRAV   = 9.80616   !acceleration due to gravity (m/s2)
  real(r4), parameter :: SB     = 5.67E-8   !Stefan-Boltzmann constant (w/m2/k4)
  real(r4), parameter :: RGAS   = 8.314     !univ. gas constant [J K-1 mol-1]
  real(r4), parameter :: KARMAN = 0.40      !von Karman constant
  real(r4), parameter :: TFRZ   = 273.15    !freezing/melting point (k)
  real(r4), parameter :: TTRI   = 273.16    !triple point (K)
  real(r4), parameter :: HSUB   = 2.8440E6  !latent heat of sublimation (j/kg)
  real(r4), parameter :: HVAP   = 2.5104E6  !latent heat of vaporization (j/kg)
  real(r4), parameter :: HFUS   = 0.3336E6  !latent heat of fusion (j/kg)
  real(r4), parameter :: CWAT   = 4.188E6   !specific heat capacity of water (j/m3/k)
  real(r4), parameter :: CICE   = 2.094E6   !specific heat capacity of ice (j/m3/k)
  real(r4), parameter :: CPAIR  = 1004.64   !heat capacity dry air at const pres (j/kg/k)
  real(r4), parameter :: TKWAT  = 0.6       !thermal conductivity of water (w/m/k)
  real(r4), parameter :: TKICE  = 2.2       !thermal conductivity of ice (w/m/k)
  real(r4), parameter :: TKAIR  = 0.023     !thermal conductivity of air (w/m/k)
  real(r4), parameter :: RAIR   = 287.04    !gas constant for dry air (j/kg/k)
  real(r4), parameter :: RVAP   = 461.269   !gas constant for  water vapor (j/kg/k)
  real(r4), parameter :: DENWAT = 1000.0    !density of water (kg/m3)
  real(r4), parameter :: DENICE = 917.0     !density of ice (kg/m3)

end module noahmp_const