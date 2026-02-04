FUNCTION GET_LUMDIST(ctx, z)

  !compute luminosity distance to redshift z
  !assumes flat universe w/ only matter and lambda
  !assumes om0,ol0,H0 set in sps_vars.f90
  
  USE fsps_context_types, ONLY: fsps_context_t
  USE fsps_constants, ONLY: SP, C_LIGHT
  USE fsps_integration, ONLY: integrate_trapezoid_array
  IMPLICIT NONE
  TYPE(fsps_context_t), INTENT(IN) :: ctx
  INTEGER :: i
  INTEGER, PARAMETER :: ii=10000
  REAL(SP), INTENT(in) :: z
  REAL(SP) :: get_lumdist, dhub
  REAL(SP), DIMENSION(ii) :: zz, hub

  !---------------------------------------------------------------!
  !---------------------------------------------------------------!

  get_lumdist = 0.0

  !Hubble distance in pc
  dhub = C_LIGHT/1E13/ctx%H0_val*1E6

  DO i=1,ii
     zz(i) = REAL(i)/ii*z
  ENDDO
  
  hub = SQRT( ctx%om0_val*(1+zz)**3 + ctx%ol0_val )

  get_lumdist = integrate_trapezoid_array(zz,1/hub) * (1+z) * dhub


END FUNCTION GET_LUMDIST
