FUNCTION GET_TUNIV(ctx, z)

  !compute age of Universe in Gyr at redshift z
  !assumes flat universe w/ only matter and lambda
  !assumes om0,ol0,H0 set in sps_vars.f90
  
   USE fsps_context_types, ONLY: fsps_context_t
   USE fsps_constants, ONLY: SP
  IMPLICIT NONE
   TYPE(fsps_context_t), INTENT(IN) :: ctx
  INTEGER :: i
  INTEGER, PARAMETER :: ii=10000
  REAL(SP), INTENT(in) :: z
  REAL(SP) :: get_tuniv, thub
  REAL(SP), DIMENSION(ii) :: lnstig, hub

  !---------------------------------------------------------------!
  !---------------------------------------------------------------!

  get_tuniv = 0.0

  !Hubble time in Gyr
   thub = 0.978E3 / ctx%H0_val

  DO i=1,ii
     lnstig(i) = REAL(i)/ii*(LOG(1E4)-LOG(1+z))+LOG(1+z)
  ENDDO
  
   hub = SQRT( ctx%om0_val*EXP(lnstig)**3 + ctx%ol0_val )

  DO i=1,ii-1
     get_tuniv = get_tuniv + 0.5*(1/hub(i)+1/hub(i+1))
  ENDDO
  get_tuniv = get_tuniv * thub * (lnstig(2)-lnstig(1))


END FUNCTION GET_TUNIV
