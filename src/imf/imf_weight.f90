SUBROUTINE IMF_WEIGHT(ctx, mini, wght, nmass)

  !weight each star by the initial mass function (IMF)
  !such that the total initial population consists of 
  !one solar mass of stars.

  !This weighting scheme assumes that the luminosity, mass, etc.
  !does not vary within the mass bin.  The point is that we
  !want each element to represent the whole bin, from 
  !mass+/-0.5dm, rather than just the values at point i.
  !Then every intergral over mass is just a sum.

   USE fsps_context_types, ONLY: fsps_context_t
   USE sps_vars
  USE sps_utils, ONLY : imf, funcint
  IMPLICIT NONE

   TYPE(fsps_context_t), INTENT(INOUT) :: ctx
  REAL(SP), INTENT(inout), DIMENSION(nm) :: wght
  REAL(SP), INTENT(in), DIMENSION(nm)    :: mini
  INTEGER, INTENT(in) :: nmass
   INTEGER  :: i
   INTEGER :: imf_type_saved
   REAL(SP), DIMENSION(3) :: imf_alpha_saved
   REAL(SP) :: imf_vdmc_saved, imf_mdave_saved
  REAL(SP) :: m1,m2

  !--------------------------------------------------------!
  !--------------------------------------------------------!

  ASSOCIATE( &
     imf_lower_limit => ctx%state%imf_lower_limit, &
     imf_upper_limit => ctx%state%imf_upper_limit, &
     imf_lower_bound => ctx%state%imf_lower_bound )

  imf_type_saved = imf_type
  imf_alpha_saved = imf_alpha
  imf_vdmc_saved = imf_vdmc
  imf_mdave_saved = imf_mdave

  imf_type = ctx%imf_type_val
  imf_alpha = ctx%state%imf_alpha
  imf_vdmc = ctx%state%imf_vdmc
  imf_mdave = ctx%state%imf_mdave

  wght = 0.0

  DO i=1,nmass

     IF (mini(i).LT.imf_lower_limit.OR.&
          mini(i).GT.imf_upper_limit) CYCLE

     IF (i.EQ.1) THEN
        !note that this is not equal to imf_lower_limit
        !only for the Geneva models, which do not extend below 1.0 Msun
        m1 = imf_lower_bound
     ELSE
        m1 = mini(i) - 0.5*(mini(i)-mini(i-1))
     ENDIF
     IF (i.EQ.nmass) THEN
        m2 = mini(i)
     ELSE
        m2 = mini(i) + 0.5*(mini(i+1)-mini(i))
     ENDIF

     IF (m2.LT.m1) THEN
        WRITE(*,*) 'IMF_WEIGHT WARNING: non-monotonic mass!',m1,m2,m2-m1
        CYCLE
     ENDIF

     IF (m2.EQ.m1) CYCLE

     wght(i) = funcint(imf,m1,m2)

  ENDDO

   !normalize the weights as an integral from lower to upper limits
    imf_type = imf_type + 10
    wght = wght / funcint(imf,imf_lower_limit,imf_upper_limit)
    imf_type = imf_type - 10

   imf_type = imf_type_saved
   imf_alpha = imf_alpha_saved
   imf_vdmc = imf_vdmc_saved
   imf_mdave = imf_mdave_saved

  RETURN

   END ASSOCIATE

END SUBROUTINE IMF_WEIGHT

