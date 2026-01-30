SUBROUTINE ADD_REMNANTS(ctx, mass, maxmass)

  !add remnant (WD, NS, BH) masses back into the total mass
  !of the SSP.  These initial-mass-dependent remnant
  !formulae are taken from Renzini & Ciotti 1993.

   USE fsps_context_types, ONLY: fsps_context_t
   USE sps_vars
  USE sps_utils, ONLY : imf, funcint
  IMPLICIT NONE

   TYPE(fsps_context_t), INTENT(INOUT) :: ctx
  REAL(SP), INTENT(inout) :: mass
  !maximum mass still alive
  REAL(SP), INTENT(in) :: maxmass
   REAL(SP) :: minmass, imfnorm
   INTEGER :: imf_type_saved
   REAL(SP), DIMENSION(3) :: imf_alpha_saved
   REAL(SP) :: imf_vdmc_saved, imf_mdave_saved

  !---------------------------------------------------------------!
  !---------------------------------------------------------------!

  ASSOCIATE( &
     imf_lower_limit => ctx%state%imf_lower_limit, &
     imf_upper_limit => ctx%state%imf_upper_limit, &
     mlim_bh => ctx%state%mlim_bh, &
     mlim_ns => ctx%state%mlim_ns )

  imf_type_saved = imf_type
  imf_alpha_saved = imf_alpha
  imf_vdmc_saved = imf_vdmc
  imf_mdave_saved = imf_mdave

  imf_type = ctx%imf_type_val
  imf_alpha = ctx%state%imf_alpha
  imf_vdmc = ctx%state%imf_vdmc
  imf_mdave = ctx%state%imf_mdave

  !normalize the weights
   imf_type = imf_type+10
   imfnorm  = funcint(imf,imf_lower_limit,imf_upper_limit)
   imf_type = imf_type-10

  !BH remnants if any
  !40<M_max<imf_up leave behind a 0.5*M BH
  !if imf_upper_limit < 40 or < maxmass, add no mass since no stars could make BH and mlo=mhi=imf_upper_limit
  minmass = MIN(MAXVAL((/mlim_bh,maxmass/)), imf_upper_limit)
  imf_type = imf_type+10
  mass = mass + 0.5*funcint(imf,minmass,imf_upper_limit)/imfnorm
  imf_type = imf_type-10

  !Add NS remnants
  !8.5<M_max<40 also eave behind 1.4 Msun NS
  !if imf_upper_limit < 8.5 , add no mass since no stars could make NS, and mlo=mhi=imf_upper_limit
  IF (maxmass.LE.mlim_bh) THEN
     minmass = MIN(MAXVAL((/mlim_ns,maxmass/)), imf_upper_limit)
     mass = mass + 1.4*funcint(imf,minmass,MIN(mlim_bh, imf_upper_limit))/imfnorm
  ENDIF

  !Add WD remnants
  !M_max<8.5 also leave behind 0.077*M+0.48 WD
  !if imf_upper_limit < 8.5, only add mass if maxmass < imf_upper_limit since otherwise no stars have evolved and mlo=mhi=imf_upper_limit
  IF (maxmass.LE.8.5) THEN
     minmass = MIN(maxmass, imf_upper_limit)
     mass = mass + 0.48*funcint(imf,minmass,MIN(mlim_ns, imf_upper_limit))/imfnorm
     imf_type = imf_type+10
     mass = mass + 0.077*funcint(imf,minmass,MIN(mlim_ns, imf_upper_limit))/imfnorm
     imf_type = imf_type-10

  ENDIF

   imf_type = imf_type_saved
   imf_alpha = imf_alpha_saved
   imf_vdmc = imf_vdmc_saved
   imf_mdave = imf_mdave_saved

    END ASSOCIATE

  RETURN

END SUBROUTINE ADD_REMNANTS
