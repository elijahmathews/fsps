SUBROUTINE ADD_REMNANTS(ctx, mass, maxmass)

  !add remnant (WD, NS, BH) masses back into the total mass
  !of the SSP.  These initial-mass-dependent remnant
  !formulae are taken from Renzini & Ciotti 1993.

   USE fsps_context_types, ONLY: fsps_context_t
  USE fsps_types, ONLY: SP
  USE fsps_imf, ONLY: get_imf_value
  USE fsps_integration, ONLY: integrate_romberg
  IMPLICIT NONE

   TYPE(fsps_context_t), INTENT(INOUT) :: ctx
  REAL(SP), INTENT(inout) :: mass
  !maximum mass still alive
  REAL(SP), INTENT(in) :: maxmass
  REAL(SP) :: minmass, imfnorm

  !---------------------------------------------------------------!
  !---------------------------------------------------------------!

  ASSOCIATE( &
     imf_lower_limit => ctx%state%imf_lower_limit, &
     imf_upper_limit => ctx%state%imf_upper_limit, &
     mlim_bh => ctx%state%mlim_bh, &
     mlim_ns => ctx%state%mlim_ns )

  !normalize the weights (mass-weighted IMF integral)
   imfnorm  = integrate_romberg(ctx, imf_mass, imf_lower_limit, imf_upper_limit)

  !BH remnants if any
  !40<M_max<imf_up leave behind a 0.5*M BH
  !if imf_upper_limit < 40 or < maxmass, add no mass since no stars could make BH and mlo=mhi=imf_upper_limit
  minmass = MIN(MAXVAL((/mlim_bh,maxmass/)), imf_upper_limit)
  mass = mass + 0.5*integrate_romberg(ctx, imf_mass, minmass, imf_upper_limit)/imfnorm

  !Add NS remnants
  !8.5<M_max<40 also eave behind 1.4 Msun NS
  !if imf_upper_limit < 8.5 , add no mass since no stars could make NS, and mlo=mhi=imf_upper_limit
  IF (maxmass.LE.mlim_bh) THEN
     minmass = MIN(MAXVAL((/mlim_ns,maxmass/)), imf_upper_limit)
  mass = mass + 1.4*integrate_romberg(ctx, imf_number, minmass, MIN(mlim_bh, imf_upper_limit))/imfnorm
  ENDIF

  !Add WD remnants
  !M_max<8.5 also leave behind 0.077*M+0.48 WD
  !if imf_upper_limit < 8.5, only add mass if maxmass < imf_upper_limit since otherwise no stars have evolved and mlo=mhi=imf_upper_limit
  IF (maxmass.LE.8.5) THEN
     minmass = MIN(maxmass, imf_upper_limit)
  mass = mass + 0.48*integrate_romberg(ctx, imf_number, minmass, MIN(mlim_ns, imf_upper_limit))/imfnorm
  mass = mass + 0.077*integrate_romberg(ctx, imf_mass, minmass, MIN(mlim_ns, imf_upper_limit))/imfnorm

  ENDIF

    END ASSOCIATE

   RETURN

CONTAINS

  pure function imf_number(ctx, x) result(res)
    type(fsps_context_t), intent(in) :: ctx
    real(SP), dimension(:), intent(in) :: x
    real(SP), dimension(size(x)) :: res
    res = get_imf_value(ctx, x, mass_weighted=.false.)
  end function imf_number

  pure function imf_mass(ctx, x) result(res)
    type(fsps_context_t), intent(in) :: ctx
    real(SP), dimension(:), intent(in) :: x
    real(SP), dimension(size(x)) :: res
    res = get_imf_value(ctx, x, mass_weighted=.true.)
  end function imf_mass

END SUBROUTINE ADD_REMNANTS
