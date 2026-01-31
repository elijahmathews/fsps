SUBROUTINE GETMAGS(ctx, zred, spec, mags, mag_compute)

  !routine to calculate magnitudes in the Vega or AB systems,
  !given an input spectrum and redshift.
  !see parameter compute_vega_mags in sps_vars.f90
  !magnitudes defined in accordance with Fukugita et al. 1996, Eqn 7
  !This routine also redshifts the spectrum, if necessary.

   USE fsps_context_types, ONLY: fsps_context_t
   USE fsps_types, ONLY: SP, tiny_number, mag2cgs
   USE fsps_interpolation, ONLY: interpolate_linear
   USE fsps_integration, ONLY: integrate_trapezoid_array
   use, intrinsic :: ieee_arithmetic, only: ieee_is_nan
  IMPLICIT NONE

   TYPE(fsps_context_t), INTENT(INOUT) :: ctx

  INTEGER  :: i
  REAL(SP), INTENT(in) :: zred
   REAL(SP), INTENT(inout), DIMENSION(:) :: spec
   REAL(SP), INTENT(inout), DIMENSION(:) :: mags
   INTEGER, DIMENSION(:), INTENT(in), OPTIONAL  :: mag_compute
   INTEGER, DIMENSION(SIZE(mags)) :: magflag
   REAL(SP), DIMENSION(SIZE(spec))  :: tspec
  REAL(SP) :: const, dm
   INTEGER :: n_spec, n_bands

  !-----------------------------------------------------------!
  !-----------------------------------------------------------!

   n_spec = SIZE(spec)
   n_bands = SIZE(mags)

   IF (n_spec.LT.2) THEN
      mags = 99.0
      RETURN
   ENDIF

   ASSOCIATE( &
     spec_lambda => ctx%state%spec_lambda, bands => ctx%state%bands, &
     cosmospl => ctx%state%cosmospl, magvega => ctx%state%magvega, &
     compute_vega_mags => ctx%compute_vega_mags_val, &
     compute_light_ages => ctx%compute_light_ages_val )

  mags = 99.
  const= 0.0

  !set up the flags determining which mags are computed
  IF (PRESENT(mag_compute)) THEN
     magflag = mag_compute
     IF (compute_vega_mags.EQ.1) &
          magflag(1) = 1  !force V band to be computed
     IF (MAXVAL(magflag).EQ.0) & !no mags being computed so exit
          RETURN
  ELSE
     magflag = 1
  ENDIF

  !redshift the spectrum
  IF (ABS(zred).GT.tiny_number) THEN
     !write(*,*) "getmags: interpolating"
   DO i=1,n_spec
    tspec(i) = interpolate_linear(spec_lambda*(1+zred),spec, spec_lambda(i))
    IF (ieee_is_nan(tspec(i))) tspec(i) = 0.0
    tspec(i) = MAX(tspec(i),0.0)
   ENDDO

     !compute additional terms for cosmological mags
   dm    = interpolate_linear(cosmospl(:,1),cosmospl(:,3),zred)
     IF (ieee_is_nan(dm) .OR. dm.LE.tiny_number) THEN
        const = 0.0
     ELSE
        dm    = 5*LOG10(dm/10.)
        const = dm - 2.5*LOG10(1+zred)
     ENDIF

  ELSE

     tspec = spec

  ENDIF

  !integrate over each filter
   DO i=1,n_bands
     IF (magflag(i).EQ.0) CYCLE
    mags(i) = integrate_trapezoid_array(spec_lambda,tspec*bands(:,i)/spec_lambda)
       IF (ieee_is_nan(mags(i)) .OR. mags(i).LE.tiny_number) THEN
        mags(i) = 99.0
     ELSE
        IF (compute_light_ages.EQ.0) THEN
           !the mag2cgs var converts from Lsun/Hz to cgs at 10pc
           mags(i) = -2.5*LOG10(mags(i)) - 48.60 - 2.5*mag2cgs + const
        ENDIF
     ENDIF
  ENDDO

  !put magnitudes in the Vega system if keyword is set
  !(V-band is the first element in the array)
      IF (compute_vega_mags.EQ.1.AND.compute_light_ages.EQ.0) &
         mags(2:n_bands) = (mags(2:n_bands)-mags(1)) - &
         (magvega(2:n_bands)-magvega(1)) + mags(1)

   END ASSOCIATE


END SUBROUTINE GETMAGS
