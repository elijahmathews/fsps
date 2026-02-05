SUBROUTINE SBF(ctx, pset, outfile)

  !routine to read in an isochrone (logt,logl,z) and produce 
  !SBFs for each point.  SBF magnitudes are a light-weighted average
  !of the stellar luminosities over stellar mass

   USE fsps_context_types, ONLY: fsps_context_t
      USE fsps_precision, ONLY: WP
      USE fsps_constants, ONLY: NM, BHB_SBS_TIME
      USE fsps_types, ONLY: PARAMS
      USE sps_utils, ONLY : getmags
      USE fsps_spectral_library, ONLY: get_stellar_spectrum
      USE fsps_stellar_modifications, ONLY: apply_blue_stragglers, modify_giant_branch, &
         modify_horizontal_branch
   USE fsps_imf, ONLY: compute_imf_weights
  IMPLICIT NONE

     TYPE(fsps_context_t), INTENT(INOUT) :: ctx
  CHARACTER(100), INTENT(in) :: outfile
  TYPE(PARAMS), INTENT(in)   :: pset
  INTEGER       :: i,j
  CHARACTER(34) :: fmt
  REAL(WP)      :: zero=0.0,hb_wght
     REAL(WP), DIMENSION(NM)     :: wght
     REAL(WP), ALLOCATABLE :: tspec(:),tspec2(:),spec1(:),spec2(:)
     REAL(WP), ALLOCATABLE :: mags(:)
     REAL(WP), ALLOCATABLE :: mini(:,:),mact(:,:),logl(:,:),logt(:,:),logg(:,:),ffco(:,:),phase(:,:),lmdot(:,:)
     INTEGER, ALLOCATABLE  :: nmass(:)
   REAL(WP), ALLOCATABLE :: time(:)

  !-----------------------------------------------------------!
  !-----------------------------------------------------------!

  !set up the format 
  ASSOCIATE( &
       nbands => ctx%state%nbands, nspec => ctx%state%nspec, nt => ctx%state%nt, &
       OUTPUT_HOME => ctx%output_home, &
       mini_isoc => ctx%state%mini_isoc, mact_isoc => ctx%state%mact_isoc, &
       logl_isoc => ctx%state%logl_isoc, logt_isoc => ctx%state%logt_isoc, &
       logg_isoc => ctx%state%logg_isoc, ffco_isoc => ctx%state%ffco_isoc, &
       lmdot_isoc => ctx%state%lmdot_isoc, phase_isoc => ctx%state%phase_isoc, &
       nmass_isoc => ctx%state%nmass_isoc, timestep_isoc => ctx%state%timestep_isoc )

  ALLOCATE(tspec(nspec),tspec2(nspec),spec1(nspec),spec2(nspec))
  ALLOCATE(mags(nbands))
  ALLOCATE(mini(nt,NM),mact(nt,NM),logl(nt,NM),logt(nt,NM),logg(nt,NM))
  ALLOCATE(ffco(nt,NM),phase(nt,NM),lmdot(nt,NM))
  ALLOCATE(nmass(nt),time(nt))

  fmt = '(F7.4,1x,3(F8.4,1x),000(F7.3,1x))'
  WRITE(fmt(21:23),'(I3,1x,I4)') nbands

  !reset arrays
  hb_wght  = 0.
  wght     = 0.
 
     OPEN(56,FILE=TRIM(OUTPUT_HOME)//'/OUTPUTS/'//TRIM(outfile)//&
       '.mags',STATUS='REPLACE')
  DO i=1,8 !write a dummy header
     WRITE(56,*) '#'
  ENDDO

  !transfer isochrones into temporary arrays
  mini  = mini_isoc(pset%zmet,:,:)  !initial mass
  mact  = mact_isoc(pset%zmet,:,:)  !actual (present) mass
  logl  = logl_isoc(pset%zmet,:,:)  !log(Lbol)
  logt  = logt_isoc(pset%zmet,:,:)  !log(Teff)
  logg  = logg_isoc(pset%zmet,:,:)  !log(g)
  ffco  = ffco_isoc(pset%zmet,:,:)  !is the TP-AGB star C-rich or O-rich?
  lmdot = lmdot_isoc(pset%zmet,:,:) !log Mdot
  phase = phase_isoc(pset%zmet,:,:) !flag indicating phase of evolution
  nmass = nmass_isoc(pset%zmet,:)   !number of elements per isochrone
  time  = timestep_isoc(pset%zmet,:)!age of each isochrone in log(yr)

  DO i=1,nt

     !compute IMF-based weights
   CALL compute_imf_weights(ctx, mini(i,:), wght, nmass(i))
     
     !modify the horizontal branch
     !need the hb weight for the blue stragglers too
   IF (pset%fbhb.GT.0.0.OR.pset%sbss.GT.1E-3) &
      CALL modify_horizontal_branch(ctx, i, pset%fbhb, time(i), hb_wght, nmass, &
      mini, mact, logl, logt, logg, phase, wght)

     !add in blue stragglers
   IF (time(i).GE.BHB_SBS_TIME.AND.pset%sbss.GT.1E-3) &
      CALL apply_blue_stragglers(ctx, i, pset%zmet, pset%sbss, hb_wght, nmass, &
      mini, mact, logl, logt, logg, phase, wght)

     !modify the TP-AGB stars and Post-AGB stars
   CALL modify_giant_branch(ctx, i, pset%zmet, time(i), nmass(i), pset%delt, pset%dell, pset%pagb, &
      pset%redgb, pset%agb, logl, logt, phase, wght)
 
     spec1 = 0.0
     spec2 = 0.0
     DO j=1,nmass(i)
           
        !get spectrum of ith star
         CALL get_stellar_spectrum(ctx, pset, mact(i,j), logt(i,j), 10**logl(i,j), logg(i,j), &
            phase(i,j), ffco(i,j), lmdot(i,j), tspec)

        !compute first and second moments of flux for
        !all stars and also by evolutionary phase
        spec2 = spec2 + wght(j)*tspec**2
        spec1 = spec1 + wght(j)*tspec

     ENDDO

     !compute the SBF
     tspec2 = spec2/spec1

   CALL GETMAGS(ctx, zero, tspec2, mags)
     WRITE(56,fmt) time(i),0.0,0.0,0.0,mags

  ENDDO

  CLOSE(56)

     END ASSOCIATE

END SUBROUTINE SBF
 
