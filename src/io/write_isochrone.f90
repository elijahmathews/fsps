SUBROUTINE WRITE_ISOCHRONE(ctx, outfile, pset)

  !routine to write all isochrones and CMDs at a given metallicity
  !note that the output age grid is the native spacing, not boosted
  !by the parameter time_res_incr

     USE fsps_context_types, ONLY: fsps_context_t
     USE fsps_types, ONLY: SP, PARAMS, nm, bhb_sbs_time, gsig4pi
     USE sps_utils, ONLY : getmags,getspec
     USE fsps_stellar_modifications, ONLY: apply_blue_stragglers, modify_giant_branch, &
          modify_horizontal_branch
     USE fsps_imf, ONLY: compute_imf_weights
  IMPLICIT NONE

     TYPE(fsps_context_t), INTENT(INOUT) :: ctx
  INTEGER :: i,tt,zz
  TYPE(PARAMS), INTENT(in) :: pset
  CHARACTER(100), INTENT(in)  :: outfile
  CHARACTER(60)  :: fmt
  REAL(SP) :: dz=0.0,loggi,hb_wght
  REAL(SP), DIMENSION(nm)     :: wght
  REAL(SP), ALLOCATABLE :: spec(:)
  REAL(SP), ALLOCATABLE :: mags(:)
  !temp arrays for the isochrone data
  REAL(SP), ALLOCATABLE :: mini(:,:),mact(:,:),logl(:,:),logt(:,:),logg(:,:),&
       ffco(:,:),phase(:,:),lmdot(:,:)
  INTEGER, ALLOCATABLE :: nmass(:)

  !---------------------------------------------------------------!
  !---------------------------------------------------------------!

  ASSOCIATE( &
       nbands => ctx%state%nbands, nspec => ctx%state%nspec, nt => ctx%state%nt, &
       OUTPUT_HOME => ctx%output_home, &
       isoc_type => ctx%state%isoc_type, &
       mini_isoc => ctx%state%mini_isoc, mact_isoc => ctx%state%mact_isoc, &
       logl_isoc => ctx%state%logl_isoc, logt_isoc => ctx%state%logt_isoc, &
       logg_isoc => ctx%state%logg_isoc, ffco_isoc => ctx%state%ffco_isoc, &
       lmdot_isoc => ctx%state%lmdot_isoc, phase_isoc => ctx%state%phase_isoc, &
       nmass_isoc => ctx%state%nmass_isoc, timestep_isoc => ctx%state%timestep_isoc, &
       zlegend => ctx%state%zlegend, mact_isoc_full => ctx%state%mact_isoc )

  ALLOCATE(spec(nspec))
  ALLOCATE(mags(nbands))
  ALLOCATE(mini(nt,nm),mact(nt,nm),logl(nt,nm),logt(nt,nm),logg(nt,nm))
  ALLOCATE(ffco(nt,nm),phase(nt,nm),lmdot(nt,nm))
  ALLOCATE(nmass(nt))

  hb_wght = 0.0
  wght    = 0.0
  zz      = pset%zmet

  fmt = '(F7.4,1x,F8.4,1x,F14.9,1x,F14.9,1x,7(F8.4,1x),000(F7.3,1x))'
  WRITE(fmt(47:49),'(I3,1x,I4)') nbands

     OPEN(40,FILE=TRIM(OUTPUT_HOME)//'/OUTPUTS/'//TRIM(outfile)//'.cmd',&
       STATUS='REPLACE')
  WRITE(40,*) '# age log(Z) mini mact logl logt logg '//&
       'phase composition log(weight) log(mdot) mags'
       
  !transfer isochrones into temporary arrays
  mini  = mini_isoc(zz,:,:)  !initial mass
  mact  = mact_isoc(zz,:,:)  !actual (present) mass
  logl  = logl_isoc(zz,:,:)  !log(Lbol)
  logt  = logt_isoc(zz,:,:)  !log(Teff)
  logg  = logg_isoc(zz,:,:)  !log(g)
  ffco  = ffco_isoc(zz,:,:)  !is the TP-AGB star C-rich or O-rich?
  lmdot = lmdot_isoc(zz,:,:) !log Mdot
  phase = phase_isoc(zz,:,:) !flag indicating phase of evolution
  nmass = nmass_isoc(zz,:)   !number of elements per isochrone

  DO tt=1,nt

     !compute IMF-based weights
     CALL compute_imf_weights(ctx, mini(tt,:), wght, nmass(tt))

     !modify the horizontal branch
     !need the hb weight for the blue stragglers too
     IF (pset%fbhb.GT.0.0.OR.pset%sbss.GT.1E-3) &
          CALL modify_horizontal_branch(ctx, tt, pset%fbhb, timestep_isoc(zz,tt), hb_wght, nmass, &
          mini, mact, logl, logt, logg, phase, wght)

     !add in blue stragglers
     IF (timestep_isoc(zz,tt).GE.bhb_sbs_time.AND.pset%sbss.GT.1E-3) &
          CALL apply_blue_stragglers(ctx, tt, pset%sbss, hb_wght, nmass, &
          mini, mact, logl, logt, logg, phase, wght)

     !modify the RGB and/or AGB stars
     CALL modify_giant_branch(ctx, tt, zz, timestep_isoc(zz,tt), nmass(tt), pset%delt, &
          pset%dell, pset%pagb, pset%redgb, pset%agb, logl, logt, phase, wght)

     DO i=1,nmass(tt)
        
        !get the spectrum
        CALL GETSPEC(ctx, pset, mact(tt,i), logt(tt,i), 10**logl(tt,i), &
             logg(tt,i), phase(tt,i), ffco(tt,i), lmdot(tt,i), wght(i), spec)
        !calculate magnitudes
     CALL GETMAGS(ctx, dz, spec, mags)

        IF (isoc_type.EQ.'bsti') THEN
           loggi = LOG10( gsig4pi*mact_isoc_full(zz,tt,i)/&
                logl(tt,i) ) + 4*logt(tt,i)
        ELSE
           loggi = logg(tt,i)
        ENDIF

        !write results to file
        WRITE(40,fmt) timestep_isoc(zz,tt),LOG10(zlegend(zz)),&
             mini(tt,i),mact(tt,i),logl(tt,i),logt(tt,i),loggi,phase(tt,i),&
             ffco(tt,i),LOG10(wght(i)),lmdot(tt,i),mags
        
     ENDDO

  ENDDO

  CLOSE(40)

     END ASSOCIATE

END SUBROUTINE WRITE_ISOCHRONE
