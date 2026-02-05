!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  Routine to calculate the evolution of a single stellar          !
!  population from a set of input theoretical isochrones and a     !
!  heterogeneous library of stellar spectra.  The code also allows !
!  for variation in the horizontal branch morphology, TP-AGB       !
!  phase, and the blue straggler population.  The output is a      !
!  time-dependent spectrum from the far-UV to the far-IR.          !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   PARAMETER RANGES:
!   1.  if (fbhb,sbs)<1E-3 then (fbhb,sbs)=0.0
!   2.  if abs(delt)>0.5 then abs(delt)=0.5
!
!-----------------------------------------------------------!
!-----------------------------------------------------------!

SUBROUTINE SSP_GEN(ctx, pset, mass_ssp, lbol_ssp, spec_ssp)

   USE fsps_precision, ONLY: WP
   USE fsps_constants, ONLY: NM, VERBOSE, BHB_SBS_TIME, TIME_RES_INCR
   USE fsps_types, ONLY: PARAMS
   USE sps_utils, ONLY: getspec
   USE fsps_smoothing, ONLY: apply_smoothing
  USE fsps_stellar_modifications, ONLY: apply_blue_stragglers, modify_giant_branch, &
     modify_horizontal_branch, add_remnant_mass, add_xray_binaries
  USE fsps_gas, only: apply_nebular_emission
  USE fsps_interpolation, ONLY: find_interval
   USE fsps_imf, ONLY: compute_imf_weights
  USE fsps_context_types, ONLY: fsps_context_t
  IMPLICIT NONE

  TYPE(fsps_context_t), INTENT(INOUT) :: ctx
  INTEGER :: i=1, j=1, stat,ii,klo,khi !,tlo,thi
  !weight given to the entire horizontal branch
     REAL(WP) :: hb_wght,dt,tco
  !array of IMF weights
     REAL(WP), DIMENSION(NM) :: wght
  !SSP spectrum
     REAL(WP), INTENT(inout), DIMENSION(:,:) :: spec_ssp
     REAL(WP), DIMENSION(SIZE(spec_ssp,1), SIZE(spec_ssp,2)) :: tspec_ssp
  !Mass and Lbol info
     REAL(WP), INTENT(inout), DIMENSION(:) :: mass_ssp, lbol_ssp

  !temp arrays for the isochrone data
     REAL(WP), ALLOCATABLE :: mini(:,:),mact(:,:),logl(:,:),logt(:,:),logg(:,:),&
     ffco(:,:),phase(:,:),lmdot(:,:)
   REAL(WP), DIMENSION(NM) :: temp_mini
  !arrays holding the number of mass elements for each
  !isochrone and the age of each isochrone
  INTEGER, ALLOCATABLE :: nmass(:)
   REAL(WP), ALLOCATABLE :: time(:)
   REAL(WP), DIMENSION(SIZE(spec_ssp,1)) :: tspec
  !structure containing all necessary parameters
  !(TYPE objects defined in sps_vars.f90)
  TYPE(PARAMS), INTENT(in) :: pset
  !CHARACTER(2) :: istr,istr2

  !-----------------------------------------------------------!
  !--------------------------Setup----------------------------!
  !-----------------------------------------------------------!

  ASSOCIATE( &
       check_sps_setup => ctx%state%check_sps_setup, &
       nz => ctx%state%nz, nt => ctx%state%nt, ntfull => ctx%state%ntfull, &
       nspec => ctx%state%nspec, isoc_type => ctx%state%isoc_type, &
       bpass_spec_ssp => ctx%state%bpass_spec_ssp, bpass_mass_ssp => ctx%state%bpass_mass_ssp, &
       imf_type => ctx%imf_type_val, imf_alpha => ctx%state%imf_alpha, &
       imf_vdmc => ctx%state%imf_vdmc, imf_mdave => ctx%state%imf_mdave, &
     n_user_imf => ctx%state%n_user_imf, imf_user_alpha => ctx%state%imf_user_alpha, &
       imf_lower_limit => ctx%state%imf_lower_limit, imf_upper_limit => ctx%state%imf_upper_limit, &
       mini_isoc => ctx%state%mini_isoc, mact_isoc => ctx%state%mact_isoc, &
       logl_isoc => ctx%state%logl_isoc, logt_isoc => ctx%state%logt_isoc, &
       logg_isoc => ctx%state%logg_isoc, ffco_isoc => ctx%state%ffco_isoc, &
       phase_isoc => ctx%state%phase_isoc, lmdot_isoc => ctx%state%lmdot_isoc, &
       nmass_isoc => ctx%state%nmass_isoc, timestep_isoc => ctx%state%timestep_isoc, &
       zlegend => ctx%state%zlegend, time_full => ctx%state%time_full, &
       lsfinfo => ctx%state%lsfinfo, spec_lambda => ctx%state%spec_lambda, &
       add_stellar_remnants => ctx%add_stellar_remnants_val, &
       add_neb_emission => ctx%add_neb_emission_val, add_xrb_emission => ctx%add_xrb_emission_val, &
       smooth_lsf => ctx%smooth_lsf_val, smooth_velocity => ctx%smooth_velocity_val, &
       sps_home => ctx%sps_home )

   ALLOCATE(mini(nt,NM),mact(nt,NM),logl(nt,NM),logt(nt,NM),logg(nt,NM))
   ALLOCATE(ffco(nt,NM),phase(nt,NM),lmdot(nt,NM))
   ALLOCATE(nmass(nt),time(nt))

  IF (check_sps_setup.EQ.0) THEN
     WRITE(*,*) 'SSP_GEN ERROR0: '//&
          'SPS_SETUP must be run once before calling SSP_GEN. '
     STOP
  ENDIF

  !reset arrays
  spec_ssp = 0.
  mass_ssp = 0.
  lbol_ssp = 0.

  !test metallicity range
  IF (pset%zmet.LT.1.OR.pset%zmet.GT.nz) THEN
     WRITE(*,*) 'SSP_GEN ERROR: metallicity outside of range',pset%zmet
     STOP
  ENDIF

  IF (isoc_type.EQ.'bpss') THEN

     !the BPASS SSPs are stored in a master array so simply
     !pull the relevant metallicity model into spec_ssp
     spec_ssp = bpass_spec_ssp(:,:,pset%zmet)
     mass_ssp = bpass_mass_ssp(:,pset%zmet)

  ELSE

     IF (imf_type.NE.0.AND.imf_type.NE.1.AND.imf_type.NE.2.&
          .AND.imf_type.NE.3.AND.imf_type.NE.4.AND.imf_type.NE.5) THEN
        WRITE(*,*) 'SSP_GEN ERROR: IMF type outside of range',imf_type
        STOP
     ENDIF

     !dump IMF parameters into common block
     imf_alpha(1) = pset%imf1
     imf_alpha(2) = pset%imf2
     imf_alpha(3) = pset%imf3
     imf_vdmc     = pset%vdmc
     imf_mdave    = pset%mdave

     !read in user-defined IMF (this needs to be done here rather than
     !in sps_setup because the user can change the IMF without having
     !to re-run the setup
     IF (imf_type.EQ.5) THEN
        IF (TRIM(pset%imf_filename).EQ.'') THEN
           OPEN(13,FILE=TRIM(sps_home)//'/data/imf.dat',ACTION='READ',STATUS='OLD')
        ELSE
           OPEN(13,FILE=TRIM(sps_home)//'/data/'//TRIM(pset%imf_filename),&
                ACTION='READ',STATUS='OLD')
        ENDIF
        DO i=1,100
           READ(13,*,IOSTAT=stat) imf_user_alpha(1,i),imf_user_alpha(2,i),&
                imf_user_alpha(3,i)
           IF (stat.NE.0) GOTO 29
        ENDDO
        WRITE(*,*) 'SSP_GEN ERROR: didnt finish reading in the imf file'
        STOP
29      CONTINUE
        CLOSE(13)
        n_user_imf = i-1
        !define the upper and lower IMF limits
        imf_lower_limit = imf_user_alpha(1,1)
        imf_upper_limit = imf_user_alpha(2,n_user_imf)
     ENDIF

     !transfer isochrones into temporary arrays
     mini  = mini_isoc(pset%zmet,:,:)  !initial mass
     mact  = mact_isoc(pset%zmet,:,:)  !actual (present) mass
     logl  = logl_isoc(pset%zmet,:,:)  !log(Lbol)
     logt  = logt_isoc(pset%zmet,:,:)  !log(Teff)
     logg  = logg_isoc(pset%zmet,:,:)  !log(g)
     ffco  = ffco_isoc(pset%zmet,:,:)  !is the TP-AGB star C-rich or O-rich?
     phase = phase_isoc(pset%zmet,:,:) !flag indicating phase of evolution
     lmdot = lmdot_isoc(pset%zmet,:,:) !log Mdot
     nmass = nmass_isoc(pset%zmet,:)   !number of elements per isochrone
     time  = timestep_isoc(pset%zmet,:)!age of each isochrone in log(yr)

     !write for control
     IF (VERBOSE.EQ.1) THEN
        WRITE(*,*)
        WRITE(*,'("   Log(Z/Zsol): ",F6.3)') LOG10(zlegend(pset%zmet)/0.019)
        WRITE(*,'("   Fraction of blue HB stars: ",F6.3)') pset%fbhb
        WRITE(*,'("   Ratio of BS to HB stars  : ",F6.3)') pset%sbss
        WRITE(*,'("   Shift to TP-AGB [log(Teff),log(Lbol)]: ",F5.2,1x,F5.2)') &
             pset%delt, pset%dell
        IF (imf_type.EQ.2) THEN
           WRITE(*,'("   IMF: ",I1,", slopes= ",3F4.1)') &
                imf_type,imf_alpha
        ELSE IF (imf_type.EQ.3) THEN
           WRITE(*,'("   IMF: ",I1,", cut-off= ",F4.2)') imf_type,imf_vdmc
        ELSE
           WRITE(*,'("   IMF: ",I1)') imf_type
        ENDIF
     ENDIF

     !-----------------------------------------------------------!
     !---------------------Generate SSPs-------------------------!
     !-----------------------------------------------------------!

     !loop over each isochrone
     DO i=1,nt

        !flag that allows us to compute only a subset of models
        IF (pset%ssp_gen_age(i).EQ.0) CYCLE

        !reset arrays
        hb_wght  = 0.
        wght     = 0.

        IF (VERBOSE.EQ.1) &
             WRITE(*,'("age=",F5.2)') time(i)

      ! Manually copy the row to a contiguous temp array
      temp_mini(1:nmass(i)) = mini(i, 1:nmass(i))

      !compute IMF-based weights
      CALL compute_imf_weights(ctx, temp_mini, wght, nmass(i))
        !modify the horizontal branch
        !need the hb weight for the blue stragglers too
         IF (pset%fbhb.GT.0.0.OR.pset%sbss.GT.1E-3) &
            CALL modify_horizontal_branch(ctx, i, pset%fbhb, time(i), hb_wght, nmass, &
            mini, mact, logl, logt, logg, phase, wght)
        !add in blue stragglers
         IF (time(i).GE.BHB_SBS_TIME.AND.pset%sbss.GT.1E-3) &
            CALL apply_blue_stragglers(ctx, i, pset%sbss, hb_wght, nmass, &
            mini, mact, logl, logt, logg, phase, wght)
        !modify the TP-AGB stars and Post-AGB stars
         CALL modify_giant_branch(ctx, i, pset%zmet, time(i), nmass(i), pset%delt, pset%dell, &
            pset%pagb, pset%redgb, pset%agb, logl, logt, phase, wght)
        ii = 1 + (i-1)*TIME_RES_INCR

        !compute IMF-weighted mass of the SSP
        mass_ssp(ii) = SUM(wght(1:nmass(i))*mact(i,1:nmass(i)))

        !add in remant masses
        IF (add_stellar_remnants.EQ.1) THEN
           CALL add_remnant_mass(ctx, mass_ssp(ii), MAXVAL(mini(i,:)))
        ENDIF

        !compute IMF-weighted bolometric luminosity (actually log(Lbol))
        lbol_ssp(ii) = LOG10(SUM(wght(1:nmass(i))*10**logl(i,1:nmass(i))))

        !compute SSP spectrum
        spec_ssp(:,ii) = 0.
        DO j=1,nmass(i)

           tco = ffco(i,j)
           IF (phase(i,j).EQ.5.AND.tco.GT.1.0) THEN
              !dilute the C star fraction
              !IF (1.0.GE.pset%fcstar) tco = 1.0
           ENDIF

          CALL GETSPEC(ctx, pset, mact(i,j), logt(i,j), &
             10**logl(i,j),logg(i,j),phase(i,j),tco,lmdot(i,j),&
             wght(j)/MAXVAL(wght(1:nmass(i))*10**logl(i,1:nmass(i))),tspec)

           !only construct SSPs for particular evolutionary
           !phases if evtype NE -1
           IF ((pset%evtype.EQ.-1.OR.pset%evtype.EQ.phase(i,j))&
                .AND.mini(i,j).LT.pset%masscut) &
                spec_ssp(:,ii) = wght(j)*tspec + spec_ssp(:,ii)

        ENDDO

     ENDDO

  ENDIF

  !-------------------------------------------------------------!
  !-now interpolate the SSPs to fill out the expanded time grid-!
  !-------------------------------------------------------------!

  IF (TIME_RES_INCR.GT.1) THEN
     DO j=1,ntfull
        IF (MOD(j-1,TIME_RES_INCR).EQ.0) CYCLE
      klo = MAX(MIN(find_interval(time,time_full(j)),nt-1),1)
        dt  = (time_full(j)-time(klo))/(time(klo+1)-time(klo))
        klo = 1+(klo-1)*TIME_RES_INCR
        khi = klo+TIME_RES_INCR
        spec_ssp(:,j) = 10**( (1-dt)*LOG10(spec_ssp(:,klo)) + &
             dt*LOG10(spec_ssp(:,khi)))
        lbol_ssp(j)   = (1-dt)*lbol_ssp(klo)   + dt*lbol_ssp(khi)
        mass_ssp(j)   = (1-dt)*mass_ssp(klo)   + dt*mass_ssp(khi)
     ENDDO
  ENDIF

  !-------------------------------------------------------------!
  !-------add the nebular emission model at the SSP level-------!
  !-------------------------------------------------------------!

  IF (add_neb_emission.EQ.2) THEN
   CALL apply_nebular_emission(ctx, pset, spec_ssp, tspec_ssp)
     spec_ssp = tspec_ssp
  ENDIF

  !-------------------------------------------------------------!
  !---------------add X-ray binaries the SSP level--------------!
  !-------------------------------------------------------------!

  IF (add_xrb_emission.EQ.1) THEN
   CALL add_xray_binaries(ctx, pset, spec_ssp, tspec_ssp)
     spec_ssp = tspec_ssp
  ENDIF

  !-------------------------------------------------------------!
  !--------now smooth by an instrumental LSF if provided--------!
  !-------------------------------------------------------------!

  IF (smooth_lsf.EQ.1) THEN
     DO j=1,ntfull
         CALL apply_smoothing(ctx, spec_lambda, spec_ssp(:,j), 99.d0, lsfinfo%minlam, &
            lsfinfo%maxlam, lsfinfo%lsf)
     ENDDO
  ENDIF

  END ASSOCIATE

END SUBROUTINE SSP_GEN


