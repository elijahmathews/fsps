SUBROUTINE ADD_DUST(ctx, pset, csp1, csp2, specdust, mdust, ncsp1, ncsp2, nebdust)

  ! Inputs
  ! ---------
  !
  ! pset:
  !   A `PARAMS` structure containing the dust parameters
  !
  !  csp1:
  !   The spectrum of the young stars
  !
  !  csp2:
  !   The spectrum of the old stars
  !
  !  ncsp1:
  !   The nebular fluxes for the young stars
  !
  !  ncsp2:
  !   The nebular fluxes for the old stars
  !
  ! Outputs
  ! ---------
  !
  ! specdust:
  !  The combined spectrum including dust attenuation and dust emission
  !
  ! nebdust:
  !   The combined nebular emission line fluxes including dust absorption
  !
  ! mdust:
  !  The dust mass required to produce the absorbed luminosity for the given dust emission parameters


   USE fsps_types, ONLY: SP, PARAMS, nemline, mypi, clight, tiny_number
   USE fsps_context_types, ONLY: fsps_context_t
   USE fsps_dust, ONLY: compute_attenuation_curve
   USE fsps_interpolation, ONLY: find_interval, interpolate_linear
   USE fsps_integration, ONLY: integrate_trapezoid_array
   use, intrinsic :: ieee_arithmetic, only: ieee_is_nan
  IMPLICIT NONE

   TYPE(fsps_context_t), INTENT(INOUT) :: ctx

   REAL(SP), DIMENSION(:), INTENT(in) :: csp1,csp2
  TYPE(PARAMS), INTENT(in) :: pset
   REAL(SP), DIMENSION(:), INTENT(out) :: specdust
  REAL(SP), INTENT(out) :: mdust
  REAL(SP), DIMENSION(nemline), INTENT(in) :: ncsp1,ncsp2
  REAL(SP), DIMENSION(nemline), INTENT(out) :: nebdust
   INTEGER :: qlo,ulo,iself=0
   REAL(SP), DIMENSION(SIZE(csp1))  :: diff_dust,tau_diff,cspi
  REAL(SP), DIMENSION(nemline)  :: diff_dust_neb,ncspi
   REAL(SP), DIMENSION(SIZE(csp1))  :: nu,dumin,dumax
    REAL(SP), DIMENSION(SIZE(csp1))  :: mduste,duste,oduste,tduste
   REAL(SP) :: lboln,lbold,labs,gamma,norm,dq,du

  !---------------------------------------------------------------!

  ASSOCIATE( &
     dust_type => ctx%dust_type_val, &
     add_dust_emission => ctx%add_dust_emission_val, &
     nebemlineinspec => ctx%nebemlineinspec_val, &
     nqpah_dustem => ctx%state%nqpah_dustem, &
     numin_dustem => ctx%state%numin_dustem, &
     ndim_dustem => ctx%state%ndim_dustem, &
     spec_lambda => ctx%state%spec_lambda, &
     nebem_line_pos => ctx%state%nebem_line_pos, &
     qpaharr => ctx%state%qpaharr, &
     uminarr => ctx%state%uminarr, &
     dustem_dustem => ctx%state%dustem_dustem, &
     dustem2_dustem => ctx%state%dustem2_dustem )
  !----------------------Test input params------------------------!
  !---------------------------------------------------------------!

  IF (dust_type.LT.0.OR.dust_type.GT.6) THEN
     WRITE(*,*) 'ADD_DUST ERROR: unknown dust_type: ', dust_type
     STOP
  ENDIF

  IF (pset%uvb.LT.0) THEN
     WRITE(*,*) 'ADD_DUST ERROR: pset%uvb<0!'
     STOP
  ENDIF

  IF (pset%wgp1.LT.1.OR.pset%wgp2.LT.1) THEN
     WRITE(*,*) 'ADD_DUST ERROR: pset%wgp1 and/or pset%wgp2 '//&
          'out of bounds: ',pset%wgp1,pset%wgp2
     STOP
  ENDIF

  IF (pset%dust_clumps.GT.tiny_number) THEN
     WRITE(*,*) 'ADD_DUST ERROR: dust_clumps feature no longer supported'
     STOP
  ENDIF

  IF (pset%frac_obrun.LT.0.0.OR.pset%frac_obrun.GT.1.0.OR.&
       pset%frac_nodust.LT.0.0.OR.pset%frac_nodust.GT.1.0) THEN
     WRITE(*,*) 'ADD_DUST ERROR: frac_obrun and/or frac_nodust out of bounds'
     STOP
  ENDIF

  !---------------------------------------------------------------!
  !----------------------Add dust absorption----------------------!
  !---------------------------------------------------------------!

  !compute attenuation curve for diffuse dust
   tau_diff = compute_attenuation_curve(spec_lambda, dust_type, pset, ctx)

  !combine old and young stars, attenuating the young
  !with a fixed power-law attn curve, and allowing a fraction
  !of the young stars to be dust-free ("OB runaways")
  !and an extra attenuation towards old stars (patchy dust)
  cspi = csp1 * EXP(-pset%dust1*(spec_lambda/5500.)**(pset%dust1_index))*&
       (1-pset%frac_obrun) + csp1*pset%frac_obrun + &
       csp2 * EXP(-pset%dust3*tau_diff)

  !normalize diffuse dust curve unless Witt & Gordon
  IF (dust_type.EQ.3) THEN
    diff_dust = EXP(-tau_diff)
  ELSE
    diff_dust = EXP(-pset%dust2 * tau_diff)
  ENDIF

  !allow a fraction of the diffuse dust spectrum to be dust-free
  specdust  = (1-pset%frac_nodust) * cspi*diff_dust + &
               cspi*pset%frac_nodust

  !as above, for nebular line luminosities
   diff_dust_neb = interpolate_linear(spec_lambda,diff_dust,nebem_line_pos)
   WHERE (ieee_is_nan(diff_dust_neb)) diff_dust_neb = 1.0
  ncspi = ncsp1 * EXP(-pset%dust1*(nebem_line_pos/5500.)**(pset%dust1_index))*&
       (1-pset%frac_obrun) + ncsp1*pset%frac_obrun + ncsp2
  nebdust  = ncspi*diff_dust_neb*(1-pset%frac_nodust) + &
       ncspi*pset%frac_nodust

  !---------------------------------------------------------------!
  !----------------------Add dust emission------------------------!
  !---------------------------------------------------------------!

  ! The dust spectrum is computed according to Draine & Li 2007
  ! we are computing the integral of the dust spectrum over P(U)dU
  ! by considering two components, a delta function at Umin and
  ! a power-law distribution from Umin to Umax=1E6 and alpha=2
  ! the relative weights of the two components are given by gamma

  IF (add_dust_emission.EQ.1) THEN

     !only add dust emission if there is absorption
     IF (pset%dust2.GT.tiny_number.OR.pset%dust1.GT.tiny_number) THEN

        iself=0
        !compute Lbol both before and after dust attenuation
        !this will determine the normalization of the dust emission
        nu    = clight/spec_lambda
         lbold = integrate_trapezoid_array(nu,specdust)
         lboln = integrate_trapezoid_array(nu,csp1+csp2)
            IF (ieee_is_nan(lbold)) lbold = 0.0
            IF (ieee_is_nan(lboln)) lboln = 0.0
        IF (nebemlineinspec.EQ.0) THEN
           lboln = lboln + SUM(ncsp1) + SUM(ncsp2) - SUM(nebdust)
        ENDIF

        !set up qpah interpolation
      qlo = MAX(MIN(find_interval(qpaharr,pset%duste_qpah),nqpah_dustem-1),1)
        dq  = (pset%duste_qpah-qpaharr(qlo))/(qpaharr(qlo+1)-qpaharr(qlo))
        dq  = MIN(MAX(dq,0.0),1.0)  !no extrapolation

        !set up Umin interpolation
        !set limits on Umin: 0.1<Umin<25.0
      ulo = MAX(MIN(find_interval(uminarr,pset%duste_umin),numin_dustem-1),1)
        du  = (pset%duste_umin-uminarr(ulo))/(uminarr(ulo+1)-uminarr(ulo))
        du  = MIN(MAX(du,0.0),1.0)  !no extrapolation

        !set limits to gamma (gamma is a fraction)
        gamma = MAX(MIN(pset%duste_gamma,1.0),0.0)

        !bi-linear interpolation over qpah and Umin
        dumin = (1-dq)*(1-du)*dustem2_dustem(:,qlo,2*ulo-1) + &
             dq*(1-du)*dustem2_dustem(:,qlo+1,2*ulo-1) + &
             dq*du*dustem2_dustem(:,qlo+1,2*(ulo+1)-1) + &
             (1-dq)*du*dustem2_dustem(:,qlo,2*(ulo+1)-1)
        dumax = (1-dq)*(1-du)*dustem2_dustem(:,qlo,2*ulo) + &
             dq*(1-du)*dustem2_dustem(:,qlo+1,2*ulo) + &
             dq*du*dustem2_dustem(:,qlo+1,2*(ulo+1)) + &
             (1-dq)*du*dustem2_dustem(:,qlo,2*(ulo+1))

        !combine both parts of P(U)dU
        mduste = (1-gamma)*dumin + gamma*dumax
        mduste = MAX(mduste,tiny_number)

        !normalize the dust emission to the luminosity absorbed by 
        !the dust, i.e., demand that Lbol remains the same
        labs = (lboln-lbold)
         norm   = integrate_trapezoid_array(nu,mduste)
            IF (ieee_is_nan(norm) .OR. norm.LE.tiny_number) THEN
                mdust = tiny_number
                RETURN
            ENDIF
            duste  = mduste/norm * labs
        duste  = MAX(duste,tiny_number)

        !include dust self-absorption
        !we need to iterate because the dust
        !will re-emit and be re-absorbed, etc.
        tduste = 0.0
        DO WHILE (((lboln-lbold).GT.1E-2).OR.iself.EQ.0)

           oduste = duste
           duste  = duste * diff_dust
           tduste = tduste + duste

           lbold = integrate_trapezoid_array(nu,duste)  !after  self-abs
           lboln = integrate_trapezoid_array(nu,oduste) !before self-abs
           IF (ieee_is_nan(lbold)) lbold = 0.0
           IF (ieee_is_nan(lboln)) lboln = 0.0

           duste = MAX(mduste/norm*(lboln-lbold),tiny_number)

           iself=1

        ENDDO

        !this factor assumes Md/Mh=0.01 (appropriate for the 
        !MW3.1 models), and includes conversion factors from 
        !Jy -> Lsun and the hydrogen mass in solar units
        mdust = 3.21E-3 / 4/mypi * labs/norm

        !add the dust emission to the stellar spectrum
        specdust = specdust + tduste

     ELSE
        mdust = tiny_number
     ENDIF

  ENDIF
  END ASSOCIATE

END SUBROUTINE ADD_DUST
