FUNCTION AGN_DUST(ctx, lam, spec, pset, lbol_csp)

     USE fsps_context_types, ONLY: fsps_context_t
     USE fsps_types, ONLY: SP, PARAMS, nagndust
     USE fsps_dust, ONLY: compute_attenuation_curve
     USE fsps_interpolation, ONLY: find_interval
     IMPLICIT NONE

     TYPE(fsps_context_t), INTENT(INOUT) :: ctx
     REAL(SP), DIMENSION(:), INTENT(in) :: lam,spec
     REAL(SP), INTENT(in)       :: lbol_csp
     TYPE(PARAMS), INTENT(in)   :: pset
     REAL(SP), DIMENSION(SIZE(lam)) :: agn_dust,agnspeci
     INTEGER  :: jlo
     REAL(SP) :: dj
 
  !--------------------------------------------------------------!

  ASSOCIATE(agndust_tau => ctx%state%agndust_tau, &
            agndust_spec => ctx%state%agndust_spec, &
            dust_type => ctx%dust_type_val)

    !interpolate in tau_agn
     jlo = MIN(MAX(find_interval(agndust_tau,pset%agn_tau),1),&
         nagndust-1)
    dj  = (pset%agn_tau-agndust_tau(jlo)) / &
         (agndust_tau(jlo+1)-agndust_tau(jlo))
    dj  = MAX(MIN(dj,1.0),0.0) !no extrapolation

    agnspeci  = (1-dj)*agndust_spec(:,jlo) + dj*agndust_spec(:,jlo+1)

    !attenuate the AGN emission by the diffuse dust
     agnspeci = agnspeci*EXP(-compute_attenuation_curve(lam, dust_type, pset, ctx))

    agn_dust = spec + 10**lbol_csp*pset%fagn*agnspeci
  END ASSOCIATE


END FUNCTION AGN_DUST
  
