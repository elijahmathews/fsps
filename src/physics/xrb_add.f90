SUBROUTINE ADD_XRB(ctx, pset, sspi, sspo)

  ! Routine to add emission from X-ray binaries

  USE fsps_context_types, ONLY: fsps_context_t
  USE fsps_types, ONLY: SP, PARAMS
  USE fsps_interpolation, ONLY: find_interval
  IMPLICIT NONE

  TYPE(fsps_context_t), INTENT(INOUT) :: ctx
  INTEGER :: t,a1,z1
  REAL(SP) :: da,dz,tmpz
  TYPE(PARAMS), INTENT(in) :: pset
  REAL(SP), INTENT(in), DIMENSION(:,:)    :: sspi
  REAL(SP), INTENT(inout), DIMENSION(:,:) :: sspo
  REAL(SP), DIMENSION(SIZE(sspi,1)) :: tmpspec

  !-----------------------------------------------------------!
  !-----------------------------------------------------------!

    ASSOCIATE( &
      zlegend => ctx%state%zlegend, zsol => ctx%state%zsol, &
      zmet_xrb => ctx%state%zmet_xrb, nz_xrb => ctx%state%nz_xrb, &
      ages_xrb => ctx%state%ages_xrb, nt_xrb => ctx%state%nt_xrb, &
      time_full => ctx%state%time_full, spec_xrb => ctx%state%spec_xrb, &
      nt => ctx%state%nt )

    !set up the interpolation variables for logZ
  tmpz = log10(zlegend(pset%zmet)/zsol)
  z1   = MAX(MIN(find_interval(zmet_xrb,tmpz),nz_xrb-1),1)
  dz   = (tmpz-zmet_xrb(z1))/(zmet_xrb(z1+1)-zmet_xrb(z1))
  dz   = MAX(MIN(dz,1.0),0.0) !no extrapolation

  sspo = sspi
 
  DO t=1,nt

     !set up age interpolant
    a1 = MAX(MIN(find_interval(ages_xrb,time_full(t)),nt_xrb-1),1)
     da = (time_full(t)-ages_xrb(a1))/(ages_xrb(a1+1)-ages_xrb(a1))
     
     IF (da.LT.0.0.OR.da.GT.1.0) CYCLE

     tmpspec = &   !interpolate in logZ and time
          (1-da)*(1-dz)* spec_xrb(:,a1,z1)+&
          da*(1-dz)* spec_xrb(:,a1+1,z1)+&
          (1-da)*dz* spec_xrb(:,a1,z1+1)+&
          da*dz* spec_xrb(:,a1+1,z1+1)

     sspo(:,t) = sspo(:,t) + pset%frac_xrb * tmpspec

  ENDDO

  END ASSOCIATE


END SUBROUTINE ADD_XRB
