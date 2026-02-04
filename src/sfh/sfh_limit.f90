function sfhlimit(ctx, tlim, sfh)
  ! Deal with all the logic for computing the limits of integration in the
  ! weight calculations.  This should return the input `tlim`, clipped to valid
  ! limits, such that if no SFR occurs between t1 and t2 then sfhlimit(t1, sfh)
  ! = sfhlimit(t2, sfh)
  !
  ! Inputs
  ! ----------
  !
  ! tlim:
  !    The proposed limit of the integration, in log(lookback_time(yrs)).
  !
  ! sfh:
  !    An sfhparams structure containing the relevant special lookback times.
  !
   use fsps_context_types, only: fsps_context_t
   use fsps_constants, only: SP
   use fsps_types, only: SFHPARAMS
  implicit none

   type(fsps_context_t), intent(in) :: ctx
  real(SP), intent(in) :: tlim
  type(SFHPARAMS), intent(in) :: sfh

  real(SP) :: sfhlimit

  real(SP) :: tlo, thi

  ! For the simha linear portion, we integrate from sf_trunc to tage or the
  ! zero crossing, whichever is smaller but still greater than sf_trunc.
  ! For everything else we integrate from 0 to tage or sf_trunc, whichever is
  ! smaller but non-zero.
  ! Of course, we are doing this using the lookback time versions of sf_trunc,
  ! the zero crossing time, and 0!!!
  if (sfh%use_simha_limits.eq.1) then
     tlo = sfh%t0
     thi = sfh%tq
  else
     tlo = sfh%tq
     thi = sfh%tage
  endif

  ! Convert to log, taking care of possible zeros.
   tlo = log10(max(tlo, 10**ctx%tiny_logt_val))
   thi = log10(max(thi, 10**ctx%tiny_logt_val))

  ! Finally, clip proposed limit to the upper and lower value
  sfhlimit = MIN(MAX(tlim, tlo), thi)

end function sfhlimit
