# csp_gen.f90: with use_given_weights = 1 the caller's SSP weights (mass formed per SSP age,
# given_weights) replace the weights computed from the SFH parameters; the SSP sum, the young /
# old split at dust_tesc, add_dust and everything after are unchanged.
s|^\( *weight_ssp, spec_young, spec_old\)$|\1, given_weights, use_given_weights|
/^  spec_young = 0\. ! young/i\
  if (use_given_weights.eq.1) then\
     total_weights = 0.\
     total_weights(:, 1) = given_weights\
     imin = 0\
     imax = ntfull\
  endif\

