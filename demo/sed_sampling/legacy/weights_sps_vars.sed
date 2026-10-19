# sps_vars.f90: the caller's SSP weights for csp_gen (weights_csp_gen.sed)
/^  REAL(SP), DIMENSION(ntfull,nz) *:: weight_ssp=0\.$/a\
  REAL(SP), DIMENSION(ntfull) :: given_weights=0.\
  INTEGER :: use_given_weights=0
