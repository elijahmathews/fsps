#ifndef FSPS_H
#define FSPS_H

#ifdef __cplusplus
extern "C" {
#endif

// Driver error/status codes
#define FSPS_STATUS_OK 0
#define FSPS_ERR_UNKNOWN_INT_PARAM 101
#define FSPS_ERR_UNKNOWN_FLOAT_PARAM 102
#define FSPS_ERR_UNKNOWN_STRING_PARAM 103
#define FSPS_WARN_ZMET_RANGE 201
#define FSPS_WARN_DUST_RANGE 202
#define FSPS_WARN_IMF_RANGE 203
#define FSPS_WARN_SF_START_GT_TAGE 204
#define FSPS_WARN_CONST_FBURST 205
#define FSPS_WARN_ZCONTINUOUS_NEB 206
#define FSPS_WARN_NMASS_RANGE 207
#define FSPS_WARN_MAG_MASK_RANGE 208
#define FSPS_WARN_SSP_AGE_RANGE 209
#define FSPS_ERR_COMPUTE_UNINIT 301
#define FSPS_ERR_CSP_UNINIT 302
#define FSPS_ERR_SSP_UNINIT 303
#define FSPS_ERR_SSP_Z_RANGE 304
#define FSPS_ERR_ZDEP_UNINIT 305
#define FSPS_ERR_ZDEP_TYPE 306
#define FSPS_ERR_INTERP_UNINIT 307
#define FSPS_ERR_STELLAR_UNINIT 308
#define FSPS_ERR_MAG_MASK_UNINIT 309
#define FSPS_ERR_SSP_AGE_UNINIT 310
#define FSPS_ERR_LOCK_HELD 401
#define FSPS_ERR_UNLOCK_NOLOCK 402

// Initialization and Teardown
// zin = -1 loads the full metallicity grid; otherwise only a single Z bin.
void fsps_initialize(int zin);
void fsps_initialize_full(int zin, int compute_vega_mags, int vactoair_flag,
						  const char *isoc_type, const char *spec_type, const char *dust_type);
void fsps_finalize();

// Dimension Queries
void fsps_get_dims(int *n_spec, int *n_time);
void fsps_get_nbands(int *n_bands);
void fsps_get_nindx(int *n_indices);
void fsps_get_nspec(int *n_spec);
void fsps_get_ntfull(int *n_time);
void fsps_get_nt(int *n_time);
void fsps_get_nm(int *n_mass);
void fsps_get_ntabmax(int *n_tabmax);
void fsps_get_nz(int *n_z);
void fsps_get_nemline(int *n_line);
void fsps_get_zsol(double *z_sol);
void fsps_get_isochrone_dimensions(int *n_age, int *n_mass);
void fsps_get_nmass_isochrone(int z_idx, int t_idx, int *n_mass);

// Parameter Control
// keys are strings (e.g., "imf_type", "dust_type", "sfh")
void fsps_set_int(const char *key, int value);
void fsps_set_float(const char *key, double value);
void fsps_set_str(const char *key, const char *value);
void fsps_set_mag_compute(int n_bands, const int *mask);
void fsps_set_ssp_gen_age(int n_age, const int *mask);
void fsps_validate_params(int *status);
// Error handling and debug logging
void fsps_set_debug(int flag);
void fsps_clear_error(void);
void fsps_get_last_error(int *status, char *message, int message_len);
void fsps_lock(int *status);
void fsps_unlock(int *status);
void fsps_get_driver_version(int *major, int *minor, int *patch);

// Computation
// spec must be a pre-allocated array of size (n_spec * n_time)
void fsps_compute(double *spec);
void fsps_compute_ssp(int zin);
void fsps_compute_ssps();
void fsps_compute_zdep(int ztype);
void fsps_interp_ssp(double zpos, double tpos, double *spec, double *mass, double *lbol);
void fsps_compute_csp(int zcontinuous);

// Photometry
// mags must be a pre-allocated array of size (n_bands * n_time)
// zred is the redshift to place the spectrum at (does not affect IGM unless zred > 0 in setup)
void fsps_get_mags(double zred, double *mags);
void fsps_get_mags_mask(double zred, double *mags, const int *mc);

// Outputs and metadata
void fsps_get_spec(double *spec);
void fsps_get_spec_peraa(double *spec);
void fsps_get_indices(const double *spec, double *indices);
void fsps_get_stats(double *age, double *mass, double *lbol, double *sfr,
					double *mdust, double *mformed, double *emlines);
void fsps_stellar_spectrum(double mact, double logt, double lbol, double logg,
						   double phase, double ffco, double lmdot, double wght,
						   double *spec);
void fsps_get_zlegend(double *zlegend);
void fsps_get_timefull(double *timefull);
void fsps_get_lambda(double *lambda);
void fsps_get_emlambda(double *emlambda);
void fsps_get_res(double *res);
void fsps_get_filter_data(double *wave_eff, double *mag_vega, double *mag_sun);
void fsps_get_ssp_weights(double *weights);
void fsps_get_csp_components(double *spec_young, double *spec_old);
void fsps_get_ssp_spec(double *spec, double *mass, double *lbol);

// Tabular SFH and LSF control
void fsps_set_sfh_tab(int ntab, const double *age, const double *sfr, const double *met);
void fsps_set_ssp_lsf(int nsv, const double *sigma, double wlo, double whi);
void fsps_smooth_spectrum(const double *wave, double *spec, double sigma_broad, double minw, double maxw);

// Misc
void fsps_write_isochrone(const char *outfile);
void fsps_get_setup_vars(int *compute_vega_mags, int *vactoair_flag);
void fsps_get_libraries(char *isoc, int isoc_len, char *spec, int spec_len, char *dust, int dust_len);

#ifdef __cplusplus
}
#endif

#endif // FSPS_H