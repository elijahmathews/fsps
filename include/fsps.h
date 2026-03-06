#ifndef FSPS_H
#define FSPS_H

#ifdef __cplusplus
extern "C" {
#endif

// ----------------------------------------------------------------------------
// Status Codes
// ----------------------------------------------------------------------------
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
#define FSPS_ERR_CSP_SCENARIO_FAILED 311
#define FSPS_ERR_CSP_ZDEP_FAILED_0 312
#define FSPS_ERR_CSP_ZDEP_FAILED_1 313
#define FSPS_ERR_CSP_ZDEP_FAILED_2 314
#define FSPS_ERR_CSP_ZDEP_FAILED_3 315
#define FSPS_ERR_LOCK_HELD 401
#define FSPS_ERR_UNLOCK_NOLOCK 402

// ----------------------------------------------------------------------------
// Driver Utilities (Logging, Version, Locking)
// ----------------------------------------------------------------------------
void fsps_set_debug(int flag);
void fsps_clear_error(void);
void fsps_get_last_error(int *status, char *message, int message_len);
void fsps_lock(int *status);
void fsps_unlock(int *status);
void fsps_get_driver_version(int *major, int *minor, int *patch);
void fsps_validate_params(int *status); // Note: Checks legacy global state context

// ----------------------------------------------------------------------------
// Context Management (Thread-Safe API)
// ----------------------------------------------------------------------------
void fsps_context_create(int *handle, int *status);
void fsps_context_destroy(int handle, int *status);
void fsps_context_setup(int zin, const char *isoc_type, const char *spec_type,
                        const char *dust_type, int handle, int *status);

// Parameter Setters
void fsps_context_set_int(int handle, const char *key, int value, int *status);
void fsps_context_set_float(int handle, const char *key, double value, int *status);
void fsps_context_set_str(int handle, const char *key, const char *value, int *status);
void fsps_context_set_fast_mode(int handle, int fast_mode, int *status);

// Dimension Queries
void fsps_context_get_dims(int handle, int *n_spec, int *n_time, int *status);
void fsps_context_get_nspec(int handle, int *n_spec, int *status);
void fsps_context_get_ntfull(int handle, int *n_time, int *status);
void fsps_context_get_nbands(int handle, int *n_bands, int *status);
void fsps_context_get_nindx(int handle, int *n_indices, int *status);
void fsps_context_get_nz(int handle, int *n_z, int *status);
void fsps_context_get_nemline(int handle, int *n_line, int *status);
void fsps_context_get_zsol(int handle, double *z_sol, int *status);

// Path Queries
void fsps_context_get_paths(int handle, char *sps_home, int sps_len,
                            char *data_home, int data_len, char *output_home, int out_len);

// Context Computation
void fsps_context_compute_ssp(int handle, double *spec, double *mass, double *lbol, int *status);

// ----------------------------------------------------------------------------
// Legacy Global State API (Backward Compatibility)
// ----------------------------------------------------------------------------

// Initialization
#ifdef FSPS_ENABLE_LEGACY
void fsps_initialize(int zin);
void fsps_initialize_full(int zin, int compute_vega_mags0, int vactoair_flag0,
                          const char *isoc_type, const char *spec_type, const char *dust_type);
#endif

// Parameter Setters (Global)
void fsps_set_int(const char *key, int value);
void fsps_set_float(const char *key, double value);
void fsps_set_str(const char *key, const char *value);

// Specific Setters
void fsps_set_sfh_tab(int ntab, double *age, double *sfr, double *met);
void fsps_set_mag_compute(int n_bands, int *mask);
void fsps_set_ssp_gen_age(int n_age, int *mask);
void fsps_set_ssp_lsf(int nsv, double *sigma, double wlo, double whi);

// Computation (Operates on Global State)
void fsps_compute(double *spec_out); // Main driver (SSP or CSP based on `sfh`)
void fsps_compute_csp(int zcontinuous);
void fsps_compute_ssp(int zin);
void fsps_compute_ssps(void);
void fsps_compute_zdep(int ztype);

// Interpolation
void fsps_interp_ssp(double zpos, double tpos, double *spec, double *mass, double *lbol);

// ----------------------------------------------------------------------------
// Data Access & Post-Processing
// ----------------------------------------------------------------------------

// Photometry
void fsps_get_mags(double zred, double *mags);
void fsps_get_mags_mask(double zred, double *mags, int *mask);

// Spectral Data
void fsps_get_spec(double *spec);
void fsps_get_spec_peraa(double *spec);
void fsps_get_indices(double *spec, double *indices);
void fsps_smooth_spectrum(double *wave, double *spec, double sigma_broad, double minw, double maxw);

// Physical Properties
void fsps_get_stats(double *age, double *mass, double *lbol, double *sfr, 
                    double *mdust, double *mformed, double *emlines);
void fsps_get_ssp_weights(double *wghts);
void fsps_get_csp_components(double *young, double *old);
void fsps_get_ssp_spec(double *spec, double *mass, double *lbol);

// Stellar Utility
void fsps_stellar_spectrum(double mact, double logt, double lbol, double logg, 
                           double phase, double ffco, double lmdot, double wght, double *spec);

// Metadata / Arrays
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

// Array Getters
void fsps_get_zsol(double *z_sol);
void fsps_get_zlegend(double *zlegend);
void fsps_get_timefull(double *timefull);
void fsps_get_lambda(double *lambda);
void fsps_get_emlambda(double *emlambda);
void fsps_get_res(double *res);
void fsps_get_filter_data(double *wave_eff, double *mag_vega, double *mag_sun);
#ifdef FSPS_ENABLE_LEGACY
void fsps_get_setup_vars(int *compute_vega_mags, int *vactoair_flag);
void fsps_get_libraries(char *isoc, int isoc_len, char *spec, int spec_len, char *dust, int dust_len);
#endif

// Isochrone Specifics
void fsps_get_isochrone_dimensions(int *n_age, int *n_mass);
void fsps_get_nmass_isochrone(int z_idx, int t_idx, int *n_mass);
void fsps_write_isochrone(const char *outfile);

// Cleanup
void fsps_finalize(void);

#ifdef __cplusplus
}
#endif

#endif // FSPS_H