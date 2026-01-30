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

// Error handling and debug logging
void fsps_set_debug(int flag);
void fsps_clear_error(void);
void fsps_get_last_error(int *status, char *message, int message_len);
void fsps_lock(int *status);
void fsps_unlock(int *status);
void fsps_get_driver_version(int *major, int *minor, int *patch);

// Context-based API (handle-oriented)
void fsps_context_create(int *handle, int *status);
void fsps_context_destroy(int handle, int *status);
void fsps_context_setup(int zin, const char *isoc_type, const char *spec_type,
						const char *dust_type, int handle, int *status);
void fsps_context_set_int(int handle, const char *key, int value, int *status);
void fsps_context_set_float(int handle, const char *key, double value, int *status);
void fsps_context_set_str(int handle, const char *key, const char *value, int *status);
void fsps_context_compute_ssp(int handle, double *spec, double *mass, double *lbol, int *status);
void fsps_context_get_paths(int handle, char *sps_home, int sps_len,
							char *data_home, int data_len, char *output_home, int out_len);
// Context dimension queries
void fsps_context_get_dims(int handle, int *n_spec, int *n_time, int *status);
void fsps_context_get_nspec(int handle, int *n_spec, int *status);
void fsps_context_get_ntfull(int handle, int *n_time, int *status);
void fsps_context_get_nbands(int handle, int *n_bands, int *status);
void fsps_context_get_nindx(int handle, int *n_indices, int *status);
void fsps_context_get_nz(int handle, int *n_z, int *status);
void fsps_context_get_nemline(int handle, int *n_line, int *status);
void fsps_context_get_zsol(int handle, double *z_sol, int *status);

#ifdef __cplusplus
}
#endif

#endif // FSPS_H