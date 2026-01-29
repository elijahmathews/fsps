#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>
#include "fsps.h"

// Helper to check for errors
void check_alloc(void *ptr) {
    if (ptr == NULL) {
        fprintf(stderr, "[C] Error: Memory allocation failed!\n");
        exit(1);
    }
}

int main() {
    printf("--- Comprehensive C Driver Test for FSPS ---\n");

    // ---------------------------------------------------------
    // 1. Initialization
    // ---------------------------------------------------------
    int metallicity_index = 10;
    printf("[C] Initializing FSPS with ZIN=%d...\n", metallicity_index);
    fsps_initialize_full(metallicity_index, 0, 0, "mist", "miles", "DL07");

    // ---------------------------------------------------------
    // 2. Query Dimensions
    // ---------------------------------------------------------
        int nspec = 0, ntfull = 0, nbands = 0, nz = 0, nemline = 0, nindx = 0;
        int nspec2 = 0, ntfull2 = 0;
        double zsol = 0.0;
    fsps_get_dims(&nspec, &ntfull);
    fsps_get_nbands(&nbands);
        fsps_get_nindx(&nindx);
        fsps_get_nspec(&nspec2);
        fsps_get_ntfull(&ntfull2);
        fsps_get_nz(&nz);
        fsps_get_nemline(&nemline);
        fsps_get_zsol(&zsol);
    
            printf("[C] Dimensions: nspec=%d, ntfull=%d, nbands=%d, nz=%d, nemline=%d, nindx=%d\n",
                nspec, ntfull, nbands, nz, nemline, nindx);
        printf("[C] Dimensions (direct): nspec=%d, ntfull=%d\n", nspec2, ntfull2);
        printf("[C] zsol=%g\n", zsol);

    if (nspec <= 0 || ntfull <= 0 || nbands <= 0) {
        fprintf(stderr, "[C] Error: Invalid dimensions.\n");
        return 1;
    }

    size_t spec_elements = (size_t)nspec * (size_t)ntfull;
    size_t mag_elements  = (size_t)nbands * (size_t)ntfull;
    size_t emline_elements = (size_t)nemline * (size_t)ntfull;

    double *spec_array = (double*)calloc(spec_elements, sizeof(double));
    double *mags_array = (double*)calloc(mag_elements, sizeof(double));
    double *mags_mask_array = (double*)calloc(mag_elements, sizeof(double));
    double *spec_peraa = (double*)calloc(spec_elements, sizeof(double));
    check_alloc(spec_array);
    check_alloc(mags_array);
    check_alloc(mags_mask_array);
    check_alloc(spec_peraa);

    double *spec_out = (double*)calloc(spec_elements, sizeof(double));
    double *age = (double*)calloc(ntfull, sizeof(double));
    double *mass = (double*)calloc(ntfull, sizeof(double));
    double *lbol = (double*)calloc(ntfull, sizeof(double));
    double *sfr = (double*)calloc(ntfull, sizeof(double));
    double *mdust = (double*)calloc(ntfull, sizeof(double));
    double *mformed = (double*)calloc(ntfull, sizeof(double));
    double *emlines = (double*)calloc(emline_elements, sizeof(double));
    check_alloc(spec_out);
    check_alloc(age);
    check_alloc(mass);
    check_alloc(lbol);
    check_alloc(sfr);
    check_alloc(mdust);
    check_alloc(mformed);
    check_alloc(emlines);

    double *zlegend = (double*)calloc(nz, sizeof(double));
    double *timefull = (double*)calloc(ntfull, sizeof(double));
    double *lambda = (double*)calloc(nspec, sizeof(double));
    double *emlambda = (double*)calloc(nemline, sizeof(double));
    double *res = (double*)calloc(nspec, sizeof(double));
    double *wave_eff = (double*)calloc(nbands, sizeof(double));
    double *mag_vega = (double*)calloc(nbands, sizeof(double));
    double *mag_sun = (double*)calloc(nbands, sizeof(double));
    double *ssp_weights = (double*)calloc((size_t)ntfull * (size_t)nz, sizeof(double));
    double *spec_young = (double*)calloc(nspec, sizeof(double));
    double *spec_old = (double*)calloc(nspec, sizeof(double));
    double *spec_star = (double*)calloc(nspec, sizeof(double));
    double *spec_smoothed = (double*)calloc(nspec, sizeof(double));
    double *indices = (double*)calloc(nindx, sizeof(double));
    check_alloc(zlegend);
    check_alloc(timefull);
    check_alloc(lambda);
    check_alloc(emlambda);
    check_alloc(res);
    check_alloc(wave_eff);
    check_alloc(mag_vega);
    check_alloc(mag_sun);
    check_alloc(ssp_weights);
    check_alloc(spec_young);
    check_alloc(spec_old);
    check_alloc(spec_star);
    check_alloc(spec_smoothed);
    check_alloc(indices);

    // ---------------------------------------------------------
    // 3. Test SSP Generation (SFH=0)
    // ---------------------------------------------------------
    printf("\n[C] --- TEST 1: SSP Generation (SFH=0) ---\n");
    
    // Set a parameter to be sure
    fsps_set_int("sfh", 0);
    fsps_set_int("imf_type", 2); // Kroupa
    fsps_set_float("imf_upper_limit", 100.0);
    fsps_set_float("imf_lower_limit", 0.1);
    fsps_set_str("imf_filename", "imf.dat");

    printf("[C] Computing SSPs...\n");
    fsps_compute(spec_array);
    
    printf("[C] Spec[0]      = %e\n", spec_array[0]);
    printf("[C] Spec[last]   = %e\n", spec_array[spec_elements - 1]);

    printf("[C] Fetching spectra and stats...\n");
    fsps_get_spec(spec_out);
    fsps_get_spec_peraa(spec_peraa);
    if (nindx > 0) {
        fsps_get_indices(spec_out, indices);
        printf("[C] Index[0]    = %f\n", indices[0]);
    }
    fsps_get_stats(age, mass, lbol, sfr, mdust, mformed, emlines);
    printf("[C] Spec_out[0]  = %e\n", spec_out[0]);
    printf("[C] Spec_perAA[0] = %e\n", spec_peraa[0]);
    printf("[C] Age[0]       = %e\n", age[0]);
    printf("[C] Emline[0,0]  = %e\n", emlines[0]);

    // ---------------------------------------------------------
    // 4. Test CSP Generation (SFH=1)
    // ---------------------------------------------------------
    printf("\n[C] --- TEST 2: CSP Generation (SFH=1, Tau Model) ---\n");
    
    fsps_set_int("sfh", 1);
    fsps_set_float("tau", 1.0);
    fsps_set_float("const", 0.0);
    fsps_set_float("sf_start", 0.0);
    fsps_set_float("tage", 13.7); // Age of universe roughly
    
    // CORRECTION: dust_type is an integer, so we must use set_int
    fsps_set_int("dust_type", 2); // Calzetti
    fsps_set_float("dust2", 0.5); // Some attenuation

    printf("[C] Computing CSP (via COMPSP)...\n");
    // Clear array to be sure we get new data
    for(size_t i=0; i<spec_elements; i++) spec_array[i] = 0.0;
    
    fsps_compute(spec_array);

    printf("[C] Spec[0]      = %e\n", spec_array[0]);
    printf("[C] Spec[last]   = %e\n", spec_array[spec_elements - 1]);

    // ---------------------------------------------------------
    // 5. Test Magnitude Retrieval
    // ---------------------------------------------------------
    printf("\n[C] --- TEST 3: Photometry ---\n");
    
    double redshift = 0.0;
    printf("[C] computing magnitudes at z=%.1f...\n", redshift);
    fsps_get_mags(redshift, mags_array);

    // Check first band (usually V or something standard depending on filters.dat) at first time step
    printf("[C] Mag[Band 0, Time 0]   = %f\n", mags_array[0]);
    // Check last band at last time step
    printf("[C] Mag[Band %d, Time %d] = %f\n", nbands-1, ntfull-1, mags_array[mag_elements-1]);

    // ---------------------------------------------------------
    // 6. Test Metadata and Helpers
    // ---------------------------------------------------------
        printf("\n[C] --- TEST 4: Metadata & Helpers ---\n");

        int param_status = 0;
        int err_status = 0;
        char err_msg[256] = {0};
        int lock_status = 0;
        int ver_major = 0, ver_minor = 0, ver_patch = 0;
        int nt = 0, nm = 0, ntabmax = 0;
        fsps_validate_params(&param_status);
        printf("[C] validate_params status=%d\n", param_status);

        fsps_get_driver_version(&ver_major, &ver_minor, &ver_patch);
        printf("[C] driver_version=%d.%d.%d\n", ver_major, ver_minor, ver_patch);

        fsps_get_nt(&nt);
        fsps_get_nm(&nm);
        fsps_get_ntabmax(&ntabmax);
        printf("[C] nt=%d nm=%d ntabmax=%d\n", nt, nm, ntabmax);

        fsps_lock(&lock_status);
        printf("[C] lock_status=%d\n", lock_status);
        fsps_unlock(&lock_status);
        printf("[C] unlock_status=%d\n", lock_status);

        fsps_clear_error();
        fsps_get_last_error(&err_status, err_msg, (int)sizeof(err_msg));
        printf("[C] last_error status=%d msg=%s\n", err_status, err_msg);

        int n_iso_age = 0, n_iso_mass = 0;
        int nmass_iso = 0;
        fsps_get_isochrone_dimensions(&n_iso_age, &n_iso_mass);
        fsps_get_nmass_isochrone(1, 1, &nmass_iso);
        printf("[C] isochrone dims: n_age=%d, n_mass=%d, nmass(1,1)=%d\n",
            n_iso_age, n_iso_mass, nmass_iso);

    fsps_get_zlegend(zlegend);
    fsps_get_timefull(timefull);
    fsps_get_lambda(lambda);
    fsps_get_emlambda(emlambda);
    fsps_get_res(res);
    fsps_get_filter_data(wave_eff, mag_vega, mag_sun);
    fsps_get_ssp_weights(ssp_weights);
    fsps_get_csp_components(spec_young, spec_old);

    printf("[C] zlegend[0]   = %e\n", zlegend[0]);
    printf("[C] timefull[0]  = %e\n", timefull[0]);
    printf("[C] lambda[0]    = %e\n", lambda[0]);
    printf("[C] emlambda[0]  = %e\n", emlambda[0]);
    printf("[C] res[0]       = %e\n", res[0]);
    printf("[C] wave_eff[0]  = %e\n", wave_eff[0]);
    printf("[C] mag_vega[0]  = %e\n", mag_vega[0]);
    printf("[C] mag_sun[0]   = %e\n", mag_sun[0]);
    printf("[C] weights[0]   = %e\n", ssp_weights[0]);
    printf("[C] young[0]     = %e\n", spec_young[0]);
    printf("[C] old[0]       = %e\n", spec_old[0]);

    int cvms = 0, vta = 0;
    char isoc_buf[32] = {0};
    char spec_buf[32] = {0};
    char dust_buf[32] = {0};
    fsps_get_setup_vars(&cvms, &vta);
    fsps_get_libraries(isoc_buf, (int)sizeof(isoc_buf), spec_buf, (int)sizeof(spec_buf),
                       dust_buf, (int)sizeof(dust_buf));
    printf("[C] compute_vega_mags=%d, vactoair_flag=%d\n", cvms, vta);
    printf("[C] libraries: isoc=%s, spec=%s, dust=%s\n", isoc_buf, spec_buf, dust_buf);

    // Set mag_compute and ssp_gen_age masks
    int *mag_mask = (int*)calloc(nbands, sizeof(int));
    int *age_mask = (int*)calloc(n_iso_age, sizeof(int));
    check_alloc(mag_mask);
    check_alloc(age_mask);
    for (int i = 0; i < nbands; i++) mag_mask[i] = (i < 5) ? 1 : 0;
    for (int i = 0; i < n_iso_age; i++) age_mask[i] = (i == 0) ? 1 : 0;
    fsps_set_mag_compute(nbands, mag_mask);
    fsps_set_ssp_gen_age(n_iso_age, age_mask);
    free(mag_mask);
    free(age_mask);

    // ---------------------------------------------------------
    // 7. Test SSP helpers, interpolation, and smoothing
    // ---------------------------------------------------------
    printf("\n[C] --- TEST 5: SSP Helpers ---\n");
    fsps_compute_ssp(metallicity_index);

    if (metallicity_index == -1) {
        double *interp_spec = (double*)calloc(nspec, sizeof(double));
        double interp_mass[1] = {0.0};
        double interp_lbol[1] = {0.0};
        check_alloc(interp_spec);
        fsps_interp_ssp(zlegend[0], timefull[0], interp_spec, interp_mass, interp_lbol);
        printf("[C] interp_spec[0] = %e\n", interp_spec[0]);
        free(interp_spec);
    } else {
        printf("[C] Skipping fsps_interp_ssp (requires full Z grid; use ZIN=-1)\n");
    }

    // Smooth a spectrum (copy first spectrum)
    for (int i = 0; i < nspec; i++) {
        spec_smoothed[i] = spec_array[i];
    }
    fsps_smooth_spectrum(lambda, spec_smoothed, 100.0, lambda[0], lambda[nspec - 1]);
    printf("[C] smoothed[0]  = %e\n", spec_smoothed[0]);

    // Stellar spectrum at a representative point
    fsps_stellar_spectrum(1.0, 3.7, 1.0, 4.5, 3.0, 0.0, 0.0, 1.0, spec_star);
    printf("[C] star_spec[0] = %e\n", spec_star[0]);

    // Tabular SFH and LSF setters
    double tab_age[3] = {0.1, 1.0, 10.0};
    double tab_sfr[3] = {1.0, 0.5, 0.1};
    double tab_met[3] = {zlegend[0], zlegend[0], zlegend[0]};
    fsps_set_sfh_tab(3, tab_age, tab_sfr, tab_met);

    double lsf_sigma[3] = {50.0, 60.0, 70.0};
    fsps_set_ssp_lsf(3, lsf_sigma, lambda[0], lambda[nspec - 1]);

    // Optional: retrieve full SSP grids if size is reasonable
    if (metallicity_index == -1) {
        size_t ssp_spec_elements = (size_t)nspec * (size_t)ntfull * (size_t)nz;
        size_t ssp_spec_bytes = ssp_spec_elements * sizeof(double);
        if (ssp_spec_bytes < (size_t)200 * 1024 * 1024) {
            double *ssp_spec = (double*)calloc(ssp_spec_elements, sizeof(double));
            double *ssp_mass = (double*)calloc((size_t)ntfull * (size_t)nz, sizeof(double));
            double *ssp_lbol = (double*)calloc((size_t)ntfull * (size_t)nz, sizeof(double));
            check_alloc(ssp_spec);
            check_alloc(ssp_mass);
            check_alloc(ssp_lbol);

            fsps_get_ssp_spec(ssp_spec, ssp_mass, ssp_lbol);
            printf("[C] ssp_spec[0] = %e\n", ssp_spec[0]);

            free(ssp_spec);
            free(ssp_mass);
            free(ssp_lbol);
        } else {
            printf("[C] Skipping fsps_get_ssp_spec (%.1f MB too large)\n",
                   (double)ssp_spec_bytes / (1024.0 * 1024.0));
        }
    } else {
        printf("[C] Skipping fsps_get_ssp_spec (requires full Z grid; use ZIN=-1)\n");
    }

    // CSP z-dependence
    fsps_compute_csp(0);

    // Write isochrone (smoke test)
    fsps_write_isochrone("test_isochrone.out");

    // Masked photometry
    int *mc = (int*)calloc(nbands, sizeof(int));
    check_alloc(mc);
    for (int i = 0; i < nbands; i++) mc[i] = (i < 3) ? 1 : 0;
    fsps_get_mags_mask(0.0, mags_mask_array, mc);
    printf("[C] Masked Mag[Band 0, Time 0]   = %f\n", mags_mask_array[0]);
    printf("[C] Masked Mag[Band 2, Time 0]   = %f\n", mags_mask_array[2]);
    free(mc);

    // ---------------------------------------------------------
    // 8. Library selection API
    // ---------------------------------------------------------
    printf("\n[C] --- TEST 6: Library selection API ---\n");
    fsps_finalize();
    fsps_initialize_full(metallicity_index, 0, 0, "", "", "");
    fsps_get_libraries(isoc_buf, (int)sizeof(isoc_buf), spec_buf, (int)sizeof(spec_buf),
                       dust_buf, (int)sizeof(dust_buf));
    printf("[C] libraries (default): isoc=%s, spec=%s, dust=%s\n", isoc_buf, spec_buf, dust_buf);

    fsps_finalize();
    fsps_initialize_full(metallicity_index, 0, 0, "prsc", "basel", "THEMIS");
    fsps_get_libraries(isoc_buf, (int)sizeof(isoc_buf), spec_buf, (int)sizeof(spec_buf),
                       dust_buf, (int)sizeof(dust_buf));
    printf("[C] libraries (custom): isoc=%s, spec=%s, dust=%s\n", isoc_buf, spec_buf, dust_buf);

    // ---------------------------------------------------------
    // 9. Full Z grid test (optional)
    // ---------------------------------------------------------
    const char *fullz = getenv("FSPS_FULL_Z_TEST");
    if (fullz && (strcmp(fullz, "1") == 0 || strcasecmp(fullz, "true") == 0)) {
        printf("\n[C] --- TEST 7: Full Z grid ---\n");
        fsps_finalize();
        fsps_initialize_full(-1, 0, 0, "mist", "miles", "DL07");

        int nz_full = 0;
        fsps_get_nz(&nz_full);
        printf("[C] Full Z grid nz=%d\n", nz_full);

        fsps_compute_ssps();
        size_t ssp_spec_elements = (size_t)nspec * (size_t)ntfull * (size_t)nz_full;
        size_t ssp_spec_bytes = ssp_spec_elements * sizeof(double);
        if (ssp_spec_bytes < (size_t)200 * 1024 * 1024) {
            double *ssp_spec = (double*)calloc(ssp_spec_elements, sizeof(double));
            double *ssp_mass = (double*)calloc((size_t)ntfull * (size_t)nz_full, sizeof(double));
            double *ssp_lbol = (double*)calloc((size_t)ntfull * (size_t)nz_full, sizeof(double));
            check_alloc(ssp_spec);
            check_alloc(ssp_mass);
            check_alloc(ssp_lbol);

            fsps_get_ssp_spec(ssp_spec, ssp_mass, ssp_lbol);
            printf("[C] fullZ ssp_spec[0] = %e\n", ssp_spec[0]);

            free(ssp_spec);
            free(ssp_mass);
            free(ssp_lbol);
        } else {
            printf("[C] Skipping full-Z get_ssp_spec (%.1f MB too large)\n",
                   (double)ssp_spec_bytes / (1024.0 * 1024.0));
        }
    } else {
        printf("\n[C] --- TEST 7: Full Z grid ---\n");
        printf("[C] Skipped (set FSPS_FULL_Z_TEST=1 to enable)\n");
    }

    // ---------------------------------------------------------
    // 10. Cleanup
    // ---------------------------------------------------------
    printf("\n[C] Finalizing FSPS...\n");
    fsps_finalize();
    
    free(spec_array);
    free(mags_array);
    free(mags_mask_array);
    free(spec_peraa);
    free(spec_out);
    free(age);
    free(mass);
    free(lbol);
    free(sfr);
    free(mdust);
    free(mformed);
    free(emlines);
    free(zlegend);
    free(timefull);
    free(lambda);
    free(emlambda);
    free(res);
    free(wave_eff);
    free(mag_vega);
    free(mag_sun);
    free(ssp_weights);
    free(spec_young);
    free(spec_old);
    free(spec_star);
    free(spec_smoothed);
    free(indices);
    
    printf("[C] Success.\n");

    return 0;
}