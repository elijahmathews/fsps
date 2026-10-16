/* The C batch entry point (fsps_context_compute_csp_batch): parameter arrays (one per
 * parameter, length n) against the same models set one at a time on the context, and the error
 * statuses for a parameter that changes the SSPs and an unknown one. Requires the FSPS data tree
 * (SPS_HOME), like the regression tests. */
#include "fsps.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define N_MOD 6
#define N_PAR 6

int main(void) {
    int status = 0, h = 0, nspec = 0, ntfull = 0, nbands = 0, failures = 0;
    const char *names[N_PAR] = {"tau", "tage", "dust2", "gas_logu", "zred", "logzsol"};
    double tau[N_MOD] = {0.5, 1.0, 2.0, 5.0, 0.3, 1.5};
    double tage[N_MOD] = {1.0, 3.0, 8.0, 12.0, 0.5, 6.0};
    double dust2[N_MOD] = {0.0, 0.3, 0.6, 1.0, 0.1, 0.4};
    double logu[N_MOD] = {-3.5, -3.0, -2.5, -2.0, -1.5, -2.2};
    double zred[N_MOD] = {0.0, 0.1, 0.5, 1.0, 2.0, 0.3};
    double logz[N_MOD] = {-1.5, -1.0, -0.5, 0.0, 0.15, -0.8};
    const double *values[N_PAR] = {tau, tage, dust2, logu, zred, logz};

    fsps_context_create(&h, &status);
    fsps_context_setup(-1, "mist", "miles", "DL07", h, &status);
    if (status != 0) { printf("setup failed\n"); return 1; }
    fsps_context_set_int(h, "sfh", 4, &status);
    fsps_context_set_int(h, "add_neb_emission", 1, &status);
    fsps_context_set_int(h, "add_dust_emission", 1, &status);
    fsps_context_set_int(h, "add_igm_absorption", 1, &status);
    fsps_context_set_float(h, "dust1", 0.5, &status);
    fsps_context_get_nspec(h, &nspec, &status);
    fsps_context_get_ntfull(h, &ntfull, &status);
    fsps_context_get_nbands(h, &nbands, &status);

    double *mags = malloc(sizeof(double) * nbands * N_MOD), *spec = malloc(sizeof(double) * nspec * N_MOD);
    double *mags1 = malloc(sizeof(double) * nbands), *spec1 = malloc(sizeof(double) * nspec);
    double mass[N_MOD], lbol[N_MOD], mass1, lbol1;

    fsps_context_compute_csp_batch(h, N_MOD, N_PAR, names, values, NULL, mags, spec, mass, lbol,
                                   NULL, NULL, NULL, NULL, &status);
    printf("batch of %d models: status %d\n", N_MOD, status);
    if (status != 0) {
        char msg[256];
        int st2 = 0;
        fsps_get_last_error(&st2, msg, 256);
        printf("error: %s\n", msg);
        return 1;
    }

    double d_mag = 0.0, d_spec = 0.0, d_scal = 0.0;
    for (int m = 0; m < N_MOD; m++) {
        /* The same model, its parameters set on the context (logzsol through the arrays, so the
         * basis is the per-metallicity one in both cases) */
        fsps_context_set_float(h, "tau", tau[m], &status);
        fsps_context_set_float(h, "tage", tage[m], &status);
        fsps_context_set_float(h, "dust2", dust2[m], &status);
        fsps_context_set_float(h, "gas_logu", logu[m], &status);
        fsps_context_set_float(h, "zred", zred[m], &status);
        const char *zname[1] = {"logzsol"};
        const double *zval[1] = {&logz[m]};
        fsps_context_compute_csp_batch(h, 1, 1, zname, zval, NULL, mags1, spec1, &mass1, &lbol1,
                                       NULL, NULL, NULL, NULL, &status);
        if (status != 0) { printf("single model %d: status %d\n", m, status); return 1; }
        double peak = 0.0;
        for (int i = 0; i < nspec; i++) peak = fmax(peak, spec1[i]);
        for (int i = 0; i < nspec; i++)
            if (spec1[i] > 1e-6 * peak)
                d_spec = fmax(d_spec, fabs(spec[(size_t)m * nspec + i] / spec1[i] - 1.0));
        double mmin = 1e30;
        for (int b = 0; b < nbands; b++)
            if (mags1[b] == mags1[b]) mmin = fmin(mmin, mags1[b]);
        for (int b = 0; b < nbands; b++)
            if (mags1[b] == mags1[b] && mags1[b] < mmin + 30.0)
                d_mag = fmax(d_mag, fabs(mags[(size_t)m * nbands + b] - mags1[b]));
        d_scal = fmax(d_scal, fmax(fabs(mass[m] / mass1 - 1.0), fabs(lbol[m] - lbol1)));
    }
    printf("max difference, arrays vs one at a time (spectrum, mag, mass/lbol): %.2e %.2e %.2e\n",
           d_spec, d_mag, d_scal);
    if (!(d_spec <= 1e-5 && d_mag <= 1e-5 && d_scal <= 1e-12)) { printf("FAIL: parameter arrays\n"); failures++; }

    /* A parameter that changes the SSPs cannot vary per model; unknown names are rejected */
    const char *bad1[1] = {"imf1"};
    const char *bad2[1] = {"not_a_parameter"};
    const double *bv[1] = {tau};
    fsps_context_compute_csp_batch(h, N_MOD, 1, bad1, bv, NULL, mags, NULL, NULL, NULL, NULL, NULL, NULL, NULL,
                                   &status);
    if (status != 2) { printf("FAIL: imf1 per model gave status %d (expected 2)\n", status); failures++; }
    fsps_context_compute_csp_batch(h, N_MOD, 1, bad2, bv, NULL, mags, NULL, NULL, NULL, NULL, NULL, NULL, NULL,
                                   &status);
    if (status != 3) { printf("FAIL: unknown parameter gave status %d (expected 3)\n", status); failures++; }

    fsps_context_destroy(h, &status);
    free(mags); free(spec); free(mags1); free(spec1);
    printf(failures ? "FAILED\n" : "PASSED\n");
    return failures ? 1 : 0;
}
