/* The single-CSP entry points of the C API (fsps_context_compute_csp / fsps_context_get_csp),
 * the context grid getters and the tabular SFH, as language bindings use them:
 *  - the API version macros equal fsps_get_driver_version;
 *  - grids: sizes, increasing wavelengths and ages, filter data;
 *  - an SSP (sfh = 0, tage = 0) gives the fsps_context_compute_ssp mass and L_bol exactly;
 *  - a tau model at tage > 0 equals the batch of one with the same parameters (to the batch
 *    test's tolerance), at zmet (zcontinuous = 0) and at logzsol (zcontinuous = 1);
 *  - the stored SSP grids are invalidated by an SSP parameter and reused otherwise;
 *  - a tabular SFH (sfh = 3) of constant SFR reproduces sfh = 1 with const = 1;
 *  - error statuses (no CSP yet, no magnitudes in fast mode, bad zcontinuous, sfh = 3 without a table, bad handle).
 * Requires the FSPS data tree (SPS_HOME), like the regression tests. */
#include "fsps.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

static int failures = 0;

static void check(int ok, const char *what) {
    printf("%s: %s\n", ok ? "ok" : "FAIL", what);
    if (!ok) failures++;
}

/* Max relative difference over pixels above 1e-6 of the peak */
static double spec_diff(const double *a, const double *b, int n) {
    double peak = 0.0, d = 0.0;
    for (int i = 0; i < n; i++) peak = fmax(peak, fabs(b[i]));
    for (int i = 0; i < n; i++)
        if (fabs(b[i]) > 1e-6 * peak) d = fmax(d, fabs(a[i] / b[i] - 1.0));
    return d;
}

int main(void) {
    int status = 0, h = 0, nspec = 0, ntfull = 0, nbands = 0, nz = 0, nline = 0, n_out = 0;
    int major = 0, minor = 0, patch = 0;

    fsps_get_driver_version(&major, &minor, &patch);
    check(major == FSPS_API_VERSION_MAJOR && minor == FSPS_API_VERSION_MINOR &&
              patch == FSPS_API_VERSION_PATCH, "API version macros = fsps_get_driver_version");

    fsps_context_create(&h, &status);
    fsps_context_setup(-1, "mist", "miles", "DL07", h, &status);
    if (status != 0) { printf("setup failed\n"); return 1; }
    fsps_context_get_nspec(h, &nspec, &status);
    fsps_context_get_ntfull(h, &ntfull, &status);
    fsps_context_get_nbands(h, &nbands, &status);
    fsps_context_get_nz(h, &nz, &status);
    fsps_context_get_nemline(h, &nline, &status);

    /* Grids */
    double *lam = malloc(sizeof(double) * nspec), *ages = malloc(sizeof(double) * ntfull);
    double *zl = malloc(sizeof(double) * nz), *eml = malloc(sizeof(double) * nline);
    double *leff = malloc(sizeof(double) * nbands), *mvega = malloc(sizeof(double) * nbands);
    double *msun = malloc(sizeof(double) * nbands);
    int st[5];
    fsps_context_get_lambda(h, lam, &st[0]);
    fsps_context_get_ages(h, ages, &st[1]);
    fsps_context_get_zlegend(h, zl, &st[2]);
    fsps_context_get_emlambda(h, eml, &st[3]);
    fsps_context_get_filter_data(h, leff, mvega, msun, &st[4]);
    int inc = 1;
    for (int i = 1; i < nspec; i++) inc &= lam[i] > lam[i - 1];
    for (int i = 1; i < ntfull; i++) inc &= ages[i] > ages[i - 1];
    for (int i = 1; i < nz; i++) inc &= zl[i] > zl[i - 1];
    check(st[0] == 0 && st[1] == 0 && st[2] == 0 && st[3] == 0 && st[4] == 0, "grid getters status 0");
    check(inc && lam[0] > 0.0 && ages[0] >= 5.0 && ages[ntfull - 1] < 10.5, "grids increasing, ages in log yr");
    check(mvega[0] == 0.0 && leff[0] > 0.0 && eml[0] > 0.0, "filter data and line wavelengths");
    fsps_context_get_lambda(-3, lam, &status);
    check(status == 1, "bad handle gives status 1");

    /* Before any CSP: status 302 */
    double age1, mass1, lbol1, sfr1, mdust1, mform1;
    fsps_context_get_csp(h, &age1, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, &status);
    check(status == 302, "get_csp before compute: 302");
    fsps_context_compute_csp(h, 7, &n_out, &status);
    check(status == 306, "zcontinuous 7: 306");

    /* SSP: all ages, mass and L_bol of compute_ssp exactly */
    fsps_context_set_int(h, "zmet", 10, &status);
    fsps_context_set_int(h, "sfh", 0, &status);
    fsps_context_set_float(h, "tage", 0.0, &status);
    double *sspec = malloc(sizeof(double) * nspec * ntfull), *smass = malloc(sizeof(double) * ntfull);
    double *slbol = malloc(sizeof(double) * ntfull);
    fsps_context_compute_ssp(h, sspec, smass, slbol, &status);
    fsps_context_compute_csp(h, 0, &n_out, &status);
    check(status == 0 && n_out == ntfull, "SSP: n_out = ntfull");
    double *cmass = malloc(sizeof(double) * ntfull), *clbol = malloc(sizeof(double) * ntfull);
    double *cage = malloc(sizeof(double) * ntfull);
    fsps_context_get_csp(h, cage, cmass, clbol, NULL, NULL, NULL, NULL, NULL, NULL, NULL, &status);
    int same = status == 0;
    for (int i = 0; i < ntfull; i++) {
        if (!(cmass[i] == smass[i] && fabs(clbol[i] - slbol[i]) < 1e-13 && cage[i] == ages[i]))
            printf("# age %d: %.17g %.17g, mass %.17g %.17g, lbol %.17g %.17g\n", i, cage[i], ages[i], cmass[i],
                   smass[i], clbol[i], slbol[i]);
        same &= cmass[i] == smass[i] && fabs(clbol[i] - slbol[i]) < 1e-13 && cage[i] == ages[i];
    }
    check(same, "SSP mass and ages = compute_ssp exactly, L_bol to 1e-13 dex");

    /* A tau model at tage = 5 Gyr with dust and nebular emission vs the batch of one */
    fsps_context_set_int(h, "sfh", 1, &status);
    fsps_context_set_float(h, "tau", 2.0, &status);
    fsps_context_set_float(h, "tage", 5.0, &status);
    fsps_context_set_float(h, "dust2", 0.4, &status);
    fsps_context_set_int(h, "dust_type", 2, &status);
    fsps_context_set_int(h, "add_neb_emission", 1, &status);
    double *spec = malloc(sizeof(double) * nspec), *mags = malloc(sizeof(double) * nbands);
    double *bspec = malloc(sizeof(double) * nspec), *bmags = malloc(sizeof(double) * nbands);
    double *lines = malloc(sizeof(double) * nline), *blines = malloc(sizeof(double) * nline);
    double bmass, blbol;
    for (int zc = 0; zc <= 1; zc++) {
        if (zc == 1) fsps_context_set_float(h, "logzsol", -0.45, &status);
        fsps_context_compute_csp(h, zc, &n_out, &status);
        check(status == 0 && n_out == 1, "tau model: one output");
        fsps_context_get_csp(h, &age1, &mass1, &lbol1, &sfr1, &mdust1, &mform1, spec, mags, lines, NULL, &status);
        check(status == 0 && fabs(age1 - log10(5e9)) < 1e-12 && mform1 > 0.99 && sfr1 > 0.0, "tau model outputs");
        const char *zname[1] = {"logzsol"};
        double lz = -0.45;
        const double *zval[1] = {&lz};
        fsps_context_compute_csp_batch(h, 1, zc, zname, zval, NULL, bmags, bspec, &bmass, &blbol, NULL, NULL,
                                       blines, NULL, &status);
        double dm = 0.0;
        for (int b = 0; b < nbands; b++)
            if (bmags[b] < 60.0) dm = fmax(dm, fabs(mags[b] - bmags[b]));
        double ds = spec_diff(bspec, spec, nspec), dl = spec_diff(blines, lines, nline);
        printf("# zcontinuous %d: single vs batch of one: spec %.2e, mag %.2e, lines %.2e, mass %.2e\n",
               zc, ds, dm, dl, fabs(bmass / mass1 - 1.0));
        check(status == 0 && ds < 1e-5 && dm < 1e-5 && dl < 1e-5 && fabs(bmass / mass1 - 1.0) < 1e-12 &&
                  fabs(blbol - lbol1) < 1e-12, "single CSP = batch of one");
    }

    /* Stored SSPs: invalidated by an SSP parameter, reused (identical) when it is restored */
    double *spec2 = malloc(sizeof(double) * nspec);
    fsps_context_set_float(h, "imf1", 1.5, &status);
    fsps_context_compute_csp(h, 1, &n_out, &status);
    fsps_context_get_csp(h, NULL, NULL, NULL, NULL, NULL, NULL, spec2, NULL, NULL, NULL, &status);
    check(spec_diff(spec2, spec, nspec) > 1e-4, "imf1 changes the CSP");
    fsps_context_set_float(h, "imf1", 1.3, &status);
    fsps_context_compute_csp(h, 1, &n_out, &status);
    fsps_context_get_csp(h, NULL, NULL, NULL, NULL, NULL, NULL, spec2, NULL, NULL, NULL, &status);
    check(spec_diff(spec2, spec, nspec) == 0.0, "imf1 restored: identical CSP");

    /* sfh = 3 without a table: error; constant SFR table = sfh 1 with const = 1 */
    fsps_context_set_int(h, "sfh", 3, &status);
    fsps_context_compute_csp(h, 1, &n_out, &status);
    check(status == 311, "sfh 3 without a table: 311");
    double tab_age[3] = {0.0, 6.0, 13.8}, tab_sfr[3] = {1.0, 1.0, 1.0};
    fsps_context_set_sfh_tab(h, 3, tab_age, tab_sfr, NULL, &status);
    fsps_context_compute_csp(h, 1, &n_out, &status);
    fsps_context_get_csp(h, NULL, &mass1, &lbol1, NULL, NULL, &mform1, spec2, NULL, NULL, NULL, &status);
    printf("# sfh 3: mass formed %.6g, mass %.6g, log L_bol %.6g\n", mform1, mass1, lbol1);
    check(status == 0, "sfh 3 with a table");
    /* The tabular SFH is in M_sun/yr, not normalized (as in FSPS): 5e9 M_sun formed by 5 Gyr.
     * (compute_csp_scenario reports mformed = 1 here; see the PR.) Per solar mass formed: */
    mass1 /= 5e9;
    for (int i = 0; i < nspec; i++) spec2[i] /= 5e9;
    fsps_context_set_int(h, "sfh", 1, &status);
    fsps_context_set_float(h, "const", 1.0, &status);
    fsps_context_compute_csp(h, 1, &n_out, &status);
    fsps_context_get_csp(h, NULL, &bmass, &blbol, NULL, NULL, NULL, spec, NULL, NULL, NULL, &status);
    printf("# tabular constant SFR vs const = 1: spec %.2e, mass %.2e\n", spec_diff(spec2, spec, nspec),
           fabs(bmass / mass1 - 1.0));
    check(spec_diff(spec2, spec, nspec) < 1e-3 && fabs(bmass / mass1 - 1.0) < 1e-3, "tabular = const SFR");
    fsps_context_set_sfh_tab(h, 0, tab_age, tab_sfr, NULL, &status);
    check(status == 3, "ntab 0: status 3");

    /* fast mode: no magnitudes */
    fsps_context_set_fast_mode(h, 1, &status);
    fsps_context_compute_csp(h, 1, &n_out, &status);
    fsps_context_get_csp(h, NULL, NULL, NULL, NULL, NULL, NULL, NULL, mags, NULL, NULL, &status);
    check(status == 309, "mags in fast mode: 309");

    fsps_context_destroy(h, &status);
    printf(failures ? "FAILED\n" : "PASSED\n");
    return failures ? 1 : 0;
}
