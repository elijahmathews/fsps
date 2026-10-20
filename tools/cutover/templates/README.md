# libopensps

libopensps is a library for stellar population synthesis: spectra, magnitudes, spectral indices
and line luminosities of simple and composite stellar populations, with dust attenuation and
emission, nebular emission, AGN dust emission and the intergalactic medium. It is a modernized
version of **FSPS** (Flexible Stellar Population Synthesis) by Charlie Conroy and contributors,
written in Fortran 2008 and built with Meson, with:

- a stable, versioned C API (`opensps.h`, `libopensps`, pkg-config `opensps`) for use from C,
  C++, Python, Julia, Rust and other languages;
- no global state: every model lives in a context, and contexts can run in parallel threads;
- batched composite populations for fitting (many models per call, one array per parameter);
- optional FP32 kernels and GPU offload (OpenACC, OpenMP target).

## Credit: FSPS

The physics of libopensps is that of FSPS, and FSPS is the work of **Charlie Conroy** and its
contributors, listed in [AUTHORS](AUTHORS). libopensps started from a snapshot of the modernized
fork [elijahmathews/fsps](https://github.com/elijahmathews/fsps) (tag `fsps-final`, which keeps
the full history); the original FSPS repository is
[cconroy20/fsps](https://github.com/cconroy20/fsps). libopensps is not an official FSPS release.

**If you use libopensps in published work, please cite the FSPS papers:**

- Conroy, Gunn & White 2009, ApJ, 699, 486
- Conroy & Gunn 2010, ApJ, 712, 833

and libopensps itself ([CITATION.cff](CITATION.cff); BibTeX entries of the FSPS papers are in
[CITATION.bib](CITATION.bib)). The stellar libraries, isochrones and other data you select have
their own references; cite those too.

## Build

Requirements: a Fortran 2008 compiler (gfortran, LLVM flang, Intel ifx or nvfortran), a C
compiler, Meson (>= 1.1) and Ninja.

```sh
meson setup build --buildtype=release                               # FP64
meson setup build32 --buildtype=release -Dkernel_precision=single   # FP32 kernels
meson compile -C build
meson install -C build        # library, opensps.h, opensps.pc
```

Options (`meson.options`): `kernel_precision` (`double` or `single`), `enable_openmp`,
`accel` (`none`, `multicore`, `gpu`, `omp_gpu`), `gpu_arch`, `cublas`, `cpu_baseline`. GPU
builds are described in `doc/GPU_PLAN.md`.

## Data

The data (isochrones, stellar spectra, dust, nebular and filter libraries; about 0.7 GB) are not
in this repository. They are distributed as versioned archives in a documented format
(`doc/DATA_FORMAT.md`) and fetched with checksums verified:

```sh
python3 tools/fetch_data.py ~/opensps-data             # TODO: the data-v1 release manifest URL
export FSPS_DATA_HOME=~/opensps-data
```

The library looks for its data in `SPS_HOME`, then `FSPS_DATA_HOME`, then
`$XDG_DATA_HOME/fsps`, `~/.local/share/fsps`, `/usr/share/fsps` and `/usr/local/share/fsps`
(these names are unchanged from FSPS, so existing setups keep working). Any library in the
standard format can be added without code changes (`doc/DATA_FORMAT.md`).

## Test

```sh
meson test -C build --suite unit                         # no external data needed
FSPS_DATA_HOME=~/opensps-data meson test -C build --suite regression
```

The unit suite runs on the small test library in `tests/data/fsps-test-data`. The regression
suite (references, a parameterized case matrix, physical identities, the batched against the
single CSP) needs the full data. Both pass in FP64 and FP32 builds before a change is merged
([CONTRIBUTING.md](CONTRIBUTING.md)).

## The C API

```c
#include <stdio.h>
#include "opensps.h"

int main(void) {
    int h, status, n_out;
    opensps_context_create(&h, &status);
    /* all metallicities (zin = -1); MIST isochrones, MILES spectra, DL07 dust emission */
    opensps_context_setup(-1, "mist", "miles", "dl07", h, &status);
    opensps_context_set_int(h, "sfh", 4, &status);           /* delayed tau */
    opensps_context_set_float(h, "tau", 2.0, &status);
    opensps_context_set_float(h, "tage", 10.0, &status);
    opensps_context_set_float(h, "logzsol", -0.3, &status);
    opensps_context_set_float(h, "dust2", 0.4, &status);
    opensps_context_compute_csp(h, 1, &n_out, &status);       /* interpolated in metallicity */
    if (status != OPENSPS_STATUS_OK) {
        char msg[256];
        opensps_get_last_error(&status, msg, (int)sizeof msg);
        fprintf(stderr, "error %d: %s\n", status, msg);
        return 1;
    }
    double mass, lbol;
    opensps_context_get_csp(h, NULL, &mass, &lbol, NULL, NULL, NULL, NULL, NULL, NULL, NULL,
                            &status);
    printf("surviving mass %.4f Msun, log L_bol %.4f\n", mass, lbol);
    opensps_context_destroy(h, &status);
    return 0;
}
```

```sh
cc example.c $(pkg-config --cflags --libs opensps)
```

Every function reports through a `status` argument (`OPENSPS_STATUS_OK` = 0; the codes are in
`opensps.h`), and `opensps_get_last_error` returns the message. The batched entry point
`opensps_context_compute_csp_batch` evaluates many models in one call. The API is versioned
(`OPENSPS_API_VERSION_MAJOR/MINOR/PATCH`, `opensps_get_driver_version`); the reference is
`doc/C_API.md`.

Coming from elijahmathews/fsps: the C names changed from `fsps_*` / `FSPS_*` to `opensps_*` /
`OPENSPS_*`, the header from `fsps.h` to `opensps.h`, the library from `libfsps` to
`libopensps` and the pkg-config name from `fsps` to `opensps`. A deprecated `fsps.h` that maps
the old names is installed for now. Environment variables, parameter names and the data
directory names are unchanged.

## Language bindings

Bindings live in their own repositories and depend on a released version of this library:

- Python: [opensps-python](https://github.com/elijahmathews/opensps-python)
  (`pip install opensps`, `import opensps`), including a compatibility layer for code written
  for python-fsps.

## Contributing and license

See [CONTRIBUTING.md](CONTRIBUTING.md). libopensps is released under the MIT license, with the
FSPS copyright notice ([LICENSE.md](LICENSE.md)).
