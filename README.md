FSPS: Flexible Stellar Population Synthesis
=====
![Version Badge](https://img.shields.io/badge/version-v3.2-blue) ![GitHub Actions Workflow Status](https://img.shields.io/github/actions/workflow/status/elijahmathews/fsps/test.yml)
 [![codecov](https://img.shields.io/codecov/c/github/elijahmathews/fsps)](https://codecov.io/github/elijahmathews/fsps)

> [!WARNING]
> This is a major refactor of the standard FSPS codebase. It introduces breaking API changes to support thread safety and C-interoperability. For the stable, single-threaded legacy version, please visit the [original repository](https://github.com/cconroy20/fsps).

## Overview

This repository contains a modernized, thread-safe implementation of the **Flexible Stellar Population Synthesis (FSPS)** code. While maintaining scientific accuracy with the original Fortran implementation, this fork overhauls the internal architecture to support:

1. **No Global State:** All simulation state is encapsulated in `fsps_context_t` handles.
2. **Thread Safety:** Multiple FSPS contexts can run simultaneously on different threads (OpenMP/pthreads compatible).
3. **Stable C ABI:** A standardized C driver (`fsps_c_driver`) allows direct linking from C, C++, Rust, Julia, and Python (via `ctypes`/`cffi`) without relying on `f2py`.
4. **Modern Build System:** Streamlined Meson build system supporting shared libraries, `pkg-config`, and standard installation paths.

## Installation

### Prerequisites

* A modern Fortran compiler (`gfortran` or `flang`).
* The FSPS data files (see below).

### Clone and Build

```sh
git clone https://github.com/elijahmathews/fsps
cd fsps

# Setup Meson
meson setup build

# Compile FSPS
meson compile -C build

# Run unit tests
meson test -v -C build
```

## C Driver API (New)

This refactor exposes the entire FSPS physical model through a clean C interface. This is the recommended way to interface with FSPS for external applications.

### Basic Usage Example (C)

```c
#include <stdio.h>
#include "fsps.h"

int main() {
    int status;
    int handle;
    
    // 1. Create a context (thread-safe independent instance)
    fsps_context_create(&handle, &status);
    
    // 2. Setup physics (MIST isochrones, MILES spectra)
    fsps_context_setup(handle, "MIST", "MILES", "DL07", &status);
    
    // 3. Configure parameters (e.g., Solar Metallicity, 1 Gyr old)
    fsps_set_param_float(handle, "logzsol", 0.0, &status);
    fsps_set_param_float(handle, "tage", 1.0, &status);
    
    // 4. Compute Spectrum
    // (Buffers would be allocated here, pointers passed to compute)
    // fsps_compute(handle, ...);
    
    // 5. Cleanup
    fsps_context_destroy(handle, &status);
    return 0;
}

```

Detailed documentation on the C API, array layouts, and error codes can be found in [`doc/FSPS_C_API.md`](doc/FSPS_C_API.md).

## Directory Structure

* **`src/`**: Modern Fortran source code. Organized by subsystem (physics, spectra, cosmology, C-driver).
* **`include/`**: C header files (`fsps.h`).
* **`data/`**: Location for physics tables (isochrones, dust models, etc).
* **`build/`**: Build artifacts, object files, and compiled modules.
* **`tests/`**: Regression suite comparing new C-driver outputs against legacy Fortran baselines.

## Citations and Acknowledgments

If you use this code in your research, you **must** cite the original FSPS papers describing the physical models:

* Conroy, Gunn, & White (2009, ApJ, 699, 486)
* Conroy & Gunn (2010, ApJ, 712, 833)

Please also refer to [`CITATION.bib`](CITATION.bib) for specific citations regarding the MIST isochrones and MILES spectral libraries if you utilize the default settings.

## License

This code is released under the MIT License. See [license file](LICENSE.md) for details.
This refactor is based on the original work by Charlie Conroy and contributors.