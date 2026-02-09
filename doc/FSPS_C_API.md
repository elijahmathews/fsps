# FSPS C Driver API

This document describes the C driver API exposed by FSPS. This API allows FSPS to be used as a shared library in C, C++, Python, Julia, Rust, and other languages supporting C interoperability.

## Build

From the repository root:
```sh
meson setup build
meson compile -C build
```

The C header is located at `include/fsps.h`.

## Architecture

The driver exposes two distinct APIs:

1.  **Context API (Recommended):** Uses opaque integer handles to manage independent, thread-safe FSPS instances.
2.  **Legacy API:** Manages a single global state (mirroring the original Fortran module behavior). Retained for backward compatibility.

## 1. Context API (Thread-Safe)

The Context API allows you to create multiple isolated instances of FSPS. This is the only way to use FSPS in a multi-threaded environment.

### Lifecycle

* **Create a context:**
    ```c
    void fsps_context_create(int *handle, int *status);
    ```
    Allocates a new FSPS instance and returns a positive integer `handle`.

* **Setup physics:**
    ```c
    void fsps_context_setup(int zin, const char *isoc, const char *spec, const char *dust, int handle, int *status);
    ```
    Loads the necessary data files for the context.
    * `zin`: Metallicity index (1 to `nz`). Use `-1` to load all metallicities (required for interpolation).
    * `isoc`, `spec`, `dust`: String identifiers for libraries (e.g., "MIST", "MILES", "DL07"). Pass empty strings to use defaults.

* **Destroy a context:**
    ```c
    void fsps_context_destroy(int handle, int *status);
    ```
    Frees memory associated with the handle.

### Configuration

Parameters are set using key-value pairs. The keys correspond to the FSPS parameter names (e.g., `tage`, `imf_type`, `logzsol`).

```c
void fsps_context_set_int(int handle, const char *key, int value, int *status);
void fsps_context_set_float(int handle, const char *key, double value, int *status);
void fsps_context_set_str(int handle, const char *key, const char *value, int *status);

```

### Computation

Currently, the Context API supports generating Simple Stellar Populations (SSPs).

```c
void fsps_context_compute_ssp(int handle, double *spec, double *mass, double *lbol, int *status);

```

* **Inputs:** `handle`
* **Outputs:**
* `spec`: Spectrum array. Dimensions: `[nspec, ntfull]`.
* `mass`: Surviving mass array. Dimensions: `[ntfull]`.
* `lbol`: Bolometric luminosity array. Dimensions: `[ntfull]`.


* **Note:** The caller is responsible for allocating the output buffers before calling this function. Use the Metadata functions below to determine the required sizes.

### Metadata & Dimensions

Use these functions to determine array sizes for allocation.

```c
// Get basic dimensions
void fsps_context_get_dims(int handle, int *n_spec, int *n_time, int *status);

// Individual getters
void fsps_context_get_nspec(int handle, int *n_spec, int *status);
void fsps_context_get_ntfull(int handle, int *n_time, int *status);
void fsps_context_get_nbands(int handle, int *n_bands, int *status);
void fsps_context_get_nindx(int handle, int *n_indices, int *status);
void fsps_context_get_nz(int handle, int *n_z, int *status);
void fsps_context_get_nemline(int handle, int *n_line, int *status);

// Get context-specific properties
void fsps_context_get_zsol(int handle, double *z_sol, int *status);
void fsps_context_get_paths(int handle, char *sps_home, int sps_len,
                            char *data_home, int data_len,
                            char *output_home, int out_len);

```

---

## 2. Legacy API (Global State)

This API manipulates the module-level global variables. It is **not** thread-safe.

### Initialization & Cleanup

```c
void fsps_initialize(int zin);
void fsps_finalize(void);

```

### Configuration

```c
void fsps_set_int(const char *key, int value);
void fsps_set_float(const char *key, double value);
void fsps_set_str(const char *key, const char *value);

// Advanced Setters
void fsps_set_sfh_tab(int ntab, double *age, double *sfr, double *met);
void fsps_set_mag_compute(int n_bands, int *mask);
void fsps_set_ssp_gen_age(int n_age, int *mask);
void fsps_set_ssp_lsf(int nsv, double *sigma, double wlo, double whi);

```

### Computation

Computes the model based on current global parameters (including `sfh` setting).

```c
void fsps_compute(double *spec_out);

```

### Data Access

These functions retrieve data calculated by the last call to `fsps_compute`.

```c
void fsps_get_mags(double zred, double *mags);
void fsps_get_stats(double *age, double *mass, double *lbol, double *sfr, 
                    double *mdust, double *mformed, double *emlines);
void fsps_get_spec(double *spec); // Retrieve internal buffer

```

---

## Error Handling

Most functions return a `status` integer.

* `0`: **FSPS_STATUS_OK**
* `>0`: Error or Warning code.

To retrieve the descriptive error message for the last failure:

```c
void fsps_get_last_error(int *status, char *message, int message_len);
void fsps_clear_error(void);

```

### Status Codes

| Range | Type | Description |
| --- | --- | --- |
| **0** | OK | Success |
| **100-199** | Error | Unknown parameter key |
| **200-299** | Warning | Parameter out of physical range (clamped or ignored) |
| **300-399** | Error | Initialization error or computation failure |
| **400-499** | Error | Lock/Threading error |

Enable debug printing to stdout:

```c
fsps_set_debug(1);

```

## Array Layout

All arrays are passed as flat C buffers. Internally, FSPS uses Fortran column-major order.

* **1D:** Continuous.
* **2D `[X, Y]`:** The `X` index varies fastest (contiguous in memory).
* **3D `[X, Y, Z]`:** `X` varies fastest, then `Y`.

**Common Dimensions:**

* `spec`: `[nspec, ntfull]`
* `mags`: `[nbands, ntfull]`
* `ssp_spec`: `[nspec, ntfull, nz]`

## Locking

The driver provides a cooperative lock. This is **not** a system mutex; it is a flag to help coordinate access if mixing Context and Legacy calls within a single process.

```c
void fsps_lock(int *status);
void fsps_unlock(int *status);

```

## Minimal Example (C)

```c
#include "fsps.h"
#include <stdlib.h>
#include <stdio.h>

int main() {
    int handle, status, nspec, nt;
    
    // 1. Create Context
    fsps_context_create(&handle, &status);
    
    // 2. Setup (Load all metallicities, default libraries)
    fsps_context_setup(-1, "", "", "", handle, &status);
    
    // 3. Get Dimensions to allocate buffers
    fsps_context_get_dims(handle, &nspec, &nt, &status);
    double *spec = (double*) malloc(nspec * nt * sizeof(double));
    double *mass = (double*) malloc(nt * sizeof(double));
    double *lbol = (double*) malloc(nt * sizeof(double));
    
    // 4. Set Parameters (e.g., Solar metallicity)
    fsps_context_set_int(handle, "zmet", 1, &status);
    
    // 5. Compute SSP
    fsps_context_compute_ssp(handle, spec, mass, lbol, &status);
    
    if (status != 0) {
        char msg[256];
        fsps_get_last_error(&status, msg, 256);
        printf("Error: %s\n", msg);
    }

    // 6. Cleanup
    fsps_context_destroy(handle, &status);
    free(spec); free(mass); free(lbol);
    
    return 0;
}

```