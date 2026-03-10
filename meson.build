project('fsps', ['fortran', 'c'],
  version : '3.2.0',
  license : 'MIT',
  license_files : ['LICENSE.md'],
  default_options : [
    'buildtype=release', 
    'warning_level=0',
    'default_library=shared'
  ],
  meson_version : '>=1.1',
)

# ===================================
# Compiler Configuration
# ===================================

fc = meson.get_compiler('fortran')
fc_id = fc.get_id()
host_cpu = host_machine.cpu_family()

# Build-time config header for CPP defines consumed by Fortran sources.
conf_data = configuration_data()

# Universal arguments
add_project_arguments('-cpp', language : 'fortran')
add_project_arguments('-I' + meson.current_build_dir(), language : 'fortran')

# GFortran-specific arguments
if fc_id == 'gcc'
    add_project_arguments('-ffree-line-length-132', language : 'fortran')
endif

# Installation paths
prefix = get_option('prefix')
libdir = get_option('libdir')
includedir = get_option('includedir')
datadir = get_option('datadir')

# -----------------------------------
# Hardware Acceleration (OpenACC)
# -----------------------------------

accel_target = get_option('accel')
accel_args = []
accel_defs = []
if accel_target != 'none'
    if fc_id != 'nvidia_hpc'
        warning('Hardware acceleration requested but compiler does not support OpenACC. Ignoring accel option.')
    else
        if accel_target == 'multicore'
            accel_args = ['-acc=multicore']
            accel_defs = ['-DACC_MULTICORE']
        elif accel_target == 'gpu'
            # Fat binary targeting Volta, Ampere, Hopper, plus PTX fallback
            accel_args = ['-acc=gpu', '-gpu=cc70,cc80,cc90,ptx']
            accel_defs = ['-DACC_GPU']
        endif
    endif
endif

# Create a dependency object that explicitly propagates these flags
# to both compile and link stages across all consuming targets.
accel_dep = declare_dependency(
    compile_args : accel_args + accel_defs,
    link_args : accel_args
)

# -----------------------------------
# CPU Parallelism (OpenMP)
# -----------------------------------

omp_dep = []
if get_option('enable_openmp')
    # Let Meson natively find OpenMP for both Fortran and C
    omp_dep = dependency('openmp', required: true)
endif

# -----------------------------------
# CPU Architecture Optimization
# -----------------------------------
cpu_baseline = get_option('cpu_baseline')
if get_option('buildtype') == 'release' and fc_id == 'gcc'
    if cpu_baseline == 'native'
        add_project_arguments('-march=native', language: 'fortran')
    elif cpu_baseline == 'x86-64-v3' and host_cpu == 'x86_64'
        add_project_arguments('-march=x86-64-v3', language: 'fortran')
    endif
endif

# -----------------------------------
# Data backend selection (legacy/hdf5/auto)
# -----------------------------------
data_backend_opt = get_option('data_backend')
hdf5_fortran_incdir = get_option('hdf5_fortran_incdir')
hdf5_fortran_link_args = get_option('hdf5_fortran_link_args')

hdf5_dep = disabler()
hdf5_backend_enabled = false

# Primary probe: Meson HDF5 Fortran discovery
hdf5_probe = dependency('hdf5', language: 'fortran', required: false)

if hdf5_probe.found()
    hdf5_dep = hdf5_probe
    hdf5_backend_enabled = true
elif hdf5_fortran_incdir != '' or hdf5_fortran_link_args.length() > 0
    # Manual fallback for systems where HDF5 Fortran modules are in non-standard paths
    hdf5_compile_args = []
    if hdf5_fortran_incdir != ''
        hdf5_compile_args += ['-I' + hdf5_fortran_incdir]
    endif

    hdf5_dep = declare_dependency(
        compile_args : hdf5_compile_args,
        link_args : hdf5_fortran_link_args
    )
    hdf5_backend_enabled = true
endif

if data_backend_opt == 'hdf5' and not hdf5_backend_enabled
    error('data_backend=hdf5 requested, but HDF5 Fortran dependency was not found. ' +
          'Provide -Dhdf5_fortran_incdir and/or -Dhdf5_fortran_link_args as fallback.')
endif

if data_backend_opt == 'auto'
    if hdf5_backend_enabled
        message('FSPS data backend: auto -> enabling HDF5 backend')
    else
        message('FSPS data backend: auto -> HDF5 unavailable, falling back to legacy backend')
    endif
endif

if data_backend_opt == 'legacy'
    hdf5_backend_enabled = false
endif

if hdf5_backend_enabled
    conf_data.set('FSPS_HAS_HDF5', 1)
else
    conf_data.set('FSPS_HAS_HDF5', 0)
endif

configure_file(
    output : 'fsps_build_config.h',
    configuration : conf_data
)

# ===================================
# Source Definitions
# ===================================

lib_sources = files(
    # Language Binding Modules
    'src/bindings/c/fsps_c_api.f90',

    # Core FSPS Modules
    'src/core/fsps_cache.f90',
    'src/core/fsps_constants.f90',
    'src/core/fsps_context_types.f90',
    'src/core/fsps_context.f90',
    'src/core/fsps_environment.f90',
    'src/core/fsps_initialization.f90',
    'src/core/fsps_strings.f90',
    'src/core/fsps_api.f90',
    'src/core/fsps_precision.f90',
    'src/core/fsps_types.f90',

    # IO Module
    'src/io/fsps_data_schema.f90',
    'src/io/fsps_data_backend.f90',
    'src/io/fsps_data_backend_legacy.f90',
    'src/io/fsps_data_registry.f90',
    'src/io/fsps_data_mapper.f90',
    'src/io/fsps_io.f90',

    # Math Modules
    'src/math/fsps_integration.f90',
    'src/math/fsps_interpolation.f90',
    'src/math/fsps_special_functions.f90',

    # Physics Modules
    'src/physics/fsps_cosmology.f90',
    'src/physics/fsps_csp.f90',
    'src/physics/fsps_dust.f90',
    'src/physics/fsps_gas.f90',
    'src/physics/fsps_imf.f90',
    'src/physics/fsps_sfh.f90',
    'src/physics/fsps_ssp.f90',
    'src/physics/fsps_stellar_modifications.f90',

    # Spectra Modules
    'src/spectra/fsps_photometry.f90',
    'src/spectra/fsps_smoothing.f90',
    'src/spectra/fsps_spectral_indices.f90',
    'src/spectra/fsps_spectral_library.f90',
)

if hdf5_backend_enabled
    lib_sources += files('src/io/fsps_data_backend_hdf5.f90')
endif

backend_deps = []
if hdf5_backend_enabled
    backend_deps += hdf5_dep
endif

# ===================================
# Main Library Target
# ===================================

# This builds libfsps.so (and handles .mod generation automatically)
libfsps = library('fsps',
    lib_sources,
    version : meson.project_version(),
    soversion : '3',
    install : true,
    dependencies : [omp_dep, accel_dep] + backend_deps,
    include_directories : include_directories('.'),
)

# Declare a dependency object for executables/tests to use
# This ensures they find the generated .mod files and link to the lib
fsps_dep = declare_dependency(
    link_with : libfsps,
    include_directories : include_directories('include'),
    dependencies : [omp_dep, accel_dep] + backend_deps
)

# ===================================
# Tests
# ===================================

# 1. Test Utilities (Helper for other tests)
test_utils_lib = static_library(
    'test_utils', 
    'tests/test_utils.f90', 
    dependencies : fsps_dep
)

test_utils_dep = declare_dependency(link_with : test_utils_lib)

# 2. Unit Tests
unit_test_sources = files(
    'tests/core/test_fsps_cache.f90',
    'tests/core/test_fsps_context.f90',
    'tests/core/test_fsps_context_types.f90',
    'tests/core/test_fsps_strings.f90',
    'tests/core/test_fsps_types.f90',

    'tests/io/test_fsps_io.f90',

    'tests/math/test_fsps_integration.f90',
    'tests/math/test_fsps_interpolation.f90',
    'tests/math/test_fsps_special_functions.f90',

    'tests/physics/test_fsps_cosmology.f90',
    'tests/physics/test_fsps_csp.f90',
    'tests/physics/test_fsps_dust.f90',
    'tests/physics/test_fsps_gas.f90',
    'tests/physics/test_fsps_imf.f90',
    'tests/physics/test_fsps_sfh.f90',
    'tests/physics/test_fsps_ssp.f90',
    'tests/physics/test_fsps_stellar_modifications.f90',

    'tests/spectra/test_fsps_photometry.f90',
    'tests/spectra/test_fsps_smoothing.f90',
    'tests/spectra/test_fsps_spectral_indices.f90',
    'tests/spectra/test_fsps_spectral_library.f90',
    
    'tests/test_fsps.f90' # Main driver
)

tests_units_exe = executable(
    'tests_units', 
    unit_test_sources, 
    dependencies : [fsps_dep, test_utils_dep]
)
test('unit_tests', tests_units_exe, suite : 'unit')

# 3. Other Test Drivers
test_runner_exe = executable('test_runner', 'tests/test_runner.f90', dependencies : fsps_dep)
executable('test_cache', 'tests/test_cache.f90', dependencies : fsps_dep)
executable('generate_test_data', 'tests/generate_test_data.f90', dependencies : fsps_dep)
executable('benchmark_hotpaths', 'tests/benchmark_hotpaths.f90', dependencies : fsps_dep)
executable('test_backend_parity', 'tests/io/test_backend_parity.f90', dependencies : fsps_dep)

if get_option('enable_stress_testing')
    stress_test_exe = executable('stress_test', 'tests/stress_test.f90', dependencies : fsps_dep)
    stress_env = environment()
    stress_env.set('FSPS_DATA_HOME', meson.project_source_root())
    test('stress_test', stress_test_exe,
        suite : 'stress',
        env : stress_env,
        timeout : 1200
    )
endif

# 3a. Regression Tests (run through Meson test framework)
regression_env = environment()
regression_env.set('FSPS_TEST_RTOL', '1.0E-5')
regression_env.set('FSPS_DATA_HOME', meson.project_source_root())

test('regression_miles_mist',
        test_runner_exe,
        args : [
            join_paths(meson.project_source_root(), 'tests/data/sps_ref_MILES-1_MIST-1.bin'),
            '--isoc', 'mist', '--spec', 'miles'
        ],
        env : regression_env,
        suite : 'regression',
        timeout : 300
)

test('regression_basel_padova',
        test_runner_exe,
        args : [
            join_paths(meson.project_source_root(), 'tests/data/sps_ref_MILES-0_MIST-0_BASEL-1_PADOVA-1.bin'),
            '--isoc', 'pdva', '--spec', 'basel'
        ],
        env : regression_env,
        suite : 'regression',
        timeout : 300
)

test('regression_themis_dust',
        test_runner_exe,
        args : [
            join_paths(meson.project_source_root(), 'tests/data/sps_ref_MILES-1_MIST-1_THEMIS-1_DL07-0.bin'),
            '--isoc', 'mist', '--spec', 'miles', '--dust', 'themis'
        ],
        env : regression_env,
        suite : 'regression',
        timeout : 300
)

test('regression_c3k_mist',
        test_runner_exe,
        args : [
            join_paths(meson.project_source_root(), 'tests/data/sps_ref_MILES-0_C3K-1_MIST-1.bin'),
            '--isoc', 'mist', '--spec', 'c3k_afe+0.0'
        ],
        env : regression_env,
        suite : 'regression',
        timeout : 300
)

test('regression_bpass',
        test_runner_exe,
        args : [
            join_paths(meson.project_source_root(), 'tests/data/sps_ref_MIST-0_BPASS-1.bin'),
            '--isoc', 'bpss'
        ],
        env : regression_env,
        suite : 'regression',
        timeout : 300
)

# 4. C Tests
executable('test_c_driver', 'tests/test_c_driver.c', 
    dependencies : fsps_dep,
    include_directories : include_directories('include')
)

# ===================================
# Installation (PkgConfig & Data)
# ===================================

# Install Headers
install_headers('include/fsps.h')

# Install Data (recursively copies the directory)
install_subdir('data', install_dir : join_paths(datadir, 'fsps'))

# Pkg-config generation
pkg = import('pkgconfig')
pkg.generate(libfsps,
    name : 'fsps',
    description : 'Flexible Stellar Population Synthesis',
    url : 'https://github.com/cconroy20/fsps'
)
