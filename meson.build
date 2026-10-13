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

# Universal arguments
add_project_arguments('-cpp', language : 'fortran')

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

# Native threads: the setup-cache lock (src/core/fsps_setup_lock.c) is a pthread mutex
threads_dep = dependency('threads')

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

# ===================================
# Source Definitions
# ===================================

lib_sources = files(
    # Language Binding Modules
    'src/bindings/c/fsps_c_api.f90',

    # Core FSPS Modules
    'src/core/fsps_cache.f90',
    'src/core/fsps_setup_lock.c',
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

# ===================================
# Main Library Target
# ===================================

# This builds libfsps.so (and handles .mod generation automatically)
libfsps = library('fsps',
    lib_sources,
    version : meson.project_version(),
    soversion : '3',
    install : true,
    dependencies : [omp_dep, accel_dep, threads_dep],
    # The library has one C source (setup-cache lock); link with the Fortran compiler so its
    # runtime (libgfortran, flang runtime, ...) is always included
    link_language : 'fortran',
)

# Declare a dependency object for executables/tests to use
# This ensures they find the generated .mod files and link to the lib
fsps_dep = declare_dependency(
    link_with : libfsps,
    include_directories : include_directories('include'),
    dependencies : [omp_dep, accel_dep, threads_dep]
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
# One TAP test per module: meson runs them in parallel and reports every assertion
unit_env = environment()
unit_env.set('FSPS_TEST_TAP', '1')
unit_modules = [
    'fsps_cache', 'fsps_context', 'fsps_context_types', 'fsps_strings', 'fsps_types', 'fsps_io',
    'fsps_integration', 'fsps_interpolation', 'fsps_special_functions', 'fsps_cosmology', 'fsps_csp',
    'fsps_dust', 'fsps_gas', 'fsps_imf', 'fsps_sfh', 'fsps_ssp', 'fsps_stellar_modifications',
    'fsps_photometry', 'fsps_smoothing', 'fsps_spectral_indices', 'fsps_spectral_library',
]
foreach unit_module : unit_modules
    test('unit_' + unit_module.substring(5), tests_units_exe,
        args : [unit_module],
        env : unit_env,
        protocol : 'tap',
        suite : 'unit'
    )
endforeach

# 3. Other Test Drivers
test_runner_exe = executable('test_runner', 'tests/test_runner.f90', dependencies : fsps_dep)
executable('test_cache', 'tests/test_cache.f90', dependencies : fsps_dep)
executable('generate_test_data', 'tests/generate_test_data.f90', dependencies : fsps_dep)
executable('generate_opensps_fixtures', 'tests/generate_opensps_fixtures.f90', dependencies : fsps_dep)
executable('benchmark_hotpaths', 'tests/benchmark_hotpaths.f90', dependencies : fsps_dep)

if get_option('enable_stress_testing')
    stress_test_exe = executable('stress_test', 'tests/stress_test.f90', dependencies : fsps_dep)
    stress_env = environment()
    stress_env.set('FSPS_DATA_HOME', meson.project_source_root())
    test('stress_test', stress_test_exe,
        suite : 'stress',
        env : stress_env,
        # Phase D serializes repeated full setups (~160 s each under ThreadSanitizer)
        timeout : 2400
    )
endif

# 3a. Regression Tests (run through Meson test framework)
# The C3K library (~600 MB) may be left out of a checkout (e.g. a sparse CI clone); its tests are
# then skipped rather than failing on missing files.
fs = import('fs')
c3k_available = fs.is_dir('data/spectra/C3K')
if not c3k_available
    message('data/spectra/C3K not found: C3K regression tests are skipped')
endif
regression_env = environment()
regression_env.set('FSPS_TEST_RTOL', '1.0E-5')
regression_env.set('FSPS_DATA_HOME', meson.project_source_root())
# Programs built on tests/test_utils.f90 report each assertion as TAP
regression_env.set('FSPS_TEST_TAP', '1')

# Integration test: contexts sharing a setup cache must not leak state into each other
context_isolation_exe = executable('test_context_isolation', 'tests/test_context_isolation.f90',
    dependencies : [fsps_dep, test_utils_dep])
test('context_isolation', context_isolation_exe,
    protocol : 'tap',
    env : regression_env,
    suite : 'regression',
    timeout : 600
)

# Integration test: outputs checked against independent physical identities and published values
physics_oracles_exe = executable('test_physics_oracles', 'tests/test_physics_oracles.f90',
    dependencies : [fsps_dep, test_utils_dep])
test('physics_oracles', physics_oracles_exe,
    protocol : 'tap',
    env : regression_env,
    suite : 'regression',
    timeout : 900
)

# Parameterized regression over the parity case matrix (shared with generate_opensps_fixtures.sh),
# sharded by case family so that meson runs the shards in parallel
param_regression_exe = executable('test_param_regression', 'tests/test_param_regression.f90',
    dependencies : [fsps_dep, test_utils_dep])
param_filter_args = c3k_available ? [] : ['--skip-spec', 'c3k']
foreach shard : ['ssp_', 'csp_', 'lib_']
    test('param_regression_' + shard.strip('_'), param_regression_exe,
        protocol : 'tap',
        args : [
            join_paths(meson.project_source_root(), 'tests/data/param_cases.psv'),
            join_paths(meson.project_source_root(), 'tests/data/param_ref'),
            '--only', shard
        ] + param_filter_args,
        env : regression_env,
        suite : 'regression',
        timeout : 1200
    )
endforeach

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

if c3k_available
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
endif

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
