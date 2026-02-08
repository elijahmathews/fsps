project('fsps', ['fortran', 'c'],
  version : '3.2.0',
  default_options : [
    'buildtype=release', 
    'warning_level=0',
    'default_library=shared'
  ]
)

# ===================================
# Compiler Configuration
# ===================================

fc = meson.get_compiler('fortran')
fc_id = fc.get_id()

# Universal arguments (supported by both gfortran and flang)
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

# ===================================
# Source Definitions
# ===================================

lib_sources = files(
    # ABI Modules
    'src/abi/fsps_c_driver.f90',

    # Core FSPS Modules
    'src/core/fsps_cache.f90',
    'src/core/fsps_constants.f90',
    'src/core/fsps_context_types.f90',
    'src/core/fsps_context.f90',
    'src/core/fsps_precision.f90',
    'src/core/fsps_types.f90',
    'src/core/sps_setup.f90',

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
    'src/spectra/spec_sbf.f90',
)

# ===================================
# Main Library Target
# ===================================

# This builds libfsps.so (and handles .mod generation automatically)
libfsps = library('fsps',
    lib_sources,
    version : meson.project_version(),
    soversion : '3',
    install : true
)

# Declare a dependency object for executables/tests to use
# This ensures they find the generated .mod files and link to the lib
fsps_dep = declare_dependency(
    link_with : libfsps,
    include_directories : include_directories('include')
)

# ===================================
# Program Executables
# ===================================

executable('autosps', 'src/programs/sps_program_autosps.f90', dependencies : fsps_dep, install : true)
executable('lesssimple', 'src/programs/sps_program_lesssimple.f90', dependencies : fsps_dep, install : true)
executable('simple', 'src/programs/sps_program_simple.f90', dependencies : fsps_dep, install : true)
executable('spec_bin', 'src/programs/sps_program_spec_bin.f90', dependencies : fsps_dep, install : true)

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
test('unit_tests', tests_units_exe)

# 3. Other Test Drivers
executable('test_runner', 'tests/test_runner.f90', dependencies : fsps_dep)
executable('test_cache', 'tests/test_cache.f90', dependencies : fsps_dep)
executable('generate_test_data', 'tests/generate_test_data.f90', dependencies : fsps_dep)

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