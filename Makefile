# ===================================
# Compiler Configuration
# ===================================

# Default to gfortran, allow for an override (e.g. make FC=ifort)
FC = gfortran

# Default flags
FFLAGS ?= -O3 -cpp -fPIC

# Directory configuration
SRC_DIR := src
TEST_DIR := tests
BUILD_DIR := build
SRC_SUBDIRS := core programs spectra sfh physics cosmology math io imf abi
SRC_DIRS := $(addprefix $(SRC_DIR)/, $(SRC_SUBDIRS))

# Module directory flags:
# Gfortran uses -J to specify where to put/find .mod files
# Intel (ifort) uses -module. We default to Gfortran syntax here.
# You can override this: make MOD_FLAG=-module
MOD_FLAG ?= -J

# Combine flags to include the build dir for .mod search
FCFLAGS := $(FFLAGS) $(MOD_FLAG)$(BUILD_DIR) -I$(BUILD_DIR)

# ===================================
# Source and Object Definitions
# ===================================

# Tell make to look for source files in these directories
VPATH = $(SRC_DIRS):$(TEST_DIR)

# The list of programs to build (executables)
PROGS = simple lesssimple autosps spec_bin

# The common object files required by the programs
# We wrap them in addprefix to place them inside the build directory
COMMON_NAMES = sps_vars.o sps_utils.o compsp.o csp_gen.o ssp_gen.o \
	spec_mags.o interp_locate.o integrate_funcint.o sps_setup.o cosmo_pz_convol.o \
	cosmo_tuniv.o integrate_sfhw.o imf.o imf_weight.o dust_add.o \
	spec_get.o spec_sbf.o blue_stragglers.o hb_mod.o remnants_add.o spec_indices.o \
	spec_smooth.o gb_mod.o nebular_add.o xrb_add.o write_isochrone.o \
	sfh_stats.o interp_linear.o integrate_tsum.o dust_agb.o interp_array.o \
	interp_zt.o vacair_conv.o igm_absorb.o cosmo_lumdist.o dust_attenuation.o \
	sfh_weight.o sfh_limit.o sfh_info.o sfh_tabular.o dust_agn.o \
	fsps_c_driver.o

COMMON_OBJS = $(addprefix $(BUILD_DIR)/, $(COMMON_NAMES))

# ===================================
# Rules
# ===================================

.PHONY: all clean shared test test_c_driver test_c

all: $(PROGS)

# --- Compilation Rules ---

# Ensure build directory exists before compiling
$(BUILD_DIR):
	@mkdir -p $(BUILD_DIR)

# Pattern rule: Compile any .f90 found in VPATH to .o in BUILD_DIR
$(BUILD_DIR)/%.o: %.f90 | $(BUILD_DIR)
	$(FC) $(FCFLAGS) -c $< -o $@

# Specific dependencies to enforce compilation order

# sps_utils.o specifically depends on sps_vars.o
$(BUILD_DIR)/sps_utils.o: $(BUILD_DIR)/sps_vars.o

# All other common objects depend on both vars and utils.
REST_OF_COMMON = $(filter-out $(BUILD_DIR)/sps_vars.o $(BUILD_DIR)/sps_utils.o, $(COMMON_OBJS))

$(REST_OF_COMMON): $(BUILD_DIR)/sps_vars.o $(BUILD_DIR)/sps_utils.o

# Main program objects also wait for modules
$(BUILD_DIR)/sps_program_simple.o $(BUILD_DIR)/sps_program_lesssimple.o $(BUILD_DIR)/sps_program_autosps.o $(BUILD_DIR)/sps_program_spec_bin.o: $(BUILD_DIR)/sps_vars.o $(BUILD_DIR)/sps_utils.o

# Dependencies for test objects
$(BUILD_DIR)/generate_test_data.o: $(BUILD_DIR)/sps_vars.o $(BUILD_DIR)/sps_utils.o
$(BUILD_DIR)/test_runner.o: $(BUILD_DIR)/sps_vars.o $(BUILD_DIR)/sps_utils.o

# --- Linking Rules ---

autosps: $(BUILD_DIR)/sps_program_autosps.o $(COMMON_OBJS)
	$(FC) $(FCFLAGS) -o $@ $^

simple: $(BUILD_DIR)/sps_program_simple.o $(COMMON_OBJS)
	$(FC) $(FCFLAGS) -o $@ $^

lesssimple: $(BUILD_DIR)/sps_program_lesssimple.o $(COMMON_OBJS)
	$(FC) $(FCFLAGS) -o $@ $^

spec_bin: $(BUILD_DIR)/sps_program_spec_bin.o $(BUILD_DIR)/sps_vars.o
	$(FC) $(FCFLAGS) -o $@ $^

# --- Shared Library Target ---

shared: $(COMMON_OBJS)
	$(FC) $(FFLAGS) -shared -o $(BUILD_DIR)/libfsps.so $^

# --- Test Targets ---

generate_test_data: $(BUILD_DIR)/generate_test_data.o $(COMMON_OBJS)
	$(FC) $(FCFLAGS) -o $@ $^

test_runner: $(BUILD_DIR)/test_runner.o $(COMMON_OBJS)
	$(FC) $(FCFLAGS) -o $@ $^

test: test_runner

test_c_driver: shared
	gcc -O3 -o test_fsps $(TEST_DIR)/test_c_driver.c -Iinclude -L$(BUILD_DIR) -lfsps -lgfortran -Wl,-rpath,$(BUILD_DIR)

test_c: test_c_driver

# --- Utilities ---

clean:
	rm -rf $(BUILD_DIR) $(PROGS) generate_test_data test_runner test_fsps
