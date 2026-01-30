# ===================================
# Compiler Configuration
# ===================================

# Default to gfortran, allow for an override (e.g. make FC=ifort)
FC = gfortran

# Default flags
FFLAGS ?= -O3 -cpp -fPIC

# Installation paths (override on make command line as needed)
PREFIX ?= /usr/local
DESTDIR ?=
LIBDIR ?= $(PREFIX)/lib
INCLUDEDIR ?= $(PREFIX)/include
DATADIR ?= $(PREFIX)/share
PKGCONFIGDIR ?= $(LIBDIR)/pkgconfig
BINDIR ?= $(PREFIX)/bin

# Shared library versioning
LIB_VERSION ?= 3.2.0
LIB_SONAME ?= 3

# Directory configuration
SRC_DIR := src
TEST_DIR := tests
BUILD_DIR := build
TEST_BUILD_DIR := $(BUILD_DIR)/tests
SRC_SUBDIRS := core programs spectra sfh physics cosmology math io imf abi
SRC_DIRS := $(addprefix $(SRC_DIR)/, $(SRC_SUBDIRS))

# Module directory flags:
# Gfortran uses -J to specify where to put/find .mod files
# Intel (ifort) uses -module. We default to Gfortran syntax here.
# You can override this: make MOD_FLAG=-module
MOD_FLAG ?= -J

# Combine flags to include the build dir for .mod search
FCFLAGS := $(FFLAGS) -fPIC $(MOD_FLAG)$(BUILD_DIR) -I$(BUILD_DIR)

# ===================================
# Source and Object Definitions
# ===================================

# Tell make to look for source files in these directories
VPATH = $(SRC_DIRS):$(TEST_DIR)

# The list of programs to build (executables)
PROGS = simple lesssimple autosps spec_bin

# The common object files required by the programs
# We wrap them in addprefix to place them inside the build directory
COMMON_NAMES = fsps_types.o sps_vars.o fsps_cache.o sps_utils.o fsps_context_types.o compsp.o csp_gen.o ssp_gen.o \
	fsps_context.o \
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


.PHONY: all clean shared test test_cache test_c_driver test_c_contexts test_c check install uninstall \
	install-lib install-headers install-bin install-pkgconfig install-data

all: $(PROGS)

# --- Compilation Rules ---

# Ensure build directory exists before compiling
$(BUILD_DIR):
	@mkdir -p $(BUILD_DIR)

# Ensure test build directory exists
$(TEST_BUILD_DIR): | $(BUILD_DIR)
	@mkdir -p $(TEST_BUILD_DIR)

# Pattern rule: Compile any .f90 found in VPATH to .o in BUILD_DIR
$(BUILD_DIR)/%.o: %.f90 | $(BUILD_DIR)
	$(FC) $(FCFLAGS) -c $< -o $@

# Specific dependencies to enforce compilation order

# Module dependency ordering
$(BUILD_DIR)/sps_vars.o: $(BUILD_DIR)/fsps_types.o

$(BUILD_DIR)/fsps_cache.o: $(BUILD_DIR)/sps_vars.o

$(BUILD_DIR)/fsps_context_types.o: $(BUILD_DIR)/sps_vars.o $(BUILD_DIR)/fsps_cache.o

# sps_utils.o specifically depends on sps_vars.o and context types
$(BUILD_DIR)/sps_utils.o: $(BUILD_DIR)/sps_vars.o $(BUILD_DIR)/fsps_cache.o $(BUILD_DIR)/fsps_context_types.o

# All other common objects depend on vars, cache, utils, and context types.
REST_OF_COMMON = $(filter-out $(BUILD_DIR)/fsps_types.o $(BUILD_DIR)/sps_vars.o $(BUILD_DIR)/fsps_cache.o \
	$(BUILD_DIR)/sps_utils.o $(BUILD_DIR)/fsps_context_types.o, $(COMMON_OBJS))

$(REST_OF_COMMON): $(BUILD_DIR)/fsps_types.o $(BUILD_DIR)/sps_vars.o $(BUILD_DIR)/fsps_cache.o \
	$(BUILD_DIR)/sps_utils.o $(BUILD_DIR)/fsps_context_types.o

# Main program objects also wait for modules
$(BUILD_DIR)/sps_program_simple.o $(BUILD_DIR)/sps_program_lesssimple.o $(BUILD_DIR)/sps_program_autosps.o $(BUILD_DIR)/sps_program_spec_bin.o: $(BUILD_DIR)/sps_vars.o $(BUILD_DIR)/sps_utils.o

# Test object compilation (keep artifacts out of tests/)
$(TEST_BUILD_DIR)/generate_test_data.o: $(TEST_DIR)/generate_test_data.f90 | $(TEST_BUILD_DIR)
	$(FC) $(FCFLAGS) -c $< -o $@

$(TEST_BUILD_DIR)/test_runner.o: $(TEST_DIR)/test_runner.f90 | $(TEST_BUILD_DIR)
	$(FC) $(FCFLAGS) -c $< -o $@

$(TEST_BUILD_DIR)/test_cache.o: $(TEST_DIR)/test_cache.f90 | $(TEST_BUILD_DIR)
	$(FC) $(FCFLAGS) -c $< -o $@

# Dependencies for test objects
$(TEST_BUILD_DIR)/generate_test_data.o: $(BUILD_DIR)/sps_vars.o $(BUILD_DIR)/sps_utils.o $(BUILD_DIR)/fsps_context.o
$(TEST_BUILD_DIR)/test_runner.o: $(BUILD_DIR)/sps_vars.o $(BUILD_DIR)/sps_utils.o $(BUILD_DIR)/fsps_context.o $(BUILD_DIR)/fsps_context_types.o
$(TEST_BUILD_DIR)/test_cache.o: $(BUILD_DIR)/sps_vars.o $(BUILD_DIR)/sps_utils.o $(BUILD_DIR)/fsps_context_types.o

# --- Linking Rules ---

autosps: $(BUILD_DIR)/sps_program_autosps.o $(COMMON_OBJS)
	$(FC) $(FCFLAGS) -o $@ $^

simple: $(BUILD_DIR)/sps_program_simple.o $(COMMON_OBJS)
	$(FC) $(FCFLAGS) -o $@ $^

lesssimple: $(BUILD_DIR)/sps_program_lesssimple.o $(COMMON_OBJS)
	$(FC) $(FCFLAGS) -o $@ $^

spec_bin: $(BUILD_DIR)/sps_program_spec_bin.o $(BUILD_DIR)/sps_vars.o $(BUILD_DIR)/sps_utils.o
	$(FC) $(FCFLAGS) -o $@ $^

# --- Shared Library Target ---

shared: $(COMMON_OBJS)
	$(FC) $(FFLAGS) -shared -Wl,-soname,libfsps.so.$(LIB_SONAME) \
		-o $(BUILD_DIR)/libfsps.so.$(LIB_VERSION) $^
	ln -sf libfsps.so.$(LIB_VERSION) $(BUILD_DIR)/libfsps.so.$(LIB_SONAME)
	ln -sf libfsps.so.$(LIB_SONAME) $(BUILD_DIR)/libfsps.so

# --- Test Targets ---

generate_test_data: $(TEST_BUILD_DIR)/generate_test_data.o $(COMMON_OBJS)
	$(FC) $(FCFLAGS) -o $@ $^

test_runner: $(TEST_BUILD_DIR)/test_runner.o $(COMMON_OBJS)
	$(FC) $(FCFLAGS) -o $@ $^

test_cache: $(TEST_BUILD_DIR)/test_cache.o $(COMMON_OBJS)
	$(FC) $(FCFLAGS) -o $@ $^

test: test_runner test_cache

test_c_driver: shared
	gcc -O3 -o test_fsps $(TEST_DIR)/test_c_driver.c -Iinclude -L$(BUILD_DIR) -lfsps -lgfortran -Wl,-rpath,$(BUILD_DIR)

test_c_contexts: shared
	gcc -O3 -o test_contexts $(TEST_DIR)/test_c_contexts.c -Iinclude -L$(BUILD_DIR) -lfsps -lgfortran -Wl,-rpath,$(BUILD_DIR)

test_c: test_c_driver test_c_contexts

check: test_c

# --- Utilities ---

clean:
	rm -rf $(BUILD_DIR) $(PROGS) generate_test_data test_runner test_cache test_fsps test_contexts

# --- Install Targets ---

install: all shared install-lib install-headers install-bin install-pkgconfig install-data

install-lib:
	install -d $(DESTDIR)$(LIBDIR)
	install -m 755 $(BUILD_DIR)/libfsps.so.$(LIB_VERSION) $(DESTDIR)$(LIBDIR)/
	ln -sf libfsps.so.$(LIB_VERSION) $(DESTDIR)$(LIBDIR)/libfsps.so.$(LIB_SONAME)
	ln -sf libfsps.so.$(LIB_SONAME) $(DESTDIR)$(LIBDIR)/libfsps.so

install-headers:
	install -d $(DESTDIR)$(INCLUDEDIR)
	install -m 644 include/fsps.h $(DESTDIR)$(INCLUDEDIR)/

install-bin:
	install -d $(DESTDIR)$(BINDIR)
	install -m 755 $(PROGS) $(DESTDIR)$(BINDIR)/

install-pkgconfig:
	install -d $(DESTDIR)$(PKGCONFIGDIR)
	sed -e "s|^prefix=.*|prefix=$(PREFIX)|" fsps.pc > $(DESTDIR)$(PKGCONFIGDIR)/fsps.pc

install-data:
	install -d $(DESTDIR)$(DATADIR)/fsps
	cp -a data $(DESTDIR)$(DATADIR)/fsps/

uninstall:
	rm -f $(DESTDIR)$(LIBDIR)/libfsps.so.$(LIB_VERSION)
	rm -f $(DESTDIR)$(LIBDIR)/libfsps.so.$(LIB_SONAME)
	rm -f $(DESTDIR)$(LIBDIR)/libfsps.so
	rm -f $(DESTDIR)$(INCLUDEDIR)/fsps.h
	rm -f $(DESTDIR)$(PKGCONFIGDIR)/fsps.pc
	rm -f $(DESTDIR)$(BINDIR)/simple $(DESTDIR)$(BINDIR)/lesssimple \
		$(DESTDIR)$(BINDIR)/autosps $(DESTDIR)$(BINDIR)/spec_bin
	rm -rf $(DESTDIR)$(DATADIR)/fsps/data
