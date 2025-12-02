# Makefile for Grid-Based FMM with Cartesian Multipole Expansion
#
# Supports:
# - gfortran (GNU Fortran)
# - ifort (Intel Fortran)

# Compiler selection (default: gfortran)
FC = gfortran

# Compiler flags
ifeq ($(FC),gfortran)
    FFLAGS = -O3 -march=native -funroll-loops -ffast-math -Wall -Wextra
    FFLAGS_DEBUG = -g -O0 -Wall -Wextra -fcheck=all -fbacktrace
else ifeq ($(FC),ifort)
    FFLAGS = -O3 -xHost -ipo -warn all
    FFLAGS_DEBUG = -g -O0 -warn all -check all -traceback
else
    FFLAGS = -O2
    FFLAGS_DEBUG = -g -O0
endif

# Directories
SRC_DIR = src
BUILD_DIR = build
BIN_DIR = bin

# Source files (in dependency order)
SOURCES = $(SRC_DIR)/types_mod.f90 \
          $(SRC_DIR)/potential_interface_mod.f90 \
          $(SRC_DIR)/cartesian_multipole_mod.f90 \
          $(SRC_DIR)/grid_fmm_mod.f90 \
          $(SRC_DIR)/md_utils_mod.f90 \
          $(SRC_DIR)/test_grid_fmm.f90

# Object files
OBJECTS = $(patsubst $(SRC_DIR)/%.f90,$(BUILD_DIR)/%.o,$(SOURCES))

# Module files
MODULES = $(BUILD_DIR)/types_mod.mod \
          $(BUILD_DIR)/potential_interface_mod.mod \
          $(BUILD_DIR)/cartesian_multipole_mod.mod \
          $(BUILD_DIR)/grid_fmm_mod.mod \
          $(BUILD_DIR)/md_utils_mod.mod

# Executable
TARGET = $(BIN_DIR)/test_grid_fmm

# Default target
all: directories $(TARGET)

# Create directories
directories:
	@mkdir -p $(BUILD_DIR)
	@mkdir -p $(BIN_DIR)

# Link executable
$(TARGET): $(OBJECTS)
	@echo "Linking $@..."
	$(FC) $(FFLAGS) -J$(BUILD_DIR) -o $@ $(OBJECTS)
	@echo "Build complete: $@"

# Compile object files
$(BUILD_DIR)/types_mod.o: $(SRC_DIR)/types_mod.f90
	@echo "Compiling $<..."
	$(FC) $(FFLAGS) -J$(BUILD_DIR) -c $< -o $@

$(BUILD_DIR)/potential_interface_mod.o: $(SRC_DIR)/potential_interface_mod.f90 \
                                        $(BUILD_DIR)/types_mod.o
	@echo "Compiling $<..."
	$(FC) $(FFLAGS) -J$(BUILD_DIR) -c $< -o $@

$(BUILD_DIR)/cartesian_multipole_mod.o: $(SRC_DIR)/cartesian_multipole_mod.f90 \
                                        $(BUILD_DIR)/types_mod.o
	@echo "Compiling $<..."
	$(FC) $(FFLAGS) -J$(BUILD_DIR) -c $< -o $@

$(BUILD_DIR)/grid_fmm_mod.o: $(SRC_DIR)/grid_fmm_mod.f90 \
                             $(BUILD_DIR)/types_mod.o \
                             $(BUILD_DIR)/potential_interface_mod.o \
                             $(BUILD_DIR)/cartesian_multipole_mod.o
	@echo "Compiling $<..."
	$(FC) $(FFLAGS) -J$(BUILD_DIR) -c $< -o $@

$(BUILD_DIR)/md_utils_mod.o: $(SRC_DIR)/md_utils_mod.f90 \
                             $(BUILD_DIR)/types_mod.o
	@echo "Compiling $<..."
	$(FC) $(FFLAGS) -J$(BUILD_DIR) -c $< -o $@

$(BUILD_DIR)/test_grid_fmm.o: $(SRC_DIR)/test_grid_fmm.f90 \
                              $(BUILD_DIR)/types_mod.o \
                              $(BUILD_DIR)/potential_interface_mod.o \
                              $(BUILD_DIR)/grid_fmm_mod.o \
                              $(BUILD_DIR)/md_utils_mod.o
	@echo "Compiling $<..."
	$(FC) $(FFLAGS) -J$(BUILD_DIR) -c $< -o $@

# Debug build
debug: FFLAGS = $(FFLAGS_DEBUG)
debug: clean all

# Clean build artifacts
clean:
	@echo "Cleaning build artifacts..."
	rm -rf $(BUILD_DIR) $(BIN_DIR)
	rm -f *.mod *.o *.csv
	@echo "Clean complete."

# Run test
run: $(TARGET)
	@echo "Running test..."
	cd $(BIN_DIR) && ./test_grid_fmm

# Phony targets
.PHONY: all directories debug clean run
