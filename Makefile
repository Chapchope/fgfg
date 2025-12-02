# Makefile for Grid-Based FMM 2D with modular architecture
#
# Version 2: Модульная архитектура с заменяемым потенциалом
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

# Module flags
MODFLAGS = -J$(BUILD_DIR)

# Default target
all: directories test_fmm_v2

# Create directories
directories:
	@mkdir -p $(BUILD_DIR)
	@mkdir -p $(BIN_DIR)

#=============================================================================
# NEW MODULAR ARCHITECTURE (v2)
#=============================================================================

# Object files for v2
OBJS_V2 = $(BUILD_DIR)/types_mod.o \
          $(BUILD_DIR)/potential_2d_mod.o \
          $(BUILD_DIR)/multipole_2d_mod.o \
          $(BUILD_DIR)/grid_fmm_2d_mod.o

# Executable for v2
test_fmm_v2: directories $(OBJS_V2)
	@echo "Linking test_fmm_v2..."
	$(FC) $(FFLAGS) $(MODFLAGS) -o $(BIN_DIR)/test_fmm_v2 $(OBJS_V2) $(SRC_DIR)/test_fmm_v2.f90
	@echo "Build complete: $(BIN_DIR)/test_fmm_v2"

# Compile rules for v2

$(BUILD_DIR)/types_mod.o: $(SRC_DIR)/types_mod.f90
	@echo "Compiling $<..."
	$(FC) $(FFLAGS) $(MODFLAGS) -c $< -o $@

$(BUILD_DIR)/potential_2d_mod.o: $(SRC_DIR)/potential_2d_mod.f90 \
                                  $(BUILD_DIR)/types_mod.o
	@echo "Compiling $<..."
	$(FC) $(FFLAGS) $(MODFLAGS) -c $< -o $@

$(BUILD_DIR)/multipole_2d_mod.o: $(SRC_DIR)/multipole_2d_mod.f90 \
                                  $(BUILD_DIR)/types_mod.o \
                                  $(BUILD_DIR)/potential_2d_mod.o
	@echo "Compiling $<..."
	$(FC) $(FFLAGS) $(MODFLAGS) -c $< -o $@

$(BUILD_DIR)/grid_fmm_2d_mod.o: $(SRC_DIR)/grid_fmm_2d_mod.f90 \
                                 $(BUILD_DIR)/types_mod.o \
                                 $(BUILD_DIR)/potential_2d_mod.o \
                                 $(BUILD_DIR)/multipole_2d_mod.o
	@echo "Compiling $<..."
	$(FC) $(FFLAGS) $(MODFLAGS) -c $< -o $@

#=============================================================================
# OLD TARGETS (for backward compatibility)
#=============================================================================

# Old simple FMM
simple: directories $(BUILD_DIR)/types_mod.o $(BUILD_DIR)/simple_fmm_2d.o
	@echo "Linking test_simple..."
	$(FC) $(FFLAGS) $(MODFLAGS) -o $(BIN_DIR)/test_simple \
		$(BUILD_DIR)/types_mod.o \
		$(BUILD_DIR)/simple_fmm_2d.o \
		$(SRC_DIR)/test_simple.f90
	@echo "Build complete: $(BIN_DIR)/test_simple"

$(BUILD_DIR)/simple_fmm_2d.o: $(SRC_DIR)/simple_fmm_2d.f90 \
                               $(BUILD_DIR)/types_mod.o
	@echo "Compiling $<..."
	$(FC) $(FFLAGS) $(MODFLAGS) -c $< -o $@

#=============================================================================
# UTILITY TARGETS
#=============================================================================

# Debug build
debug: FFLAGS = $(FFLAGS_DEBUG)
debug: clean all

# Clean build artifacts
clean:
	@echo "Cleaning build artifacts..."
	rm -rf $(BUILD_DIR) $(BIN_DIR)
	rm -f *.mod *.o *.csv
	@echo "Clean complete."

# Run new test (v2)
run: test_fmm_v2
	@echo ""
	@echo "============================================"
	@echo "  Running Grid-based FMM 2D (v2)..."
	@echo "============================================"
	@echo ""
	./$(BIN_DIR)/test_fmm_v2

# Run old simple test
run-simple: simple
	@echo "Running simple FMM test..."
	./$(BIN_DIR)/test_simple

# Help
help:
	@echo "Available targets:"
	@echo "  all          - Build test_fmm_v2 (new modular version)"
	@echo "  test_fmm_v2  - Build new modular FMM"
	@echo "  simple       - Build old simple FMM"
	@echo "  run          - Build and run test_fmm_v2"
	@echo "  run-simple   - Build and run simple test"
	@echo "  debug        - Build with debug flags"
	@echo "  clean        - Remove all build artifacts"
	@echo "  help         - Show this help message"

# Phony targets
.PHONY: all directories test_fmm_v2 simple debug clean run run-simple help
