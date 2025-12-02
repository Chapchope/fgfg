# Makefile for Grid-Based FMM 2D - ACCURATE VERSION
#
# Точный и быстрый FMM с произвольным порядком разложения

FC = gfortran
FFLAGS = -O3 -march=native -funroll-loops -ffast-math -Wall
FFLAGS_DEBUG = -g -O0 -Wall -fcheck=all -fbacktrace

SRC_DIR = src
BUILD_DIR = build
BIN_DIR = bin

MODFLAGS = -J$(BUILD_DIR)

# Default target
all: directories accurate

directories:
	@mkdir -p $(BUILD_DIR)
	@mkdir -p $(BIN_DIR)

#=============================================================================
# ACCURATE VERSION (Main)
#=============================================================================

accurate: directories $(BUILD_DIR)/types_mod.o $(BUILD_DIR)/grid_fmm_accurate.o
	@echo "Linking test_accurate..."
	$(FC) $(FFLAGS) $(MODFLAGS) -o $(BIN_DIR)/test_accurate \
		$(BUILD_DIR)/types_mod.o \
		$(BUILD_DIR)/grid_fmm_accurate.o \
		$(SRC_DIR)/test_accurate.f90
	@echo "Build complete: $(BIN_DIR)/test_accurate"

$(BUILD_DIR)/types_mod.o: $(SRC_DIR)/types_mod.f90
	@echo "Compiling $<..."
	$(FC) $(FFLAGS) $(MODFLAGS) -c $< -o $@

$(BUILD_DIR)/grid_fmm_accurate.o: $(SRC_DIR)/grid_fmm_accurate.f90 \
                                   $(BUILD_DIR)/types_mod.o
	@echo "Compiling $<..."
	$(FC) $(FFLAGS) $(MODFLAGS) -c $< -o $@

#=============================================================================
# SIMPLE VERSION (Reference)
#=============================================================================

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

# Clean
clean:
	@echo "Cleaning..."
	rm -rf $(BUILD_DIR) $(BIN_DIR)
	rm -f *.mod *.o *.csv
	@echo "Clean complete."

# Run main accurate version
run: accurate
	@echo ""
	@echo "============================================"
	@echo "  Running Accurate Grid-based FMM 2D"
	@echo "============================================"
	@echo ""
	./$(BIN_DIR)/test_accurate

# Run simple version for comparison
run-simple: simple
	@echo "Running simple FMM..."
	./$(BIN_DIR)/test_simple

# Help
help:
	@echo "Available targets:"
	@echo "  all          - Build accurate FMM (default)"
	@echo "  accurate     - Build accurate FMM"
	@echo "  simple       - Build simple FMM (reference)"
	@echo "  run          - Build and run accurate FMM"
	@echo "  run-simple   - Build and run simple FMM"
	@echo "  debug        - Build with debug flags"
	@echo "  clean        - Remove all build artifacts"
	@echo "  help         - Show this help message"

.PHONY: all directories accurate simple debug clean run run-simple help
