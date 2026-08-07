# Prefer NVIDIA HPC compiler if available, otherwise fall back
CXX := $(shell \
	if command -v nvc++ >/dev/null 2>&1; then echo nvc++; \
	elif command -v clang++ >/dev/null 2>&1; then echo clang++; \
	else echo g++; fi)

# Common flags
COMMON_FLAGS := -std=c++17 -Wall

# Parallel backend for the ray_step kernels. Defaults match the previous
# hardcoded behavior per compiler (openacc for nvc++/g++, none for clang++,
# since plain clang++ has neither -fopenacc nor a bundled -fopenmp runtime).
# Override with e.g. `make BACKEND=openmp` to try CPU multicore instead.
BACKEND ?= $(if $(findstring nvc++,$(CXX)),openacc,$(if $(findstring clang++,$(CXX)),none,openacc))

ifeq ($(BACKEND),openacc)
    BACKEND_DEFINE := -DUSE_OPENACC
    ifeq ($(findstring nvc++,$(CXX)),nvc++)
        BACKEND_FLAGS := -acc=gpu -gpu=ccall
    else
        BACKEND_FLAGS := -fopenacc
    endif
else ifeq ($(BACKEND),openmp)
    BACKEND_DEFINE := -DUSE_OPENMP
    ifeq ($(shell uname -s),Darwin)
        # Apple clang accepts neither plain -fopenmp nor bundles an OpenMP
        # runtime; libomp (MacPorts or Homebrew) supplies both the headers
        # and libomp.dylib, wired in via -Xpreprocessor per Apple's own
        # (and MacPorts' libomp port) documented recipe.
        LIBOMP_PREFIX := $(firstword \
            $(wildcard /opt/homebrew/opt/libomp) \
            $(wildcard /usr/local/opt/libomp) \
            $(wildcard /opt/local/lib/libomp))
        ifeq ($(LIBOMP_PREFIX),/opt/local/lib/libomp)
            LIBOMP_INC := /opt/local/include/libomp
            LIBOMP_LIB := /opt/local/lib/libomp
        else ifneq ($(LIBOMP_PREFIX),)
            LIBOMP_INC := $(LIBOMP_PREFIX)/include
            LIBOMP_LIB := $(LIBOMP_PREFIX)/lib
        else
            $(error BACKEND=openmp needs libomp: install it via `port install libomp +top_level` or `brew install libomp`)
        endif
        BACKEND_FLAGS := -Xpreprocessor -fopenmp -I$(LIBOMP_INC) -L$(LIBOMP_LIB) -lomp
    else
        BACKEND_FLAGS := -fopenmp
    endif
else ifeq ($(BACKEND),none)
    BACKEND_DEFINE :=
    BACKEND_FLAGS :=
else
    $(error Unknown BACKEND '$(BACKEND)': use openacc, openmp, or none)
endif

# Compiler-specific optimization flags
ifeq ($(findstring nvc++,$(CXX)),nvc++)
    # NVHPC's combined gpu,multicore lowering can fail for the adaptive
    # integrator; GPU-only mode is both the intended HPC target and reliable.
    CFLAGS := $(COMMON_FLAGS) -O4 -fast $(BACKEND_FLAGS) $(BACKEND_DEFINE)
else
    CFLAGS := $(COMMON_FLAGS) -O3 $(BACKEND_FLAGS) $(BACKEND_DEFINE)
endif

HEADERS := $(wildcard *.h)
# ffi_bridge.cpp is built separately (as a shared lib, with -fPIC) - keep it
# out of the main CLI executable's object/link list.
SOURCES := $(filter-out ffi_bridge.cpp, $(wildcard *.cpp))
OBJECTS := $(SOURCES:.cpp=.o)

TARGET := main
LIBTARGET := libblackhole.so

# Default target
all: $(TARGET) $(LIBTARGET)

# Link
$(TARGET): $(OBJECTS)
	$(CXX) $(CFLAGS) -o $@ $^ -lm

# The persistent FFI worker (web_server.py) dlopens this once and keeps
# rendering through it, instead of paying CUDA context init per frame.
$(LIBTARGET): ffi_bridge.cpp $(HEADERS)
	$(CXX) $(CFLAGS) -shared -fPIC -o $@ ffi_bridge.cpp -lm

# Compile with header dependencies
%.o: %.cpp $(HEADERS)
	$(CXX) $(CFLAGS) -MMD -MP -c $< -o $@

# Auto-generated dependency files
-include $(OBJECTS:.o=.d)

# Header-only unit tests. Built without the parallel backend: they exercise
# per-ray helpers directly rather than the ray_step kernels, so there is nothing
# to offload and no reason to need nvc++ to run them.
TESTS := $(basename $(wildcard tests/test_*.cpp))

test: $(TESTS)
	@for t in $(TESTS); do echo "== $$t"; ./$$t || exit 1; done

tests/test_%: tests/test_%.cpp $(HEADERS)
	$(CXX) $(COMMON_FLAGS) -O2 -o $@ $< -lm

# Cleanup
clean:
	$(RM) $(OBJECTS) $(OBJECTS:.o=.d) $(TARGET) $(LIBTARGET) $(TESTS)

.PHONY: all test clean

# Example NVIDIA flags:
# -gpu=cc61,fastmath -Mipa=fast -Minfo=all
# SDL2 include example:
# -I/usr/include/SDL2
