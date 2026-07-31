# Prefer NVIDIA HPC compiler if available, otherwise fall back
CXX := $(shell \
	if command -v nvc++ >/dev/null 2>&1; then echo nvc++; \
	elif command -v clang++ >/dev/null 2>&1; then echo clang++; \
	else echo g++; fi)

# Common flags
COMMON_FLAGS := -std=c++17 -Wall

# Compiler-specific optimization flags
ifeq ($(findstring nvc++,$(CXX)),nvc++)
    # NVHPC's combined gpu,multicore lowering can fail for the adaptive
    # integrator; GPU-only mode is both the intended HPC target and reliable.
    CFLAGS := $(COMMON_FLAGS) -O4 -fast -acc=gpu -gpu=ccall
else ifeq ($(findstring clang++,$(CXX)),clang++)
    CFLAGS := $(COMMON_FLAGS) -O3
else
    CFLAGS := $(COMMON_FLAGS) -O3 -fopenacc
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

# Cleanup
clean:
	$(RM) $(OBJECTS) $(OBJECTS:.o=.d) $(TARGET) $(LIBTARGET)

# Example NVIDIA flags:
# -gpu=cc61,fastmath -Mipa=fast -Minfo=all
# SDL2 include example:
# -I/usr/include/SDL2
