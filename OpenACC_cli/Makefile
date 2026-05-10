# Prefer NVIDIA HPC compiler if available, otherwise fall back
CXX := $(shell \
	if command -v nvc++ >/dev/null 2>&1; then echo nvc++; \
	elif command -v clang++ >/dev/null 2>&1; then echo clang++; \
	else echo g++; fi)

# Common flags
COMMON_FLAGS := -std=c++17 -Wall

# Compiler-specific optimization flags
ifeq ($(findstring nvc++,$(CXX)),nvc++)
    CFLAGS := $(COMMON_FLAGS) -O4 -fast -acc=gpu,multicore -gpu=ccall
else ifeq ($(findstring clang++,$(CXX)),clang++)
    CFLAGS := $(COMMON_FLAGS) -O3
else
    CFLAGS := $(COMMON_FLAGS) -O3 -fopenacc
endif

HEADERS := $(wildcard *.h)
SOURCES := $(wildcard *.cpp)
OBJECTS := $(SOURCES:.cpp=.o)

TARGET := main

# Default target
all: $(TARGET)

# Link
$(TARGET): $(OBJECTS)
	$(CXX) $(CFLAGS) -o $@ $^ -lm

# Compile with header dependencies
%.o: %.cpp $(HEADERS)
	$(CXX) $(CFLAGS) -MMD -MP -c $< -o $@

# Auto-generated dependency files
-include $(OBJECTS:.o=.d)

# Cleanup
clean:
	$(RM) $(OBJECTS) $(OBJECTS:.o=.d) $(TARGET)

# Example NVIDIA flags:
# -gpu=cc61,fastmath -Mipa=fast -Minfo=all
# SDL2 include example:
# -I/usr/include/SDL2
