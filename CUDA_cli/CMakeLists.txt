cmake_minimum_required(VERSION 3.18) # Minimum version for good CUDA support
project(MyCUDAProject LANGUAGES CUDA CXX)

# Set CUDA compilation flags
set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -use_fast_math -m64")

# Build a "fat binary" targeting Pascal through Hopper architectures
# This ensures the binary is portable and performs optimally across different GPUs
if(NOT DEFINED CMAKE_CUDA_ARCHITECTURES)
    set(CMAKE_CUDA_ARCHITECTURES 60 61 70 75 80 86 89 90)
endif()

# Collect all source files
file(GLOB SOURCES "*.cu" "*.cpp")
file(GLOB HEADERS "*.h")

# Create executable
add_executable(main ${SOURCES} ${HEADERS})

# Link math library
if(NOT MSVC)
    target_link_libraries(main m)
endif()

# Installation: install binary to bin/
install(TARGETS main DESTINATION bin)

# Packaging support (optional)
include(CPack)

