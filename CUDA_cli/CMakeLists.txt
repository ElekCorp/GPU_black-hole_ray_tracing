cmake_minimum_required(VERSION 3.18) # Minimum version for good CUDA support
project(MyCUDAProject LANGUAGES CUDA CXX)

# Set CUDA compilation flags
set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -use_fast_math -m64")

# Collect all source files
file(GLOB SOURCES "*.cu")
file(GLOB HEADERS "*.h")

# Create executable
add_executable(main ${SOURCES} ${HEADERS})

# Link math library
target_link_libraries(main m)

# Installation: install binary to bin/
install(TARGETS main DESTINATION bin)

# Packaging support (optional)
include(CPack)

