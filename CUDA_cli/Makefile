CXX := nvcc
CFLAGS := -use_fast_math -m64 -std=c++17

HEADERS := $(wildcard *.h)
SOURCES := $(wildcard *.cu)

 
main: $(SOURCES) cli_parser.cpp cli_parser.h cuda_ray.cuh
	$(CXX) $(CFLAGS) -o main $(SOURCES) cli_parser.cpp -lm
 

 
clean:
	$(RM) main
