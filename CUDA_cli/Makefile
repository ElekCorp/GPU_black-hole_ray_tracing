CXX := nvcc
CFLAGS := -use_fast_math -m64

HEADERS := $(wildcard *.h)
SOURCES := $(wildcard *.cu)

 
main: $(SOURCES) cli_parser.cpp
	$(CXX) $(CFLAGS) -o main $(SOURCES)  cli_parser.cpp -lm
 

 
clean:
	$(RM) main
