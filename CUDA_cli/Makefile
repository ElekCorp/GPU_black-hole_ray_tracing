CXX := nvcc
CFLAGS := -use_fast_math -m64

HEADERS := $(wildcard *.h)
SOURCES := $(wildcard *.cu)

 
main: $(SOURCES)
	$(CXX) $(CFLAGS) -o main $(SOURCES) -lm
 

 
clean:
	$(RM) main
