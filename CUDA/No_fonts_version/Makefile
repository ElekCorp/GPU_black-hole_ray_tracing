CXX := nvcc
CFLAGS := -use_fast_math -m64

HEADERS := $(wildcard *.h)
SOURCES := $(wildcard *.cu)

 
main:
	$(CXX) $(CFLAGS) -o main $(SOURCES) `sdl2-config --cflags --libs` -lSDL2_gfx -lSDL2_ttf -lSDL2_image -lSDL2_mixer -lm
 

 
clean:
	$(RM) main
