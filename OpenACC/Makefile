CXX := nvc++
CFLAGS := -O4 -fast -acc=gpu -gpu=cc61 -Mipa=fast -Minfo=all  -I/usr/include/SDL2

HEADERS := $(wildcard *.h)
SOURCES := $(wildcard *.cpp)
OBJECTS := $(SOURCES:%.cpp=%.o)

 
main: $(OBJECTS)
	$(CXX) $(CFLAGS) -o main $(OBJECTS) `sdl2-config --cflags --libs` -lSDL2_gfx -lSDL2_ttf -lSDL2_image -lSDL2_mixer -lm
 
%.o: %.cpp
	$(CXX) -c $(CFLAGS) $< -o $@
 
clean:
	$(RM) $(OBJECTS)
	#O4 -fast -acc=gpu -gpu=cc61,fastmath -Mipa=fast -Minfo=all  -I/usr/include/SDL2
