CXX := g++
CFLAGS :=-Ofast -Wall `/Users/elek/SDL2/bin/sdl2-config --cflags --libs` -I/Users/elek/SDL2/include #-I/usr/include/SDL2

HEADERS := $(wildcard *.h)
SOURCES := $(wildcard *.cpp)
OBJECTS := $(SOURCES:%.cpp=%.o)

 
main: $(OBJECTS)
	$(CXX) $(CFLAGS) -o main $(OBJECTS) -lm
%.o: %.cpp
	$(CXX) -c $(CFLAGS) $< -o $@
 
clean:
	$(RM) $(OBJECTS)
	#O4 -fast -acc=gpu -gpu=cc61,fastmath -Mipa=fast -Minfo=all  -I/usr/include/SDL2
