CXX := clang++ #nvc++
CFLAGS := -Wall -O3 -std=c++17#-O4 -fast -acc=gpu -gpu=cc61 -Mipa=fast #-Minfo=all

HEADERS := $(wildcard *.h)
SOURCES := $(wildcard *.cpp)
OBJECTS := $(SOURCES:%.cpp=%.o)

 
main: $(OBJECTS)
	$(CXX) $(CFLAGS) -o main $(OBJECTS) -lm
 
%.o: %.cpp
	$(CXX) -c $(CFLAGS) $< -o $@
 
clean:
	$(RM) $(OBJECTS) main
	#O4 -fast -acc=gpu -gpu=cc61,fastmath -Mipa=fast -Minfo=all  -I/usr/include/SDL2
