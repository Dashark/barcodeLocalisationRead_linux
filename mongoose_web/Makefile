CXX = g++
CXXFLAGS += -std=c++11 -Wall -g
CPPFLAGS += `pkg-config --cflags hiredis glib-2.0 opencv4`
LDFLAGS  += `pkg-config --libs   hiredis glib-2.0 opencv4` -lpthread

SRC += main.cpp mongoose.cpp

OBJS += $(patsubst %.cpp, %.o, $(SRC))

all: main

%.o: %.cpp
	@echo ">>> Compile $@ begin >>>"
	$(CXX) -c -o $(@) $(CXXFLAGS) $(CPPFLAGS) $<
	@echo "<<< Compile $@ end <<<\n"

main: $(OBJS)
	# @echo ">>> Link $(OBJS) to main begin >>>"
	$(CXX) -o $(@) $(CXXFLAGS) $(OBJS) $(LDFLAGS)
	# @echo "<<< Link $(OBJS) to main end <<<"

clean:
	@echo "*** Removing linked and compiled files ***"
	rm -f $(OBJS)
	rm -f main

