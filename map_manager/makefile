CXX = g++

MODULES = main map_manager
OBJECTS = $(addsuffix .o, $(MODULES))

CV_LIB_DIR = /usr/local/lib
CV_INC_DIR = /usr/local/include
CV_FLAGS = -I$(CV_INC_DIR) -L$(CV_LIB_DIR) -lopencv_core -lopencv_highgui

map_manager: $(OBJECTS)
	$(CXX) $(CV_FLAGS) $(OBJECTS) -o map_manager

%.o: %.cpp
	$(CXX) -c $< -o $@

.PHONY: clean
clean:
	rm *.o
