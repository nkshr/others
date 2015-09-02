CXX = g++
MODULES = primal_form steepest_descent
OBJECTS = $(addsuffix .o, $(MODULES)) 

svm: $(OBJECTS)	
	$(CXX) $(OBJECTS) -o svm

%.o: %.cpp
	$(CXX) -c $< -o $@

.PHONEY:clean
clean:
	rm *.o
