CCC = g++
#CC11 = gnu++11
CC11 = gnu++0x
EIGEN = -I $EIGEN_DIR

CCFLAGS = -O -std=$(CC11) 

PROGRAMS = testfunctions cmaes_basic cmaes 

all:
	make $(PROGRAMS)

clean :
	/bin/rm -rf $(PROGRAMS)
	/bin/rm -rf *.o *~ 

cmaes_basic: auxiliaries.o cmaes_basic.o
	$(CCC) $(CCFLAGS) $(EIGEN) auxiliaries.o cmaes_basic.o -o cmaes_basic
cmaes: auxiliaries.o cmaes.o
	$(CCC) $(CCFLAGS) $(EIGEN) auxiliaries.o cmaes.cpp -o cmaes
testfunctions: testfunctions.o 
	$(CCC) $(CCFLAGS) testfunctions.o -o testfunctions
testfunctions.o: testfunctions.cpp
	$(CCC) -c $(CCFLAGS) $(EIGEN) testfunctions.cpp -o testfunctions.o
auxiliaries.o: auxiliaries.cpp
	$(CCC) -c $(CCFLAGS) $(EIGEN) auxiliaries.cpp -o auxiliaries.o
cmaes_basic.o: cmaes_basic.cpp
	$(CCC) -c $(CCFLAGS) $(EIGEN) cmaes_basic.cpp -o cmaes_basic.o
cmaes.o: cmaes.cpp
	$(CCC) -c $(CCFLAGS) $(EIGEN) cmaes.cpp -o cmaes.o
