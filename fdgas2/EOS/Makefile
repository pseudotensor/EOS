CXX=g++
CFLAGS=-Wall -O3 --fast-math




default: eos_tester

clean:
	rm eos_tester *.o

eos_tester: main.o numericalrecipes.o
	$(CXX) $(CFLAGS) -o eos_tester main.o numericalrecipes.o

main.o: main.cpp fermidirac.h fermidiracintegrator.h numericalrecipes.h constants.h
	$(CXX) $(CFLAGS) -c main.cpp

numericalrecipes.o: numericalrecipes.h
	$(CXX) $(CFLAGS) -c numericalrecipes.cpp

