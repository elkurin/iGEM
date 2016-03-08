CXX_OPT=-std=c++11 -Wall -Wextra -Wconversion

a.out: main.o equation.o eular_runge.o
	g++ $(CXX_OPT) -o a.out main.o equation.o eular_runge.o

main.o: main.cpp
	g++ $(CXX_OPT) -c main.cpp

equation.o: equation.cpp
	g++ $(CXX_OPT) -c equation.cpp

eular_runge.o: eular_runge.cpp
	g++ $(CXX_OPT) -c eular_runge.cpp

.PHONY: clean
clean:
	rm -f a.out main.o equation.o eular_runge.o
