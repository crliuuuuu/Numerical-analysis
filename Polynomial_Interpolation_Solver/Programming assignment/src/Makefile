all:	interpolationsolver.h problemB.cpp problemC.cpp problemD.cpp problemE.cpp 
	g++ -std=c++11 -o B problemB.cpp -O3
	g++ -std=c++11 -o C problemC.cpp -O3
	g++ -std=c++11 -o D problemD.cpp -O3
	g++ -std=c++11 -o E problemE.cpp -O3
debug:	equationsolver.h problemB.cpp problemC.cpp problemD.cpp problemE.cpp 
	g++ -std=c++11 -o B problemB.cpp -g
	g++ -std=c++11 -o C problemC.cpp -g
	g++ -std=c++11 -o D problemD.cpp -g
	g++ -std=c++11 -o E problemE.cpp -g
clean:
	rm -f B C D E 
