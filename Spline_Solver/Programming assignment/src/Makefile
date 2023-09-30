all:	spline.h ProblemA.cpp ProblemCD.cpp ProblemE.cpp ProblemA.json ProblemCD.json ProblemE.json report.tex
	g++ -o A ProblemA.cpp -ljsoncpp -std=c++11 -O3
	g++ -o CD ProblemCD.cpp -ljsoncpp -std=c++11 -O3
	g++ -o E ProblemE.cpp -ljsoncpp -std=c++11 -O3
	xelatex report.tex
	rm -f report.log report.aux
debug:	spline.h ProblemA.cpp ProblemCD.cpp ProblemE.cpp ProblemA.json ProblemCD.json ProblemE.json
	g++ -o A ProblemA.cpp -ljsoncpp -std=c++11 -g
	g++ -o CD ProblemCD.cpp -ljsoncpp -std=c++11 -g
	g++ -o E ProblemE.cpp -ljsoncpp -std=c++11 -g
report:	report.tex
	xelatex report.tex
	rm -f report.log report.aux
clean:
	rm -f A CD E ProblemA_result_error.txt ProblemA_result_points.txt ProblemC_result_points.txt ProblemD_result_error.txt ProblemD_result_errorpoints.txt ProblemE_result_points.txt 
