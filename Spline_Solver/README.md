# Spline Solver Project

This project implements an spline solver using piecewise-polynomial splines and B-splines. 

## File Structure

- `spline.h`: This header file defines two classes, `Bspline_interpolation` and `ppForm_interpolation`, which respectively implement B-splines and piecewise-polynomial splines for interpolation. Both classes provide methods for performing first and third degree spline interpolation. For cubic spline interpolation, it supports three types of boundary conditions: complete cubic spline (CCS), cubic spline with specified second derivatives at its end points (CS-SSD), and natural cubic spline (NCS).
- `ProblemA`,`ProblemCD`,`ProblemE`: Three test programs.

## Test Programs

- `ProblemA`: Tests all the spline interpolation methods of the function f(x)=1/(1+25x^2) in the interval [-1,1]. Visualizes the results and reports the convergence rates with respect to the number of subintervals.
- `ProblemCD`: Tests B-spline interpolation of the function f(x)=1/(1+x^2) in the interval [-5,5] and analyzes interpolation errors at different knots.
- `ProblemE`: Uses B-splines for interpolation of a closed planar curve in the shape of heart.

## Compilation and Execution

This project utilizes the `eigen3` library for solving systems of linear equations and uses the `jsoncpp` library for parameter input. These libraries are included in the project as `<eigen3/Eigen/...>` and `<jsoncpp/json/json.h>` respectively. Therefore, please ensure the correct file relationships while compiling.

You can compile the code by typing `make` in the directory where the `Makefile` is located. This will create the executable files `A`, `CD`, `E` for `problemA` through `problemE` respectively. Run each one to get the output for each sub-question.

```bash
# Compile
make

# Run
./A
./CD
./E
