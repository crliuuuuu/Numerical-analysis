# Spline Solver Project

This project implements a spline solver using piecewise-polynomial splines and B-splines. 

## File Structure

- `spline.h`: This header file defines two classes, `Bspline_interpolation` and `ppForm_interpolation`, which respectively implement B-splines and piecewise-polynomial splines for interpolation. Both classes provide methods for performing first and third degree spline interpolation. For cubic spline interpolation, it supports three types of boundary conditions: complete cubic spline (CCS), cubic spline with specified second derivatives at its end points (CS-SSD), and natural cubic spline (NCS).
- `ProblemA.cpp`,`ProblemCD.cpp`,`ProblemE.cpp`: Three test programs.
- `ProblemA.json`, `ProblemCD.json`, `ProblemE.json`: Input parameter files for the corresponding test programs.

## Test Programs

- `ProblemA.cpp`: Tests all the spline interpolation methods of the function f(x)=1/(1+25x^2) in the interval [-1,1]. Visualizes the results and reports the convergence rates with respect to the number of subintervals.
- `ProblemCD.cpp`: Tests B-spline interpolation of the function f(x)=1/(1+x^2) in the interval [-5,5] and analyzes interpolation errors at different knots.
- `ProblemE.cpp`: Uses B-splines for interpolation of a closed planar curve in the shape of heart.

## Input Parameter Files

These JSON files provide parameters for the tests. Here are the descriptions of the parameters in each file:

**ProblemA.json**

- `interval_left`: The left endpoint of the interpolation interval. 
- `interval_right`: The right endpoint of the interpolation interval.
- `n`: Number of interpolation points.
- `method`: Input 1, 2, or 3, representing the interpolation boundary conditions. For cubic spline, 1 represents complete cubic spline, 2 represents cubic spline with specified second derivatives, 3 represents natural cubic spline; for linear spline, input 1-3 does not make a difference.
- `order`: The order of the spline. Input 1 or 3, representing linear and cubic splines, respectively.

**ProblemCD.json**

- `point_1`: Points for cubic Bspline interpolation. 
- `point_2`: Points for linear Bspline interpolation. 
- `output_point`: Points where the error function needs to output values. N

**ProblemE.json**

- `n`: Number of interpolation points. 
- `order`: The order of the spline. Input 1 or 3, representing linear and cubic splines, respectively.
- `spline_form`: Interpolation method. Input "Bspline" or "PP" to choose Bspline interpolation or ppForm interpolation respectively.

## Compilation and Execution

This project utilizes the `eigen3` library for solving systems of linear equations and uses the `jsoncpp` library for parameter input. These libraries are included in the project as `<eigen3/Eigen/...>` and `<jsoncpp/json/json.h>` respectively. Therefore, please ensure the correct file relationships while compiling.

You can compile the code by typing `make` in the directory where the `Makefile` is located. This will create the executable files `A`, `CD`, `E` for `problemA.cpp` through `problemE.cpp` respectively. Run each one to get the output for each sub-question.

```bash
# Compile
make

# Run
./A
./CD
./E
```

## Dependencies

This project depends on the following libraries:

1. [Eigen3]: A C++ template library for linear algebra. To install Eigen3, you can follow the instructions on their official website.

2. [JsonCpp]: A C++ library for interacting with JSON. You can install it using the following command:

    For Ubuntu:

    ```bash
    sudo apt-get install libjsoncpp-dev
    ```
