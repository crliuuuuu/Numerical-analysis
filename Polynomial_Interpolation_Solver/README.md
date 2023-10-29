# Polynomial Interpolation Solver Project

This project implements an polynomial interpolation solver using Newton polynomial and Hermite polynomial. 

## File Structure

- `interpolationsolver.h`: Defines `Function` and `Interpolation` classes for handling different types of functions and interpolation formulas. Includes derived classes for Newton and Hermite polynomials, along with classes for calculating coefficients using DIVIDED DIFFERENCES.
- `ProblemB.cpp` - `ProblemE.cpp`: Four test programs.

## Test Programs

- `ProblemB.cpp`: Tests the interpolation of the function f(x)=1/(1+x^2) at uniformly distributed points in the interval [-5,5] to demonstrate the Runge's phenomenon.
- `ProblemC.cpp`: Uses Chebyshev interpolation to avoid the Runge's phenomenon observed in `ProblemB`.
- `ProblemD.cpp`: Uses Hermite interpolation to predict the position of a vehicle on the road.
- `ProblemE.cpp`: Uses Newton's formula to predict the survival time of larvae, and analyzes the prediction results.

## Compilation and Execution

You can compile the code by typing `make` in the directory where the `Makefile` is located. This will create the executable files `B`, `C`, `D`, `E` for `problemB.cpp` through `problemE.cpp` respectively. Run each one to get the output for each sub-question.

```bash
# Compile
make

# Run
./B
./C
./D
./E
