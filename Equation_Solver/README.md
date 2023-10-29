# Equation Solver Project

This project implements an equation solver using three different methods: the bisection method, Newton's method, and the secant method. 

## File Structure

- `equationsolver.h`: Includes the derived classes of `EquationSolver`, implementing the bisection method, Newton's method, and the secant method.
- `ProblemB.cpp` - `ProblemF.cpp`: Five test programs.

## Test Programs

- `ProblemB.cpp`: Performs root finding tests using the bisection method.
- `ProblemC.cpp`: Performs root finding tests using Newton's method.
- `ProblemD.cpp`: Performs root finding tests using the secant method.
- `ProblemE.cpp`: Solves a problem of finding the depth of water in a trough using all three methods.
- `ProblemF.cpp`: Solves a problem of finding the maximum angle to avoid hang-up failure in all-terrain vehicle design using all three methods.

## Compilation and Execution

You can compile the code by typing `make` in the directory where the `Makefile` is located. This will create the executable files `B`, `C`, `D`, `E`, `F` for `problemB.cpp` through `problemF.cpp` respectively. Run each one to get the output for each sub-question.

```bash
# Compile
make

# Run
./B
./C
./D
./E
./F
