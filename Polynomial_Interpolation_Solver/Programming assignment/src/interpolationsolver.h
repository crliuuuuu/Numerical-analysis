#include <iostream>
#include <cmath>
#include <algorithm>
#include <limits>
#include <vector>

double eps = std::numeric_limits<float>::epsilon();


class Function
{
public:
    virtual double operator()(double _x) = 0;
};


class Interpolation
{
public:
    virtual double solve(double _x) = 0; 
};

/*
    Newton polynomial
*/
class Newton_polynomial : public Function
{
private:
    std::vector<double> C, P;
    int n;
public:
    Newton_polynomial(std::vector<double> _C, std::vector<double> _P): C(_C), P(_P) {}
    double operator()(double _x)
    {
        double y = C[0];
        n = P.size();
        for (int i = 1; i < n; i++)
        {
            double delta_y = C[i];
            for (int j = 0; j < i; j++)
            {
                delta_y = delta_y*(_x - P[j]);
            }
            y = y + delta_y;
        }
        return y;
    }
};

/*
    Hermite polynomial
*/
class Hermite_polynomial : public Function
{
private:
    std::vector<double> C, L;
    int n;
public:
    Hermite_polynomial(std::vector<double> _C, std::vector<double> _L): C(_C), L(_L) {}
    double operator()(double _x)
    {
        double y = C[0];
        n = L.size();
        for (int i = 1; i < n; i++)
        {
            double delta_y = C[i];
            for (int j = 0; j < i; j++)
            {
                delta_y = delta_y*(_x - L[j]);
            }
            y = y + delta_y;
        }
        return y;
    }
};

/*
    Interpolation using Newton's Interpolation Formula
*/
class Newton_Formula : public Interpolation
{
private:
    std::vector<double> X, Coef;
    Function &f;
    int n;
public:
    // Constructor function, after inputting the function _f to be interpolated and the interpolation point _X, 
    // use the difference quotient to find the coefficients of the Newton interpolation polynomial
    Newton_Formula(Function &_f, std::vector<double> _X): f(_f), X(_X) 
    {
        n = X.size();
        double table[n][n];
        for (int i = 0; i < n; i++)
        {
            table[i][0] = f(X[i]);
        }
        for (int i = 1; i < n; i++)
        {
            for (int j = i; j < n; j++)
            {
                table[j][i] = (table[j][i-1] - table[j-1][i-1])/(X[j] - X[j-i]);
            }
        }
        for (int i = 0; i < n; i++)
        {
            Coef.push_back(table[i][i]);
        }
    }
    // solve() function returns the polynomial value corresponding to any point of the Newton interpolation polynomial
    double solve(double _x)
    {
        Newton_polynomial Pn(Coef, X);
        return Pn(_x);
    }
    /*
    If needed, the Get_coef() function returns the corresponding coefficients of the Newton interpolation polynomial
    */
    std::vector<double> Get_coef()
    {
        return Coef;
    }
};

/*
    Interpolation using Hermite Interpolation Formula, inherits from Interpolation class
*/
class Hermite_Formula : public Interpolation
{
private:
    std::vector<double> X, Label, Coef;
    int n;
public:
    /*
    Constructor function, after inputting the interpolation point _Label and the function (or its derivative value) _X corresponding to each interpolation point,
    use the difference quotient to find the coefficients of the Hermite interpolation polynomial
    Note that the input order of the interpolation points in _Label and _X needs to meet: 
    1. The same interpolation points need to be adjacent; 
    2. The function or its derivative values corresponding to the same interpolation points should be arranged in ascending order of derivative order.
    */
    Hermite_Formula(std::vector<double> _Label, std::vector<double> _X): Label(_Label), X(_X)
    {
        n = X.size();
        double table[n][n] = {{0}};
        int count = 1;
        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < count; j++)
            {
                table[i][j] = X[i+j+1-count];
            }
            if (Label[i] == Label[i+1])
                count = count + 1;
            else
                count = 1;
        }
        for (int i = 1; i < n; i++)
        {
            for (int j = i; j < n; j++)
            {
                if (table[j][i] == 0)
                    table[j][i] = (table[j][i-1] - table[j-1][i-1])/(Label[j] - Label[j-i]);
            }
        }
        for (int i = 0; i < n; i++)
        {
            Coef.push_back(table[i][i]);
        }
    }
    // the solve() function returns the polynomial value corresponding to any point of the Hermite interpolation polynomial
    double solve(double _x)
    {
        Newton_polynomial Pn(Coef, Label);
        return Pn(_x);
    }
    // Due to the problem requirements, the solvediff() function returns the derivative of the interpolation polynomial at the point _x
    double solvediff(double _x)
    {
        Newton_polynomial Pn(Coef, Label);
        return (Pn(_x + eps) - Pn(_x - eps))/(2.0*eps);
    }
    /*
    If needed, the Get_coef() function returns the corresponding coefficients of the Hermite interpolation polynomial
    */
    std::vector<double> Get_coef()
    {
        return Coef;
    }
};
