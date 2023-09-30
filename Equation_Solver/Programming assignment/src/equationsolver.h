#include <iostream>
#include <limits>
#include <cmath>
#include <algorithm>

double eps = std::numeric_limits<double>::epsilon();

class Function
{
public:
    virtual double operator()(double _x) = 0;
    virtual double diff(double _x) 
    {
        return 0;  //对于某些不需要用到导数的函数设置默认值为0
    }
};

class EquationSolver
{
public:
  virtual double solve() = 0; 
};

class Bisection : public EquationSolver
{
private:
    double a, b, delta;
    int M;
    Function &f;
public:
      Bisection(double _a, double _b, double _delta, int _M, Function &_f ):
            a(_a),
			b(_b),
			delta(_delta),
			M(_M),
			f(_f){};  

    double solve()
    {
        double h = b - a;
        double u = f(a);
        double c, w;
        for(int k = 1; k <= M; k++)
        {
            h = h/2;
            c = a + h;
            w = f(c);
            if (fabs(h) < delta || fabs(w) < eps)
	            break;
            else if (w*u >= 0)
	            a = c;
            else {;}				    
        }
        return c;
    }
};

class Newton : public EquationSolver
{
private:
    double x0, delta;
    int M;
    Function &f;
public:
      Newton(double _x0, double _delta, int _M, Function &_f ):
            x0(_x0),
			delta(_delta),
			M(_M),
			f(_f){};  

    double solve()
    {
        double x = x0;
        double u;
        for(int k = 0; k <= M; k++)
        {
            u = f(x);
            if (fabs(u) < eps)
	            break;
			x = x - u/f.diff(x);  
        }
        return x;
    }
};

class Secant : public EquationSolver
{
private:
    double x0, x1, delta;
    int M;
    Function &f;
public:
      Secant(double _x0, double _x1, double _delta, int _M, Function &_f):
            x0(_x0),
			x1(_x1),
			delta(_delta),
			M(_M),
			f(_f){};  

    double solve()
    {
        double m = x1;
        double n = x0;
        double u = f(m);
        double v = f(n);
        double temp;
        for(int k = 2; k <= M; k++)
        {
            if (fabs(u) > fabs(v))
            {
                temp = m; m = n; n = temp;
                temp = u; u = v; v = temp;
            }
            double s = (m - n)/(u - v);
            n = m;
            v = u;
            m = m - u*s;
            u = f(m);
            if (fabs(m-n) < delta || fabs(u) < eps)
                break;
        }
        return m;
    }
};