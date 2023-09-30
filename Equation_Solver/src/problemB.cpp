#include "equationsolver.h"
#define pi acos(-1)

class Fun1 : public Function
{
public:
  double operator()(double _x)
    {
        double y;
        y = 1/_x - tan(_x);
        return y;
    }
};

class Fun2 : public Function
{
public:
  double operator()(double _x)
    {
        double y;
        y= 1/_x - pow(2,_x);
        return y;
    }
};

class Fun3: public Function
{
public:
  double operator()(double _x)
    {
        double y;
        y= pow(2,-_x) + exp(_x) + 2*cos(_x) - 6;
        return y;
    }
};

class Fun4: public Function
{
public:
  double operator()(double _x)
    {
        double y;
        y= (pow(_x,3) + 4*pow(_x,2) + 3*_x + 5)/(2*pow(_x,3) - 9*pow(_x,2) + 18*_x - 2);
        return y;
    }
};

int main()
{
    Fun1 f1;
    Fun2 f2;
    Fun3 f3;
    Fun4 f4;
    Bisection b1(0, pi/2, eps, 100, f1);
    Bisection b2(0, 1, eps, 100, f2);
    Bisection b3(1, 3, eps, 100, f3);
    Bisection b4(0, 4, eps, 100, f4);

    double ans1 = b1.solve();
    double ans2 = b2.solve();
    double ans3 = b3.solve();
    double ans4 = b4.solve();

    std::cout << "root of function 1 is r1 =  " << ans1 << ", f(r) = " << (f1)(ans1) << std::endl;
    std::cout << "root of function 2 is r2 =  " << ans2 << ", f(r) = " << (f2)(ans2) << std::endl;
    std::cout << "root of function 3 is r3 =  " << ans3 << ", f(r) = " << (f3)(ans3) << std::endl;
    std::cout << "root of function 4 is r4 =  " << ans4 << ", f(r) = " << (f4)(ans4) << std::endl;

    return 0;
}
