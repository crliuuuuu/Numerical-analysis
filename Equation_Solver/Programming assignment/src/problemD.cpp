#include "equationsolver.h"
#define pi acos(-1)

class Fun1 : public Function
{
public:
  double operator()(double _x)
    {
        double y;
        y = sin(0.5*_x) - 1;
        return y;
    }
};

class Fun2 : public Function
{
public:
  double operator()(double _x)
    {
        double y;
        y= exp(_x) - tan(_x);
        return y;
    }
};

class Fun3: public Function
{
public:
  double operator()(double _x)
    {
        double y;
        y= pow(_x,3) - 12*pow(_x,2) + 3*_x + 1;
        return y;
    }
};


int main()
{
    Fun1 f1;
    Fun2 f2;
    Fun3 f3;
    Secant s1(0, pi/2, eps, 50, f1);
    Secant s12(0, pi*4.5, eps, 50, f1);
    Secant s2(1, 1.4, eps, 50, f2);
    Secant s22(1, 5.4, eps, 50, f2);
    Secant s3(0, -0.5, eps, 50, f3);
    Secant s32(0.35, 0.45, eps, 50, f3);

    double ans1 = s1.solve();
    double ans12 = s12.solve();
    double ans2 = s2.solve();
    double ans22 = s22.solve();
    double ans3 = s3.solve();
    double ans32 = s32.solve();


    std::cout << "root of function 1 is r1 =  " << ans1 << ", f(r) = " << (f1)(ans1) << std::endl;
    std::cout << "root of function 1 with different initial values is r1 =  " << ans12 << ", f(r) = " << (f1)(ans12) << std::endl;
    std::cout << "root of function 2 is r2 =  " << ans2 << ", f(r) = " << (f2)(ans2) << std::endl;
    std::cout << "root of function 2 with different initial values is r2 =  " << ans22 << ", f(r) = " << (f2)(ans22) << std::endl;
    std::cout << "root of function 3 is r3 =  " << ans3 << ", f(r) = " << (f3)(ans3) << std::endl;
    std::cout << "root of function 3 with different initial values is r3 =  " << ans32 << ", f(r) = " << (f3)(ans32) << std::endl;

    return 0;
}
