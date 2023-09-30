#include "equationsolver.h"
#define pi acos(-1)
const double L = 10, r = 1, V = 12.4;

class Fun : public Function
{
public:
    double operator()(double _x)
    {
        double y;
        y = V - L*(0.5*pi*r*r - r*r*asin(_x/r) - _x*sqrt(r*r - _x*_x));
        return y;
    }
    double diff(double _x)
    {
        double y;
        y = -L*(-r/sqrt(1 - (_x/r)*(_x/r)) - sqrt(r*r - _x*_x ) + _x*_x/sqrt(r*r - _x*_x));
        return y;
    }
};

int main()
{
    Fun f;
    Bisection b(0, 0.2*r, 0.01, 50, f);
    Newton n(0.2*r, 0.01, 10, f);
    Secant s(0, 0.2*r, 0.01, 10, f);

    double ans1 = b.solve();
    double ans2 = n.solve();
    double ans3 = s.solve();

    std::cout << "root using Bisection Method is h =  " << ans1 << ", f(h) = " << (f)(ans1) << " and the depth is: " << r - ans1 << std::endl;
    std::cout << "root using Newton Method is h =  " << ans2 << ", f(h) = " << (f)(ans2) << " and the depth is: " << r - ans2 <<std::endl;
    std::cout << "root using Secant Method is h =  " << ans3 << ", f(h) = " << (f)(ans3) << " and the depth is: " << r - ans3 << std::endl;

    return 0;
}