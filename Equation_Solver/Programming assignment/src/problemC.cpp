#include "equationsolver.h"

class Fun : public Function
{
public:
    double operator()(double _x)
    {
        double y;
        y = _x - tan(_x);
        return y;
    }
    double diff(double _x)
    {
        double y;
        y = 1 - 1/pow(cos(_x),2);
        return y;
    }
};

int main()
{
    Fun f;
    Newton n1(4.5, eps, 15, f);
    Newton n2(7.7, eps, 15, f);

    double ans1 = n1.solve();
    double ans2 = n2.solve();

    std::cout << "root near 4.5 is r1 =  " << ans1 << ", f(r) = " << (f)(ans1) << std::endl;
    std::cout << "root near 7.7 is r2 =  " << ans2 << ", f(r) = " << (f)(ans2) << std::endl;

    return 0;
}
