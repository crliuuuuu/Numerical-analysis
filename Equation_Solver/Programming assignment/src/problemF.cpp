#include "equationsolver.h"
#define pi acos(-1)
const double l = 89, h = 49, beta1 = 11.5/180*pi;

class Fun1 : public Function
{
public:
    double operator()(double _x)
    {
        double y;
        double D = 55;
        double A = l*sin(beta1);
        double B = l*cos(beta1);
        double C = (h + 0.5*D)*sin(beta1) - 0.5*D*tan(beta1);
        double E = (h + 0.5*D)*cos(beta1) - 0.5*D;
        y = A*sin(_x)*cos(_x) + B*pow(sin(_x),2) - C*cos(_x) - E*sin(_x);
        return y;
    }
    double diff(double _x)
    {
        double y;
        double D = 55;
        double A = l*sin(beta1);
        double B = l*cos(beta1);
        double C = (h + 0.5*D)*sin(beta1) - 0.5*D*tan(beta1);
        double E = (h + 0.5*D)*cos(beta1) - 0.5*D;
        y = -A*pow(sin(_x),2) + A*pow(cos(_x),2) + 2*B*sin(_x)*cos(_x) + C*sin(_x) - E*cos(_x);
        return y;
    }
};

class Fun2 : public Function
{
public:
    double operator()(double _x)
    {
        double y;
        double D = 30;
        double A = l*sin(beta1);
        double B = l*cos(beta1);
        double C = (h + 0.5*D)*sin(beta1) - 0.5*D*tan(beta1);
        double E = (h + 0.5*D)*cos(beta1) - 0.5*D;
        y = A*sin(_x)*cos(_x) + B*pow(sin(_x),2) - C*cos(_x) - E*sin(_x);
        return y;
    }
    double diff(double _x)
    {
        double y;
        double D = 30;
        double A = l*sin(beta1);
        double B = l*cos(beta1);
        double C = (h + 0.5*D)*sin(beta1) - 0.5*D*tan(beta1);
        double E = (h + 0.5*D)*cos(beta1) - 0.5*D;
        y = -A*pow(sin(_x),2) + A*pow(cos(_x),2) + 2*B*sin(_x)*cos(_x) + C*sin(_x) - E*cos(_x);
        return y;
    }
};

int main()
{
    Fun1 f1;
    Newton n1(33*pi/180, eps, 20, f1);
    double ans1 = n1.solve();
    std::cout << "root using Newton Method in (a) is alpha = " << ans1*180/pi << "°, f(alpha) = " << (f1)(ans1) << std::endl;

    Fun2 f2;
    Newton n2(33*pi/180, eps, 20, f2);
    double ans2 = n2.solve();
    std::cout << "root using Newton Method in (b) is alpha = " << ans2*180/pi << "°, f(alpha) = " << (f2)(ans2) << std::endl;

    Secant s1(32*pi/180, 33*pi/180, eps, 20, f1);
    double ans3 = s1.solve();
    std::cout << "root using Secand Method in (c) is alpha = " << ans3*180/pi << "°, f(alpha) = " << (f1)(ans3) << std::endl;   

    Secant s2(160*pi/180, 180*pi/180, eps, 20, f1);
    double ans4 = s2.solve();
    std::cout << "root using Secand Method with different initial value in (c) is alpha = " << ans4*180/pi << "°, f(alpha) = " << (f1)(ans4) << std::endl;   

    return 0;
}
