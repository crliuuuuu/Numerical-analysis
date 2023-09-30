#include "interpolationsolver.h"
#include <iostream>
#include <fstream>

class Fun : public Function
{
public:
    double operator()(double _x)
    {
        double y;
        y = 1/(1 + _x*_x);
        return y;
    }
};


int main()
{
    std::ofstream fout;
    fout.open("probB_result.txt");
    Fun f;

    for (int k = 1; k < 5; k++)
    {
        std::vector<double> v;
        int n = 2*k;
        v.push_back(0.0);
        for (int i = 1; i <= k; i++)
        {
            v.push_back(10.0*i/n);
            v.push_back(-10.0*i/n);
        }
        Newton_Formula b(f, v);
        fout << "n = " << n << ", sampling with interval 0.01 in [-5,5]" << std::endl;
        for (int i = 0; i < 1000; i++)
        {
            fout << b.solve(-5.0 + i*0.01) << ", ";
        }
        fout << b.solve(5.0) <<std::endl;
    }

    fout.close();
    return 0;
}
