#include "interpolationsolver.h"
#include <iostream>
#include <fstream>
#define pi acos(-1)

class Fun : public Function
{
public:
    double operator()(double _x)
    {
        double y;
        y = 1/(1 + 25.0*_x*_x);
        return y;
    }
};


int main()
{
    std::ofstream fout;
    fout.open("probC_result.txt");
    Fun f;

    /*
        分别对不同n的值对应的插值函数进行描点采样，并将采样数据存储至probC_result.txt文件以便导入至Python画图
    */
    for (int k = 1; k < 5; k++)
    {
        std::vector<double> v;
        int n = 5*k;
        for (int i = 1; i <= n ; i++)
        {
            v.push_back(cos(pi*(2*i - 1)/(2*n)));
        }
        Newton_Formula b(f, v);
        fout << "n = " << n << ", sampling with interval 0.002 in [-1,1]" << std::endl;
        for (int i = 0; i < 1000; i++)
        {
            fout << b.solve(-1.0 + 2.0*i/1000) << ", ";
        }
        fout << b.solve(1.0) <<std::endl;
    }

    fout.close();
    return 0;
}
