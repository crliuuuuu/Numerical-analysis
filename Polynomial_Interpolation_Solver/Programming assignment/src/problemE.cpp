#include "interpolationsolver.h"
#include <iostream>
#include <fstream>

class Fun1 : public Function
{
public:
    double operator()(double _x)
    {
        double y;
        switch ((int)_x)
        {
            case 0:
                y = 6.67;
                break;
            case 6:
                y = 17.3;
                break;
            case 10:
                y = 42.7;
                break;
            case 13:
                y = 37.3;
                break;
            case 17:
                y = 30.1;
                break;
            case 20:
                y = 29.3;
                break;
            case 28:
                y = 28.7;
                break;
            default:
                break;
        }
        return y;
    }
};

class Fun2 : public Function
{
public:
    double operator()(double _x)
    {
        double y;
        switch ((int)_x)
        {
            case 0:
                y = 6.67;
                break;
            case 6:
                y = 16.1;
                break;
            case 10:
                y = 18.9;
                break;
            case 13:
                y = 15.0;
                break;
            case 17:
                y = 10.6;
                break;
            case 20:
                y = 9.44;
                break;
            case 28:
                y = 8.89;
                break;
            default:
                break;
        }
        return y;
    }
};

int main()
{
    std::ofstream fout;
    fout.open("probE_result.txt");
    Fun1 Sp1;
    Fun2 Sp2;
    
    std::vector<double> v1{0,6,10,13,17,20,28};
    std::vector<double> v2{0,6,10,13,17,20,28};

    Newton_Formula b1(Sp1, v1);
    Newton_Formula b2(Sp2, v2);
    std::vector<double> Coef1,Coef2;

    /*
        为方便报告展示，记录两类样本对应的插值函数系数和15天后预测的Average Weight值
    */
    Coef1 = b1.Get_coef();
    Coef2 = b2.Get_coef();
    std::cout << "The Coefficients of curve for Sp1 are: ";
    for (int i = 0; i < v1.size(); i++)
    {
        std::cout << Coef1[i] << " ";
    }
    std::cout << "" << std::endl;
    std::cout << "The Coefficients of curve for Sp2 are: ";
    for (int i = 0; i < v1.size(); i++)
    {
        std::cout << Coef2[i] << " ";
    }
    std::cout << "" << std::endl;
    std::cout << "Average weight of Sp1 after another 15 days is: " << b1.solve(43.0) << std::endl;
    std::cout << "Average weight of Sp2 after another 15 days is: " << b2.solve(43.0) << std::endl;

    /*
        分别对插值函数进行描点采样，并将采样数据存储至probE_result.txt文件以便导入至Python画图
    */
    fout << "Sampling average weight of Sp1 with interval 0.01 in [0,28]" << std::endl;
    for (int i = 0; i < 2800; i++)
    {
        fout << b1.solve(0.01*i) << ", ";
    }
        fout << b1.solve(28.0) <<std::endl;

    fout << "Sampling average weight of Sp1 with interval 0.01 in [28,43]" << std::endl;
    for (int i = 0; i < 1500; i++)
    {
        fout << b1.solve(28 + 0.01*i) << ", ";
    }
        fout << b1.solve(43.0) <<std::endl;

    fout << "Sampling average weight of Sp2 with interval 0.01 in [0,28]" << std::endl;
    for (int i = 0; i < 2800; i++)
    {
        fout << b2.solve(0.01*i) << ", ";
    }
        fout << b2.solve(28.0) <<std::endl;

    fout << "Sampling average weight of Sp2 with interval 0.01 in [28,43]" << std::endl;
    for (int i = 0; i < 1500; i++)
    {
        fout << b2.solve(28 + 0.01*i) << ", ";
    }
        fout << b2.solve(43.0) <<std::endl;

    fout.close();
    return 0;
}
