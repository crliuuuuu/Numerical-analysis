#include <iostream>
#include <cmath>
#include <vector>
#include <fstream>


int main()
{
    std::ofstream fout;
    fout.open("point_value.txt");

    std::vector<double> v1,v2,v3;
    for (int i = 0; i <= 100; i++)
    {
        double x = 0.99 + i*0.02/100;
        double y1, y2, y3;
        y1 = x*x*x*x*x*x*x*x - 8*x*x*x*x*x*x*x + 28*x*x*x*x*x*x - 56*x*x*x*x*x + 70*x*x*x*x - 56*x*x*x + 28*x*x -8*x + 1;
        y2 = (((((((x - 8)*x + 28)*x - 56)*x + 70)*x - 56)*x + 28)*x - 8)*x + 1;
        y3 = (x - 1)*(x - 1)*(x - 1)*(x - 1)*(x - 1)*(x - 1)*(x - 1)*(x - 1);
        v1.push_back(y1);
        v2.push_back(y2);
        v3.push_back(y3);
    }
    fout << "Point values of function a: " << std::endl;
    for (int i = 0; i < 100; i++)
    {
        fout << v1[i] << ", ";
    }
    fout << v1[100] <<std::endl;
    fout << "Point values of function b: " << std::endl;
    for (int i = 0; i < 100; i++)
    {
        fout << v2[i] << ", ";
    }
    fout << v2[100] <<std::endl;
    fout << "Point values of function c: " << std::endl;
    for (int i = 0; i < 100; i++)
    {
        fout << v3[i] << ", ";
    }
    fout << v3[100] <<std::endl;
    fout.close();

    std::cout<< "UFL = " << pow(2,-1) << std::endl;
    std::cout<< "OFL = " << pow(2,1)*(2-pow(2,1-3))<< std::endl;
    std::vector<double> F;
    for (int i = -1; i <= 1; i++)
    {
        for (int j = 0; j < 4; j++)
        {
            F.push_back((1+j*0.25)*pow(2,i));
            F.push_back(-(1+j*0.25)*pow(2,i));
        }
    }
    F.push_back(0);
    std::cout << "All numbers in F are: " << std::endl;
    for (int i = 1; i < F.size(); i++)
    {
        std::cout << F[i-1] << ", ";
    }
    std::cout << F[F.size()-1] << ". " << std::endl;
    std::cout << F.size() <<" numbers in total." << std::endl;

    std::vector<double> S;
    for (int j = 1; j < 4; j++)
    {
        S.push_back((j*0.25)*pow(2,-1));
        S.push_back(-(j*0.25)*pow(2,-1));
    }
    std::cout << "All subnormal numbers are: " << std::endl;
    for (int i = 1; i < S.size(); i++)
    {
        std::cout << S[i-1] << ", ";
    }
    std::cout << S[S.size()-1] << ". " << std::endl;
    std::cout << S.size() <<" numbers in total." << std::endl;
    return 0;
}
