#include <iostream>
#include <cmath>
#include <algorithm>
#include <limits>
#include <vector>

double eps = std::numeric_limits<float>::epsilon();

/*
    定义函数类，对不同函数进行运算符()的重载
*/
class Function
{
public:
    virtual double operator()(double _x) = 0;
};

/*
    定义插值类，抽象Newton插值和Hermite插值
*/
class Interpolation
{
public:
    virtual double solve(double _x) = 0; 
};

/*
    Newton插值多项式继承函数类，给定系数a_k即可唯一表示对应的函数
*/
class Newton_polynomial : public Function
{
private:
    std::vector<double> C, P;
    int n;
public:
    Newton_polynomial(std::vector<double> _C, std::vector<double> _P): C(_C), P(_P) {}
    double operator()(double _x)
    {
        double y = C[0];
        n = P.size();
        for (int i = 1; i < n; i++)
        {
            double delta_y = C[i];
            for (int j = 0; j < i; j++)
            {
                delta_y = delta_y*(_x - P[j]);
            }
            y = y + delta_y;
        }
        return y;
    }
};

/*
    Hermite插值多项式继承函数类，给定系数a_k即可唯一表示对应的函数
*/
class Hermite_polynomial : public Function
{
private:
    std::vector<double> C, L;
    int n;
public:
    Hermite_polynomial(std::vector<double> _C, std::vector<double> _L): C(_C), L(_L) {}
    double operator()(double _x)
    {
        double y = C[0];
        n = L.size();
        for (int i = 1; i < n; i++)
        {
            double delta_y = C[i];
            for (int j = 0; j < i; j++)
            {
                delta_y = delta_y*(_x - L[j]);
            }
            y = y + delta_y;
        }
        return y;
    }
};

/*
    利用Newton插值公式进行插值，继承插值类
*/
class Newton_Formula : public Interpolation
{
private:
    std::vector<double> X, Coef;
    Function &f;
    int n;
public:
    //构造函数，输入待插值函数_f和插值点_X后,利用差商求出Newton插值多项式系数
    Newton_Formula(Function &_f, std::vector<double> _X): f(_f), X(_X) 
    {
        n = X.size();
        double table[n][n];
        for (int i = 0; i < n; i++)
        {
            table[i][0] = f(X[i]);
        }
        for (int i = 1; i < n; i++)
        {
            for (int j = i; j < n; j++)
            {
                table[j][i] = (table[j][i-1] - table[j-1][i-1])/(X[j] - X[j-i]);
            }
        }
        for (int i = 0; i < n; i++)
        {
            Coef.push_back(table[i][i]);
        }
    }
    //根据题意，solve()函数返回Newton插值多项式任意一点对应的多项式值
    double solve(double _x)
    {
        Newton_polynomial Pn(Coef, X);
        return Pn(_x);
    }
    /*
        如有需要，Get_coef()函数返回Newton插值多项式的相应系数
        注意：根据题意，希望能将系数封装在类内，因此除了实验报告展示之外，应尽量避免使用该函数
    */
    std::vector<double> Get_coef()
    {
        return Coef;
    }
};

/*
    利用Hermite插值公式进行插值，继承插值类
*/
class Hermite_Formula : public Interpolation
{
private:
    std::vector<double> X, Label, Coef;
    int n;
public:
    /*
        构造函数，输入待插值点_Label以及每一个插值点对应的函数（或其导数值）_X后,利用差商求出Hermite插值多项式系数
        注意，待插值点在_Label和_X中的输入排序需要满足：1.相同的插值点需要相邻；2.相同的插值点对应的函数或其导数值按导数阶数升序相邻排列。
    */
    Hermite_Formula(std::vector<double> _Label, std::vector<double> _X): Label(_Label), X(_X)
    {
        n = X.size();
        double table[n][n] = {{0}};
        int count = 1;
        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < count; j++)
            {
                table[i][j] = X[i+j+1-count];
            }
            if (Label[i] == Label[i+1])
                count = count + 1;
            else
                count = 1;
        }
        for (int i = 1; i < n; i++)
        {
            for (int j = i; j < n; j++)
            {
                if (table[j][i] == 0)
                    table[j][i] = (table[j][i-1] - table[j-1][i-1])/(Label[j] - Label[j-i]);
            }
        }
        for (int i = 0; i < n; i++)
        {
            Coef.push_back(table[i][i]);
        }
    }
    //根据题意，solve()函数返回Hermite插值多项式任意一点对应的多项式值
    double solve(double _x)
    {
        Newton_polynomial Pn(Coef, Label);
        return Pn(_x);
    }
    //由于题目要求，solvediff()函数返回插值多项式在_x点处的导数
    double solvediff(double _x)
    {
        Newton_polynomial Pn(Coef, Label);
        return (Pn(_x + eps) - Pn(_x - eps))/(2.0*eps);
    }
    /*
        如有需要，Get_coef()函数返回Hermite插值多项式的相应系数
        注意：根据题意，希望能将系数封装在类内，因此除了实验报告展示之外，应尽量避免使用该函数
    */
    std::vector<double> Get_coef()
    {
        return Coef;
    }
};
