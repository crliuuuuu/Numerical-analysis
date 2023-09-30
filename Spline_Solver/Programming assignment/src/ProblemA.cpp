#include <iostream>
#include <fstream>
#include "spline.h"
#include <cmath>
#include <algorithm>
#include <vector>
#include <jsoncpp/json/json.h>

#define Spline_Form Bspline_interpolation 

class Fun : public Function
{
public:
    double operator()(double _x)
    {
        double y;
        y = 1.0/(1.0 + 25.0*_x*_x);
        return y;
    }
    double diff(double _x) 
    {
        double y;
        y = -50.0*_x/pow(1 + 25.0*_x*_x, 2);
        return y;  
    }
    double diff_2(double _x)
    {
        double y;
        y = 50.0*(75.0*_x*_x - 1)/pow(1.0 + 25.0*_x*_x, 3);
        return y;  
    }
};

int main()
{
    Json::Reader reader;
	Json::Value root;
    std::ofstream fout, fout1;
 
	std::ifstream in("ProblemA.json", std::ios::binary);
	if (!in.is_open())
	{
		std::cout << "error: cannot open file." << std::endl;
		return -1;
    }

    std::vector<int> point_Num;
    double left, right;
    int mth, ord;
    Fun f;

    if (reader.parse(in, root))
	{
        left = root["interval_left"].asDouble();
        right = root["interval_right"].asDouble();
        if (left > right)
        {
            std::cerr<< "error: invalid interval input." <<std::endl;
		    return -1;
        }
        for (int i = 0; i < root["n"].size(); i++)
	    {
	        int temp_n = root["n"][i].asInt();
	        point_Num.push_back(temp_n);
	    }
        mth = root["method"].asInt();
        ord = root["order"].asInt();
    }

    fout.open("ProblemA_result_points.txt");
    fout1.open("ProblemA_result_error.txt");
    fout << "order = " << ord << ", method = " << mth << std::endl;
    fout1 << "max interpolation error with order = " << ord << ", method = " << mth << std::endl;

    for (int i = 0; i < point_Num.size(); i++)
	{
        std::vector<double> T;
        for (int j = 0; j < point_Num[i]; j++)
		{
            T.push_back(left + j*(right-left)/(point_Num[i]-1));
		}

        Spline_Form P(f, T, mth, ord);
        P.solve();
        fout << "N = " << point_Num[i] << ", sampling with interval 0.01 in [-1,1]" << std::endl;
        for (int j = 0; j < 200; j++)
		{
            fout << P(left + j*0.01) << ", ";
		}
        fout << P(right) <<std::endl;

        fout1 << "N = " << point_Num[i] << ": ";
        double max_error = 0;
        for (int j = 0; j < point_Num[i] - 1; j++)
		{
            if (fabs(P(0.5*(T[j] + T[j+1])) - f(0.5*(T[j] + T[j+1]))) > max_error)
                max_error = fabs(P(0.5*(T[j] + T[j+1])) - f(0.5*(T[j] + T[j+1])));
		}
        fout1 << max_error << std::endl;

	}
    fout.close();
    fout1.close();

    return 0;
}

