#include <iostream>
#include <fstream>
#include "spline.h"
#include <cmath>
#include <algorithm>
#include <vector>
#include <jsoncpp/json/json.h>


class Fun : public Function
{
public:
    double operator()(double _x)
    {
        double y;
        y = 1/(1 + _x*_x);
        return y;
    }
    double diff(double _x) 
    {
        double y;
        y = -2.0*_x/pow(1.0 + _x*_x, 2);
        return y;  
    }
    double diff_2(double _x)
    {
        double y;
        y = 2.0*(3.0*_x*_x - 1.0)/pow(1.0 + _x*_x, 3);
        return y;  
    }
};

int main()
{
    Json::Reader reader;
	Json::Value root;
    std::ofstream fout, fout1, fout2;
 
	std::ifstream in("ProblemCD.json", std::ios::binary);
	if (!in.is_open())
	{
		std::cout << "error: cannot open file." << std::endl;
		return -1;
    }

    std::vector<double> point_1, point_2, output_point;
    Fun f;

    if (reader.parse(in, root))
	{
        for (int i = 0; i < root["point_1"].size(); i++)
	    {
	        double temp_n = root["point_1"][i].asDouble();
	        point_1.push_back(temp_n);
	    }
        for (int i = 0; i < root["point_2"].size(); i++)
	    {
	        double temp_n = root["point_2"][i].asDouble();
	        point_2.push_back(temp_n);
	    }
        for (int i = 0; i < root["output_point"].size(); i++)
	    {
	        double temp_n = root["output_point"][i].asDouble();
	        output_point.push_back(temp_n);
	    }
    }

    fout.open("ProblemC_result_points.txt");
    fout1.open("ProblemD_result_error.txt");
    fout2.open("ProblemD_result_errorpoints.txt");

    Bspline_interpolation B1(f, point_1, 1, 3);
    B1.solve();
    fout << "Cubic spline sampling with interval 0.01 in [-5,5]" << std::endl;
    fout2 << "Error function of cubic spline sampling with interval 0.01 in [-5,5]" << std::endl;
    for (int i = 0; i < 1000; i++)
	{
        fout << B1(point_1[0] + i*0.01) << ", ";
        fout2 << fabs(B1(point_1[0] + i*0.01) - f(point_1[0] + i*0.01)) << ", ";
	}
    fout << B1(point_1[point_1.size()-1]) <<std::endl;
    fout2 << fabs(B1(point_1[point_1.size()-1]) - f(point_1[point_1.size()-1])) <<std::endl;

    Bspline_interpolation B2(f, point_2, 1, 1);
    B2.solve();
    fout << "Linear spline sampling with interval 0.01 in [-4.5,4.5]" << std::endl;
    fout2 << "Error function of linear spline sampling with interval 0.01 in [-4.5,4.5]" << std::endl;
    for (int i = 0; i < 900; i++)
	{
        fout << B2(point_2[0] + i*0.01) << ", ";
        fout2 << fabs(B2(point_2[0] + i*0.01) - f(point_2[0] + i*0.01)) << ", ";
	}
    fout << B2(point_2[point_2.size()-1]) <<std::endl;
    fout2 << fabs(B2(point_2[point_2.size()-1]) - f(point_2[point_2.size()-1])) <<std::endl;

    std::vector<double> err_1, err_2;
    fout1 << "Error values of cubic spline at required points" << std::endl;
    for (int i = 0; i < output_point.size() - 1; i++)
	{
	        fout1 << fabs(B1(output_point[i]) - f(output_point[i])) << ", ";
	}
    fout1 << fabs(B1(output_point[output_point.size() - 1]) - f(output_point[output_point.size() - 1])) << std::endl;
    fout1 << "Error values of linear spline at required points" << std::endl;
    for (int i = 0; i < output_point.size() - 1; i++)
	{
	        fout1 << fabs(B2(output_point[i]) - f(output_point[i])) << ", ";
	}
    fout1 << fabs(B2(output_point[output_point.size() - 1]) - f(output_point[output_point.size() - 1])) << std::endl;
    
    fout.close();
    fout1.close();
    fout2.close();

    return 0;
}

