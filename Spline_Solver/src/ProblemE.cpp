#include <iostream>
#include <fstream>
#include "spline.h"
#include <cmath>
#include <algorithm>
#include <vector>
#include <jsoncpp/json/json.h>


int main()
{
    Json::Reader reader;
	Json::Value root;
    std::ofstream fout;
 
	std::ifstream in("ProblemE.json", std::ios::binary);
	if (!in.is_open())
	{
		std::cout << "error: cannot open file." << std::endl;
		return -1;
    }

    std::vector<int> point_Num;
    std::string Form;
    int ord, mth;

    if (reader.parse(in, root))
	{
        for (int i = 0; i < root["n"].size(); i++)
	    {
	        int temp_n = root["n"][i].asInt();
	        point_Num.push_back(temp_n);
	    }
        ord = root["order"].asInt();
        Form = root["spline_form"].asString();
    }
    if (Form == "Bspline")
        mth = 1;
    else if (Form == "PP")
        mth = 2;
    else
    {
        std::cout << "error: invalid spline form." << std::endl;
		return -1;
    }
 

    fout.open("ProblemE_result_points.txt");
    fout << "order = " << ord << ", method = " << Form << std::endl;

    for (int i = 0; i < point_Num.size(); i++)
	{
        std::vector<std::vector<double> > T;
        std::vector<double> Tc;
        int num_1 = floor((point_Num[i]- 4.0)/4.0);
        int num_2 = floor((point_Num[i]- 4.0)/2.0) - num_1;
        Tc.push_back(0.0);
        Tc.push_back(sqrt(3.0)*2.0/3.0);
        T.push_back(Tc);
        Tc.clear();
        for (int j = 1; j <= num_1; j++)
		{
            std::vector<double> temp_T;
            double x = j*sqrt(3.0)/(num_1+1.0);
            temp_T.push_back(x);
            temp_T.push_back((sqrt(3.0-x*x)+sqrt(fabs(x)))*2.0/3.0);
            T.push_back(temp_T);
		}
        Tc.push_back(sqrt(3.0));
        Tc.push_back(sqrt(sqrt(3.0))*2.0/3.0);
        T.push_back(Tc);
        Tc.clear();
        for (int j = 1; j <= num_2; j++)
		{
            std::vector<double> temp_T;
            double x = sqrt(3.0) - j*sqrt(3.0)/(num_2+1.0);
            temp_T.push_back(x);
            temp_T.push_back((-sqrt(3.0-x*x)+sqrt(fabs(x)))*2.0/3.0);
            T.push_back(temp_T);
		}
        Tc.push_back(0.0);
        Tc.push_back(-sqrt(3.0)*2.0/3.0);
        T.push_back(Tc);
        Tc.clear();
        for (int j = 1; j <= num_2; j++)
		{
            std::vector<double> temp_T;
            double x = -j*sqrt(3.0)/(num_2+1.0);
            temp_T.push_back(x);
            temp_T.push_back((-sqrt(3.0-x*x)+sqrt(fabs(x)))*2.0/3.0);
            T.push_back(temp_T);
		}
        Tc.push_back(-sqrt(3.0));
        Tc.push_back(sqrt(sqrt(3.0))*2.0/3.0);
        T.push_back(Tc);
        Tc.clear();
        for (int j = 1; j <= num_1; j++)
		{
            std::vector<double> temp_T;
            double x = -sqrt(3.0) + j*sqrt(3.0)/(num_1+1.0);
            temp_T.push_back(x);
            temp_T.push_back((sqrt(3.0-x*x)+sqrt(fabs(x)))*2.0/3.0);
            T.push_back(temp_T);
		}
        Tc.push_back(0.0);
        Tc.push_back(sqrt(3.0)*2.0/3.0);
        T.push_back(Tc);
        Tc.clear();

        curve_spline C(T, mth, ord);
        C.solve();
        std::vector<std::vector<double> > FitP = C.Get_Point();
    
        fout << "n = " << point_Num[i] << std::endl;
        for (int j = 0; j < FitP.size(); j++)
		{
            fout << "component " << j+1 << std::endl;
            for (int k = 0; k < FitP[j].size() - 1; k++)
		    {
                fout << FitP[j][k] << ", ";
		    }
            fout << FitP[j][FitP[j].size() - 1] <<std::endl;
		}
	}
    fout.close();

    return 0;
}

