#include "interpolationsolver.h"
#include <iostream>
#include <fstream>

int main()
{
    std::vector<double> T{0,0,3,3,5,5,8,8,13,13};
    std::vector<double> D_and_S{0,75,225,77,383,80,623,74,993,72};
    Hermite_Formula h(T, D_and_S);

    std::cout << "The car's distance at t = 10s is " << h.solve(10) << " feet" << std::endl;
    std::cout << "The car's speed at t = 10s is " << h.solvediff(10) << " feet/s" << std::endl;

    /*
    Determine the maximum speed in the observed process, sample the interpolation function by plotting points, 
    and store the sampled data in the probD_result.txt file for importing and plotting in Python.
    */
    std::ofstream fout;
    fout.open("probD_result.txt");
    fout << "Sampling the derivative of Hermite polynomial with interval 0.01 in [0,13]" << std::endl;
    double max = 0;
    for (int i = 0; i < 100*13; i++)
    {
        fout << h.solvediff(i*0.01) << ", ";
        if (h.solvediff(i*0.01) > max)
            max = h.solvediff(i*0.01);
    }
    fout << h.solvediff(13.0) << std::endl;
    std::cout << "The max speed during the observation is " << max << " feet per second" << std::endl;


    fout.close();
    return 0;
}
