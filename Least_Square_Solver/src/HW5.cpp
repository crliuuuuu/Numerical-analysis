#include <iostream>
#include <eigen3/Eigen/Dense>
#include "DLS_solver.h"

int main()
{
    Eigen::MatrixXd x = (Eigen::MatrixXd(21, 1) <<
        0,0.5,1,1.5,2,2.5,3,3.5,4,4.5,5,5.5,6,6.5,7,7.5,8,8.5,9,9.5,10).finished();
    Eigen::MatrixXd y = (Eigen::MatrixXd(21, 1) <<
    2.9,2.7,4.8,5.3,7.1,7.6,7.7,7.6,9.4,9.0,9.6,10.0,10.2,9.7,8.3,8.4,9.0,8.3,6.6,6.7,4.1).finished();

    // 使用正规方程组求解
    auto solver1 = Normal_equations(x, y);

    Eigen::VectorXd a_N = solver1.solve();

    // 输出系数
    std::cout << "The results using normal equations are: " << std::endl;
    std::cout << "a0 = " << a_N[0] << std::endl;
    std::cout << "a1 = " << a_N[1] << std::endl;
    std::cout << "a2 = " << a_N[2] << std::endl;

    // 使用QR分解求解
    auto solver2 = QR_factorization(x, y);

    Eigen::VectorXd a_QR = solver2.solve();

    // 输出系数
    std::cout << "The results using QR factorization are: " << std::endl;
    std::cout << "a0 = " << a_QR[0] << std::endl;
    std::cout << "a1 = " << a_QR[1] << std::endl;
    std::cout << "a2 = " << a_QR[2] << std::endl;

    // 输出G的条件数
    Eigen::MatrixXd G = solver1.get_G();
    std::cout << "The condition number of G is: " << G.lpNorm<2>() * G.inverse().lpNorm<2>() << std::endl;

    // 输出R1的条件数
    Eigen::MatrixXd R1 = solver2.get_R1();
    std::cout << "The condition number of R1 is: " << R1.lpNorm<2>() * R1.inverse().lpNorm<2>() << std::endl;

    return 0;
}




