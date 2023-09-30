#include <iostream>
#include <eigen3/Eigen/Dense>

// DLS_solver类
class DLS_solver
{
public:
    DLS_solver() = default;
    virtual ~DLS_solver() = default;

    virtual Eigen::VectorXd solve() = 0;
};

// Normal_equations类
class Normal_equations : public DLS_solver
{
private:
    Eigen::MatrixXd X;
    Eigen::VectorXd y;
    Eigen::MatrixXd G;
public:
    Normal_equations() = default;

    // 构造时输入题目中所需的x和y
    Normal_equations(Eigen::MatrixXd& x, Eigen::MatrixXd& y)
    {
        // 初始化点集
        X.resize(x.rows(), 3);
        X.col(0) = Eigen::VectorXd::Ones(x.rows());
        X.col(1) = x;
        X.col(2) = x.array().square();
        this->y = y;
    }

    ~Normal_equations() override = default;

    // 求解正规方程组
    Eigen::VectorXd solve() override
    {
        // 计算Gram matrix G，并求解
        G = X.transpose() * X;
        Eigen::VectorXd c = G.ldlt().solve(X.transpose() * y);

        // 返回系数
        return c;
    }

    // 返回Gram matrix G
    Eigen::MatrixXd get_G() const
    {
        return G;
    }
};

// QR_factorization类
class QR_factorization : public DLS_solver
{
private:
    Eigen::MatrixXd X;
    Eigen::VectorXd y;
    Eigen::MatrixXd R;
public:
    QR_factorization() = default;

    // 构造时输入题目中所需的x和y
    QR_factorization(Eigen::MatrixXd& x, Eigen::MatrixXd& y)
    {
        // 初始化点集
        X.resize(x.rows(), 3);
        X.col(0) = Eigen::VectorXd::Ones(x.rows());
        X.col(1) = x;
        X.col(2) = x.array().square();
        this->y = y;
    }

    ~QR_factorization() override = default;

    // QR分解求解方程组
    Eigen::VectorXd solve() override
    {
        // QR分解求解
        Eigen::HouseholderQR<Eigen::MatrixXd> qr(X);
        Eigen::VectorXd c = qr.solve(y);

        // 存储得到的R
        R = qr.matrixQR().triangularView<Eigen::Upper>();

        // 返回系数
        return c;
    }

    // 返回矩阵 R1
    Eigen::MatrixXd get_R1() const 
    { 
        return R.block(0, 0, 3, 3); 
    }
};
