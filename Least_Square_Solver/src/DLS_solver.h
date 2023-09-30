#include <iostream>
#include <eigen3/Eigen/Dense>


class DLS_solver
{
public:
    DLS_solver() = default;
    virtual ~DLS_solver() = default;

    virtual Eigen::VectorXd solve() = 0;
};


class Normal_equations : public DLS_solver
{
private:
    Eigen::MatrixXd X;
    Eigen::VectorXd y;
    Eigen::MatrixXd G;
public:
    Normal_equations() = default;

    Normal_equations(Eigen::MatrixXd& x, Eigen::MatrixXd& y)
    {
        // initialize
        X.resize(x.rows(), 3);
        X.col(0) = Eigen::VectorXd::Ones(x.rows());
        X.col(1) = x;
        X.col(2) = x.array().square();
        this->y = y;
    }

    ~Normal_equations() override = default;

    // solve normal equations
    Eigen::VectorXd solve() override
    {
        // calculate Gram matrix G
        G = X.transpose() * X;
        Eigen::VectorXd c = G.ldlt().solve(X.transpose() * y);

        // return coefficient
        return c;
    }

    // return Gram matrix G
    Eigen::MatrixXd get_G() const
    {
        return G;
    }
};


class QR_factorization : public DLS_solver
{
private:
    Eigen::MatrixXd X;
    Eigen::VectorXd y;
    Eigen::MatrixXd R;
public:
    QR_factorization() = default;

    QR_factorization(Eigen::MatrixXd& x, Eigen::MatrixXd& y)
    {
        // initialize
        X.resize(x.rows(), 3);
        X.col(0) = Eigen::VectorXd::Ones(x.rows());
        X.col(1) = x;
        X.col(2) = x.array().square();
        this->y = y;
    }

    ~QR_factorization() override = default;

    // QRfactorization
    Eigen::VectorXd solve() override
    {
        Eigen::HouseholderQR<Eigen::MatrixXd> qr(X);
        Eigen::VectorXd c = qr.solve(y);

        R = qr.matrixQR().triangularView<Eigen::Upper>();

        // return coefficient
        return c;
    }

    // return R1
    Eigen::MatrixXd get_R1() const 
    { 
        return R.block(0, 0, 3, 3); 
    }
};
