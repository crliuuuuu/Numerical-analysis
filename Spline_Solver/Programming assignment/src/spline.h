#include <iostream>
#include <cmath>
#include <algorithm>
#include <vector>
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Sparse>



class Function
{
public:
    virtual double operator()(double _x) = 0; 
    virtual double diff(double _x) 
    {
        return 0;  
    }
        virtual double diff_2(double _x) 
    {
        return 0;  
    }
};


class Interpolation
{
public:
    virtual void solve() = 0;  
    virtual double operator()(double _x) 
    {
        return 0; 
    }  
};


class B_spline : public Function
{
private:
    std::vector<double> t;
    int dim;
public:
    B_spline(int _dim, std::vector<double> _t): dim(_dim), t(_t)  // dim: Order of the B-spline; t: [t_{i-1} to t_{i+dim}]
    {
        if (dim != t.size() - 2)
	    {
		    std::cerr<< "error: the order of B-spline does not match the number of intepolation points on function's support." <<std::endl;
		    exit(-1);
	    }
    } 
    double operator()(double _x)
    {
        std::vector<double> t1(t.begin(), t.end() - 1);
        std::vector<double> t2(t.begin() + 1, t.end());
        B_spline B1(dim - 1, t1);
        B_spline B2(dim - 1, t2);
        if (dim == 1)  // For 1st order B-spline basis functions, directly return the expression value, otherwise, use the recursive formula of the B-spline for recursion
        {
            if (_x > t[0] && _x < t[1])
                return (_x - t[0])/(t[1] - t[0]);
            else if (_x >= t[1] && _x < t[2])
                return (t[2] - _x)/(t[2] - t[1]);
            else
                return 0;
        }
        else
        {
            return (_x - t[0])*B1(_x)/(t[dim] - t[0]) + (t[dim+1] - _x)*B2(_x)/(t[dim+1] - t[1]);
        }
    }
};


class Polynomial : public Function
{
private:
    std::vector<double> Coef;
    double x_0;
    int n;
public:
    /*
    Input Coef = [a_0,...a_n], base = x_0
    Construct a polynomial function f(x) = a_0 + ... + a_n*(x - x_0)^n
    */
    Polynomial(std::vector<double> _Coef, double _base): Coef(_Coef), x_0(_base) {}
    double operator()(double _x)
    {
        double y = Coef[0];
        n = Coef.size();
        for (int i = 1; i < n; i++)
        {
            double delta_y = Coef[i];
            for (int j = 0; j < i; j++)
            {
                delta_y = delta_y*(_x - x_0);
            }
            y = y + delta_y;
        }
        return y;
    }
};

/*
    Discrete form of the function (used for unified input form during curve fitting)
*/
class Discrete_function : public Function
{
private:
    std::vector<double> Point, Value, dF, sec_dF;
    int n;
public:
    /*
    Point: Discrete point set of the function;
    Value, dF, sec_dF: Function values, first derivative and second derivative at the point set
    */
    Discrete_function(std::vector<double> _Point, std::vector<double> _Value, std::vector<double> _dF, std::vector<double> _sec_dF): Point(_Point), Value(_Value), dF(_dF), sec_dF(_sec_dF) 
    {
        if ((Point.size() != Value.size()) || (Value.size() != dF.size()) || (dF.size() != sec_dF.size()))
        {
            std::cerr<< "error: mismatch of input vectors' size." <<std::endl;
		    exit(-1);
        }

        n = Point.size();
    }

    double operator()(double _x)
    {
        for (int i = 0; i < n; i++)
        {
            if (Point[i] == _x)
                return Value[i]; 
        }

        std::cerr<< "error: invalid input." <<std::endl;
		exit(-1);

        return 0;
    }

    double diff(double _x)
    {
        for (int i = 0; i < n; i++)
        {
            if (Point[i] == _x)
                return dF[i]; 
        }

        std::cerr<< "error: invalid input." <<std::endl;
		exit(-1);

        return 0;
    }

    double diff_2(double _x)
    {
        for (int i = 0; i < n; i++)
        {
            if (Point[i] == _x)
                return sec_dF[i]; 
        }

        std::cerr<< "error: invalid input." <<std::endl;
		exit(-1);

        return 0;
    }
};

/*
    Carry out interpolation using B-spline basis functions
*/
class Bspline_interpolation : public Interpolation
{
private:
    std::vector<double> Point, add_point1, add_point2;
    std::vector<B_spline> B;
    Function &f;
    int n, method, order; // n: number of sample points 
    Eigen::VectorXd coef_x; // coef_x: coefficients obtained after solving the equation
public:
    /*
    f: function to be interpolated; Point: interpolation points t_1,...t_n;
    method: interpolation boundary conditions (1 for complete cubic spline, 2 for cubic spline with specified derivatives, 3 for natural cubic spline),
            no difference for linear spline when method input is 1-3;
    order: order of interpolation spline basis function (1 for linear spline, 3 for cubic spline)
    */
    Bspline_interpolation(Function &_f, std::vector<double> _Point, int _method, int _order): f(_f), Point(_Point), method(_method), order(_order)
    {
        if ((method != 1) && (method != 2) && (method != 3))
	    {
		    std::cerr<< "error: the method is not available." <<std::endl;
		    exit(-1);
	    }
        if ((order != 1) && (order != 3))
	    {
		    std::cerr<< "error: the order is not available." <<std::endl;
		    exit(-1);
	    }

        n = Point.size();

        if (order == 3)
        {
            add_point1 = {Point[0]-3*(Point[1]-Point[0]), Point[0]-2*(Point[1]-Point[0]), Point[0]-1*(Point[1]-Point[0])};
            add_point2 = {Point[n-1]+1*(Point[1]-Point[0]), Point[n-1]+2*(Point[1]-Point[0]), Point[n-1]+3*(Point[1]-Point[0])};
        }
        if (order == 1)
        {
            add_point1 = {Point[0]-1*(Point[1]-Point[0])};
            add_point2 = {Point[n-1]+1*(Point[1]-Point[0])};
        }


        Point.insert(Point.begin(), add_point1.begin(), add_point1.end());
        Point.insert(Point.end(), add_point2.begin(), add_point2.end());
        for (int i = 2-order; i < n + 1; i++)
            {
                std::vector<double> v(Point.begin()+i-2+order, Point.begin()+i+2*order);
                B_spline b(order, v);
                B.push_back(b); 
            }
    }

    void solve()
    {
        if (order == 3)
        {
            Eigen::SparseMatrix<double> A(n+2, n+2);
            std::vector<Eigen::Triplet<double> > tripletlist;
            Eigen::MatrixXd y(n+2,1); 
            y = Eigen::MatrixXd::Zero(n+2, 1); // Solving the equation Ax = y

            for (int i = 1; i < n + 1; i++)
		    {
                tripletlist.push_back(Eigen::Triplet<double>(i, i-1, B[i-1](Point[i+2])));
                tripletlist.push_back(Eigen::Triplet<double>(i, i, B[i](Point[i+2])));
                tripletlist.push_back(Eigen::Triplet<double>(i, i+1, B[i+1](Point[i+2])));
                y(i, 0) = f(Point[i+2]);
		    }

            if (method ==  1)
            {
                std::vector<double> v0(Point.begin()+1, Point.begin()+5);
                std::vector<double> v1(Point.begin()+2, Point.begin()+6);
                std::vector<double> vn_1(Point.end()-6, Point.end()-2);
                std::vector<double> vn(Point.end()-5, Point.end()-1);
                B_spline b0(2, v0), b1(2, v1), bn_1(2, vn_1), bn(2, vn);

                tripletlist.push_back(Eigen::Triplet<double>(0, 0, (-3.0*b0(Point[3]))/(Point[4]-Point[1])));
                tripletlist.push_back(Eigen::Triplet<double>(0, 1, (3.0*b0(Point[3]))/(Point[4]-Point[1]) - (3.0*b1(Point[3]))/(Point[5]-Point[2])));
                tripletlist.push_back(Eigen::Triplet<double>(0, 2, (3.0*b1(Point[3]))/(Point[5]-Point[2])));
                tripletlist.push_back(Eigen::Triplet<double>(n+1, n-1, (-3.0*bn_1(Point[n+2]))/(Point[n+3]-Point[n])));
                tripletlist.push_back(Eigen::Triplet<double>(n+1, n, (3.0*bn_1(Point[n+2]))/(Point[n+3]-Point[n]) - (3.0*bn(Point[n+2]))/(Point[n+4]-Point[n+1])));
                tripletlist.push_back(Eigen::Triplet<double>(n+1, n+1, (3.0*bn(Point[n+2]))/(Point[n+4]-Point[n+1])));

                y(0, 0) = f.diff(Point[3]);
                y(n+1, 0) = f.diff(Point[n+2]);
            }
            else
            {
                tripletlist.push_back(Eigen::Triplet<double>(0, 0, 6.0/((Point[4]-Point[1])*(Point[4]-Point[2]))));
                tripletlist.push_back(Eigen::Triplet<double>(0, 1, -6.0/((Point[4]-Point[1])*(Point[4]-Point[2])) - 6.0/((Point[5]-Point[2])*(Point[4]-Point[2]))));
                tripletlist.push_back(Eigen::Triplet<double>(0, 2, 6.0/((Point[5]-Point[2])*(Point[4]-Point[2]))));
                tripletlist.push_back(Eigen::Triplet<double>(n+1, n-1, 6.0/((Point[n+3]-Point[n])*(Point[n+3]-Point[n+1]))));
                tripletlist.push_back(Eigen::Triplet<double>(n+1, n, -6.0/((Point[n+3]-Point[n])*(Point[n+3]-Point[n+1])) - 6.0/((Point[n+4]-Point[n+1])*(Point[n+3]-Point[n+1]))));
                tripletlist.push_back(Eigen::Triplet<double>(n+1, n+1, 6.0/((Point[n+4]-Point[n+1])*(Point[n+3]-Point[n+1]))));
                if (method ==  2)
                {
                    y(0, 0) = f.diff_2(Point[3]);
                    y(n+1, 0) = f.diff_2(Point[n+2]);
                }
                if (method ==  3)
                {
                    y(0, 0) = 0.0;
                    y(n+1, 0) = 0.0;
                }
            }

            A.setFromTriplets(tripletlist.begin(), tripletlist.end());
            A.makeCompressed();
            Eigen::SparseLU<Eigen::SparseMatrix<double> > Solver_sparse;
            Solver_sparse.compute(A);
            coef_x = Solver_sparse.solve(y);
        }
        if (order == 1)
        {
            Eigen::SparseMatrix<double> A(n, n);
            std::vector<Eigen::Triplet<double> > tripletlist;
            Eigen::MatrixXd y(n, 1); 
            y = Eigen::MatrixXd::Zero(n, 1); // Solving the equation Ax = y

            for (int i = 0; i < n; i++)
		    {
                tripletlist.push_back(Eigen::Triplet<double>(i, i, 1.0));
                y(i, 0) = f(Point[i+1]);
		    }

            A.setFromTriplets(tripletlist.begin(), tripletlist.end());
            A.makeCompressed();
            Eigen::SparseLU<Eigen::SparseMatrix<double> > Solver_sparse;
            Solver_sparse.compute(A);
            coef_x = Solver_sparse.solve(y);
        }
    }

    /*
    Return the function value of the interpolation function at _x
    */
    double operator()(double _x)
    {
        if (order == 3)
        {
            if ((_x < Point[3]) || (_x > Point[n+2]))
	        {
		        std::cerr<< "error: the input is out of range." <<std::endl;
		        exit(-1);
	        }
            if (_x == Point[3])
                return f(Point[3]);

            int pos = 1;
            while (_x > Point[pos+3])
            {
                pos = pos + 1;
            }
            
            return coef_x(pos-1)*B[pos-1](_x) + coef_x(pos)*B[pos](_x) + coef_x(pos+1)*B[pos+1](_x) + coef_x(pos+2)*B[pos+2](_x);    
        }
        if (order == 1)
        {
            if ((_x < Point[1]) || (_x > Point[n]))
	        {
		        std::cerr<< "error: the input is out of range." <<std::endl;
		        exit(-1);
	        }
            if (_x == Point[1])
                return f(Point[1]);

            int pos = 1;
            while (_x > Point[pos+1])
            {
                pos = pos + 1;
            }

            return coef_x(pos-1)*B[pos-1](_x) + coef_x(pos)*B[pos](_x);      
        }
        return 0;    
    }

    /*
    If necessary, the Get_coef() function returns the corresponding coefficients of the interpolation polynomial
    */
    Eigen::VectorXd Get_coef()
    {
        return coef_x;
    }
};

/*
    Carry out interpolation using piecewise polynomials
*/
class ppForm_interpolation : public Interpolation
{
private:
    std::vector<double> Point;
    Function &f;
    std::vector<Polynomial> P; // Store n-1 piecewise polynomials obtained from interpolation 
    std::vector<std::vector<double> > C; // Store coefficients of the piecewise polynomial
    int n, method, order; // n: number of sample points 
    Eigen::VectorXd coef_m; // coef_m: coefficients obtained after solving the equation
public:
    /*
    f: function to be interpolated; Point: interpolation points t_1,...t_n;
    method: interpolation boundary conditions (1 for complete cubic spline, 2 for cubic spline with specified derivatives, 3 for natural cubic spline);
            no difference for linear spline when method input is 1-3;
    order: order of interpolation spline basis function (1 for linear spline, 3 for cubic spline)
    */
    ppForm_interpolation(Function &_f, std::vector<double> _Point, int _method, int _order): f(_f), Point(_Point), method(_method), order(_order)
    {
        if ((method != 1) && (method != 2) && (method != 3))
	    {
		    std::cerr<< "error: the method is not available." <<std::endl;
		    exit(-1);
	    }
        if ((order != 1) && (order != 3))
	    {
		    std::cerr<< "error: the order is not available." <<std::endl;
		    exit(-1);
	    }

        n = Point.size();
    }

    void solve()
    {
        if (order == 3)
        {
            Eigen::SparseMatrix<double> A(n, n);
            std::vector<Eigen::Triplet<double> > tripletlist;
            Eigen::MatrixXd y(n, 1); 
            y = Eigen::MatrixXd::Zero(n, 1); 

            for (int i = 1; i < n - 1; i++)
		    {
                tripletlist.push_back(Eigen::Triplet<double>(i, i-1, (Point[i+1]-Point[i])/(Point[i+1]-Point[i-1])));
                tripletlist.push_back(Eigen::Triplet<double>(i, i, 2.0));
                tripletlist.push_back(Eigen::Triplet<double>(i, i+1, (Point[i]-Point[i-1])/(Point[i+1]-Point[i-1])));
                y(i, 0) = 3.0*(Point[i]-Point[i-1])/(Point[i+1]-Point[i-1])*(f(Point[i+1])-f(Point[i]))/(Point[i+1]-Point[i]) + 3*(Point[i+1]-Point[i])/(Point[i+1]-Point[i-1])*(f(Point[i])-f(Point[i-1]))/(Point[i]-Point[i-1]);
		    }

            if (method ==  1)
            {
                tripletlist.push_back(Eigen::Triplet<double>(0, 0, 1.0));
                tripletlist.push_back(Eigen::Triplet<double>(n-1, n-1, 1.0));
                y(0, 0) = f.diff(Point[0]);
                y(n-1, 0) = f.diff(Point[n-1]);
            }
            else
            {
                tripletlist.push_back(Eigen::Triplet<double>(0, 0, 4.0));
                tripletlist.push_back(Eigen::Triplet<double>(0, 1, 2.0));
                tripletlist.push_back(Eigen::Triplet<double>(n-1, n-2, 2.0));
                tripletlist.push_back(Eigen::Triplet<double>(n-1, n-1, 4.0));
                if (method ==  2)
                {
                    y(0, 0) = 6.0*(f(Point[1])-f(Point[0]))/(Point[1]-Point[0]) - f.diff_2(Point[0])*(Point[1]-Point[0]);
                    y(n-1, 0) = 6.0*(f(Point[n-1])-f(Point[n-2]))/(Point[n-1]-Point[n-2]) + f.diff_2(Point[n-1])*(Point[n-1]-Point[n-2]);
                }
                if (method ==  3)
                {
                    y(0, 0) = 6.0*(f(Point[1])-f(Point[0]))/(Point[1]-Point[0]);
                    y(n-1, 0) = 6.0*(f(Point[n-1])-f(Point[n-2]))/(Point[n-1]-Point[n-2]);
                }
            }

            A.setFromTriplets(tripletlist.begin(), tripletlist.end());
            A.makeCompressed();
            Eigen::SparseLU<Eigen::SparseMatrix<double> > Solver_sparse;
            Solver_sparse.compute(A);
            coef_m = Solver_sparse.solve(y);

            for (int i = 1; i < n; i++)
		    {
                std::vector<double> pp;
                double K; // Ki
                K = (f(Point[i])-f(Point[i-1]))/(Point[i]-Point[i-1]);
                pp.push_back(f(Point[i-1]));
                pp.push_back(coef_m(i-1));
                pp.push_back((3.0*K-2.0*coef_m(i-1)-coef_m(i))/(Point[i]-Point[i-1]));
                pp.push_back((coef_m(i-1)+coef_m(i)-2.0*K)/pow(Point[i]-Point[i-1],2));
                C.push_back(pp);

                Polynomial p(pp, Point[i-1]);
                P.push_back(p); 
		    }
        }
        if (order == 1)
        {
            for (int i = 1; i < n; i++)
		    {
                std::vector<double> pp;
                double K; // Ki
                K = (f(Point[i])-f(Point[i-1]))/(Point[i]-Point[i-1]);
                pp.push_back(f(Point[i])-K*Point[i]);
                pp.push_back(K);
                C.push_back(pp);
                
                Polynomial p(pp, 0);
                P.push_back(p); 
		    }
        }
    }

    /*
    Return the function value of the interpolation function at _x
    */
    double operator()(double _x)
    {
        if ((_x < Point[0]) || (_x > Point[n-1]))
	    {
		    std::cerr<< "error: the input is out of range." <<std::endl;
		    exit(-1);
	    }
        if (_x == Point[0])
            return f(Point[0]);

        int pos = 1;
        while (_x > Point[pos])
        {
            pos = pos + 1;
        }

        if (_x == Point[pos])
            return f(Point[pos]);
        else
            return P[pos-1](_x);
    }

    /*
    If necessary, the Get_coef() function returns the corresponding coefficients of the piecewise interpolation polynomial
    */
    std::vector<std::vector<double> > Get_coef()
    {
        return C;
    }
};

/*
    Carry out spline interpolation for curves
*/
class curve_spline : public Interpolation
{
private:
    const double interval = 0.01;
    std::vector<std::vector<double> > Point, Fit_Point; 
    std::vector<double> cc_length; // Store cumulative chordal length
    int n, d, method, order; 
public:
    /*
    Point: sample points on the curve (if there are n points, each point is d-dimensional, then input a 2D vector of n*d);
    method: interpolation method (1 for B-spline interpolation, 2 for ppForm)
    order: order of interpolation spline basis function (1 for linear spline, 3 for cubic spline)
    */
    curve_spline(std::vector<std::vector<double> > _Point, int _method, int _order): Point(_Point), method(_method), order(_order)
    {
        if ((method != 1) && (method != 2))
	    {
		    std::cerr<< "error: the method is not available." <<std::endl;
		    exit(-1);
	    }
        if ((order != 1) && (order != 3))
	    {
		    std::cerr<< "error: the order is not available." <<std::endl;
		    exit(-1);
	    }

        n = Point.size();
        d = Point[0].size();

        for (int i = 0; i < n; i++)
		{
            if (Point[i].size() != d)
            {
                std::cerr<< "error: mismatch of sample points' dimension." <<std::endl;
		        exit(-1);
            }
		}
    }


    double norm_2(std::vector<double> x, std::vector<double> y)
    {
        double dist = 0.0;
        for (int i = 0; i < d; i++)
		{
            dist = dist + (x[i]-y[i])*(x[i]-y[i]);
		}

        return sqrt(dist);
    }

    void solve()
    {
        // calculate the cumulative chordal length
        cc_length.push_back(0.0);
        for (int i = 1; i < n; i++)
		{
            cc_length.push_back(cc_length[i-1] + norm_2(Point[i], Point[i-1]));
		}

        std::vector<double> EMPTY(n);
        for (int i = 0; i < d; i++)
		{
            std::vector<double> temp_p, temp_fitp;
            for (int j = 0; j < n; j++)
		    {
                temp_p.push_back(Point[j][i]);
		    }
            Discrete_function temp_f(cc_length, temp_p, EMPTY, EMPTY);
            if (method == 1)
            {
                Bspline_interpolation BSI(temp_f,cc_length, 3, order);
                BSI.solve();
                double length = 0;
                while (length <= cc_length[n-1])
                {
                    temp_fitp.push_back(BSI(length));
                    length = length + interval;
                }
            }
            if (method == 2)
            {
                ppForm_interpolation ppF(temp_f, cc_length, 3, order);
                ppF.solve();
                double length = 0;
                while (length <= cc_length[n-1])
                {
                    temp_fitp.push_back(ppF(length));
                    length = length + interval;
                }
            }
            Fit_Point.push_back(temp_fitp);
		}
    }

    /*
    Get_Point() function returns d vectors, where the i-th vector represents the i-th component of the points sampled at equal intervals after the curve spline interpolation
    */
    std::vector<std::vector<double> > Get_Point()
    {
        return Fit_Point;
    }
};
