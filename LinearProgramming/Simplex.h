//
//  Simplex.h
//  LinearProgramming
//
//  Created by Guoxin Jin on 2/3/13.
//  Copyright (c) 2013 Guoxin Jin. All rights reserved.
//

#ifndef __LinearProgramming__Simplex__
#define __LinearProgramming__Simplex__

#include <iostream>
#include <Eigen/Dense>
using Eigen::MatrixXd;
using Eigen::VectorXd;
using namespace std;
class Simplex
{
public:
    Simplex();
    Simplex(MatrixXd A,VectorXd b,VectorXd c, VectorXd base_index);
    Simplex(MatrixXd A,VectorXd b,VectorXd c);
    ~Simplex();
    VectorXd& SolveLP(int method=0);//0 Simplex method, 1 Interior point method
    VectorXd getX();
    double   getZ();
    VectorXd getIdx();
    Eigen::MatrixXd& getB();
    Eigen::MatrixXd& getN();
private:
    int c_num;//constraint num;
    int v_num;//variable num;
    Eigen::VectorXd B_index, N_index;
    //private:
    void pivoting(MatrixXd& A, int row, int col);
    void reduceCost(Eigen::VectorXd& cB, Eigen::MatrixXd& B,
                    Eigen::MatrixXd& N, Eigen::VectorXd& cN,
                    Eigen::VectorXd& rst);
    void getSlice(const MatrixXd& src, const Eigen::VectorXd& col_index, Eigen::MatrixXd& rst);
    void getSlice(const Eigen::VectorXd& src, const Eigen::VectorXd& col_index, Eigen::VectorXd& rst);
    void setSlice(const MatrixXd& src, const Eigen::VectorXd& col_index, Eigen::MatrixXd& rst);
    
    double affineScaling(VectorXd& x, double alpha);
    //Eigen::MatrixXd& getB();
    //Eigen::MatrixXd& getN();
    void initTable();
    //VectorXd& SolveLP();
    double updateTable();
    
private:
    Eigen::MatrixXd _A;
    Eigen::VectorXd _c;
    Eigen::VectorXd _b;
    Eigen::VectorXd _x;
    Eigen::VectorXd _z;
    Eigen::MatrixXd _B;
    Eigen::MatrixXd _Binv;
    Eigen::MatrixXd _N;
    Eigen::MatrixXd _Tab;
    double J;//target function
    
};

#endif /* defined(__LinearProgramming__Simplex__) */
