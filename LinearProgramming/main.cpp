//
//  main.cpp
//  LinearProgramming
//
//  Created by Guoxin Jin on 2/3/13.
//  Copyright (c) 2013 Guoxin Jin. All rights reserved.
//

#include <iostream>
#include <Eigen/Dense>
#include "Simplex.h"
using std::cout;
using std::endl;
using Eigen::MatrixXd;
using Eigen::VectorXd;
int main(int argc, const char * argv[])
{

    // insert code here...
    int m = 2;
    int n = 4;
    MatrixXd A(m,n);
    VectorXd c(n);
    VectorXd b(m);
    VectorXd index(m);
    /*
    A << 4, 3, 1, 0, 0,
    2, 3, 0, 1, 0,
    4, 2, 0, 0, 1;
    c << 9,12,0,0,0;
    b << 180, 150, 160;
    
    index<<2,3,4;
    */
    A<<-1, 2,1,0,
    1,0,0,1;
    c<<-1,4,0,0;
    b<<30,30;
    Simplex sim(A,b,c);
    VectorXd sol = sim.SolveLP(1);
    index = sim.getIdx();
    double J = sim.getZ();
    cout<<"the solution is:\n"<<index.transpose()<<endl<<sol.transpose()<<endl;
    cout<<"the optimal value is:  "<<J<<endl;
    return 0;
}

