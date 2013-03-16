//
//  Simplex.cpp
//  LinearProgramming
//
//  Created by Guoxin Jin on 2/3/13.
//  Copyright (c) 2013 Guoxin Jin. All rights reserved.
//

#include "Simplex.h"
Simplex::Simplex()
{
    c_num=0;
    v_num=0;
}

Simplex::~Simplex()
{
}
Simplex::Simplex(MatrixXd A,VectorXd b,VectorXd c, VectorXd base_index)
{
  _A = A;
  _b = b;
  _c = c;
  J = 0;
  c_num = (int) A.rows();
  v_num = (int) A.cols();
  B_index = base_index;
  N_index = VectorXd(v_num-c_num);
  int count=0;
  int idx=0;
  while(idx<v_num)
  {
    bool find = false;
    for (int i=0;i<c_num;i++)
    {
      if(B_index(i)==idx)
      {
        find = true;
        break;
      }
    }
    if (!find)
    {
      N_index(count)=idx;
      count++;
    }
    idx++;
  }
  //cout<<N_index<<endl;
  _B = MatrixXd(c_num,c_num);
  _N = MatrixXd(c_num,v_num-c_num);
  _x = VectorXd(v_num);
  getB();
  getN();
  _Binv  = _B.inverse();
  
}
Simplex::Simplex(MatrixXd A,VectorXd b,VectorXd c)
{
  _A = A;
  _b = b;
  _c = c;
  c_num = (int) A.rows();
  v_num = (int) A.cols();
  B_index = VectorXd(c_num);
  N_index = VectorXd(v_num-c_num);
  for (int i=0; i<c_num;i++)
    B_index(i)=i;
  int count=0;
  int idx=0;
  while(idx<v_num)
  {
    bool find = false;
    for (int i=0;i<c_num;i++)
    {
      if(B_index(i)==idx)
      {
        find = true;
        break;
      }
    }
    if (!find)
    {
      N_index(count)=idx;
      count++;
    }
    idx++;
  }
  _x = VectorXd(v_num);
  _B = MatrixXd(c_num,c_num);
  _N = MatrixXd(c_num,v_num-c_num);
  _x = VectorXd(v_num);
  getB();
  getN();
  _Binv  = _B.inverse();
  
}
void Simplex::pivoting(MatrixXd& A, int row, int col)
{
  double pivot = A(row,col);
  A.row(row)/=pivot;
  double ratio = pivot;
  for (int i=0;i<A.rows();i++)
  {
    if(i==row)
    {
        continue;
    }
    ratio = A(i,col);
    A.row(i)-=A.row(row)*ratio;
  }
  return;
    
}

void  Simplex::reduceCost(Eigen::VectorXd& cB, Eigen::MatrixXd& Binv,
                Eigen::MatrixXd& N, Eigen::VectorXd& cN,
                Eigen::VectorXd& rst)
{
  //cout<<"Binv:\n"<<Binv<<endl;
  //cout<<"N:\n"<<N<<endl;
  VectorXd temp = (cB.transpose()*Binv)*N ;
  rst = temp - cN;
  //cout<<"z-c:\n"<<rst<<endl;
  //inverse()*N-cN;
  //rst = rst.transpose();
  return;
}

void Simplex::getSlice(const Eigen::MatrixXd& src, const Eigen::VectorXd& col_index, Eigen::MatrixXd& rst)
{
  int count = 0;
  while (count< rst.cols()) {
    rst.col(count) = src.col(col_index(count));
    count++;
  }
  return;
}

void Simplex::getSlice(const Eigen::VectorXd& src, const Eigen::VectorXd& col_index, Eigen::VectorXd& rst)
{
  int count = 0;
  while (count< rst.size()) {
    rst(count) = src(col_index(count));
    count++;
  }
  return;
}
void Simplex::setSlice(const Eigen::MatrixXd& src, const Eigen::VectorXd& col_index, MatrixXd& rst)
{
  int count = 0;
  while (count< src.cols()) {
    rst.col(col_index(count)) = src.col(count);
    count++;
  }
  return;
}

Eigen::MatrixXd& Simplex::getB()
{
  getSlice(_A, B_index, _B);
  return _B;
}

Eigen::MatrixXd& Simplex::getN()
{
  getSlice(_A, N_index, _N);
  return _N;
}

void Simplex::initTable()
{
  _Tab = Eigen::MatrixXd(c_num+1, v_num+1);
  Eigen::VectorXd cost(v_num);
  VectorXd cB(c_num);
  VectorXd cN(v_num-c_num);
  getSlice(_c, B_index, cB);
  getSlice(_c, N_index, cN);
  //reduceCost(cB, getB(), getN(), cN, cost);
  //VectorXd cost(v_num);
  //cout<<"cB:\n"<<cB<<endl;
  //cout<<"cN:\n"<<cN<<endl;
  _Binv = _B.inverse();
  reduceCost(cB, _Binv, _A, _c,cost );
  //cout<<cost<<endl;
  
  _Tab.block(1,0,c_num,v_num) = _Binv*_A;
  
  _Tab.block(0,0,1,v_num) = cost.transpose();
  J = cB.transpose()*_Binv*_b;
  _Tab(0,v_num)=J;
  _Tab.block(1,v_num,c_num,1) = _Binv*_b;
  for (int i=0; i< c_num;i++)
  {
    _x(B_index(i))=_Tab(i+1,v_num);
  }
  /*
  MatrixXd BxN = _B.inverse()*_N;
  MatrixXd temp = _Tab.block(1,0,c_num,v_num);
  setSlice(BxN, N_index, temp);
  setSlice(MatrixXd::Identity(c_num,c_num),B_index,temp);
  _Tab.block(1,0,c_num,v_num)=temp;
  for (int i=0;i<c_num;i++)
  {
    _Tab(0,B_index(i)) = 0;
    _Tab(i)
  }
  */
  
  cout<<"_Tab:\n"<<_Tab<<endl;
  return;
}

double Simplex::updateTable()
{
  Eigen::VectorXd cost = _Tab.row(0);
  int i,j;
  double min = cost.minCoeff(&j);//col, row
  Eigen::VectorXd theta;
  //j is the next basis select in
  if (min<0)
  {
    theta = _Tab.block(1,v_num,c_num,1).cwiseQuotient(_Tab.block(1,j,c_num,1));
    cout<<"theta:\n"<<theta<<endl;
    theta.minCoeff(&i);
    B_index(i)= j;
    cout<<"\n B_index: "<<B_index.transpose();
    i++;
    pivoting(_Tab, i, j);
    J = _Tab(0,v_num);
  }
  cout<<"_Tab:\n"<<_Tab<<endl;
  return J;
}

Eigen::VectorXd& Simplex::SolveLP(int method)
{
  double ep=0.001;
  double delta = 0.0001;
  double alpha = 0.98;
  double old_J=-INFINITY;
  initTable();
  _x<<4,10,14,26;
  cout<<"\nInit x: "<<_x.transpose()<<endl;
  
  J = _c.transpose()*_x;
  do {
    old_J = J;
    if (method ==0)
    {
      J=updateTable();
    }
    else
    {
      J=affineScaling(_x, alpha);
    }
    //cout<<"J: "<<J<<", Old_J: "<<old_J<<endl;
  } while (abs(J-old_J)/(abs(old_J)+delta)>ep);
  if (method==0)
  {  _x = _Tab.block(1,v_num,c_num,1);
    cout<<"index of basis:\n"<<B_index.transpose()<<endl;
  }
  cout<<"solution:\n"<<_x;
  return _x;
}

VectorXd Simplex::getX()
{
  return _x;
}
double   Simplex::getZ()
{
  return J;
}
VectorXd Simplex::getIdx()
{
  return B_index;
}
double Simplex::affineScaling(VectorXd &x, double alpha)
{
  // the return value indicat the status of solution
  // <0 unbounded
  // =0 not converged yet
  // >0 converged
  MatrixXd D = x.asDiagonal();
  //MatrixXd iD = D.inverse();
  Eigen::MatrixXd Ak = _A*D;
  Eigen::MatrixXd L = Ak*(Ak.transpose());
  Eigen::VectorXd wk = L.inverse()*Ak*D*_c;
  cout<<"\n wk=: "<<wk.transpose()<<endl;
  Eigen::VectorXd rk = _c- _A.transpose()*wk;
  cout<<"\n rk=: "<<rk.transpose()<<endl;
  Eigen::VectorXd dk = D*rk;
  cout<<"\n dk=: "<<dk.transpose()<<endl;
  double ak=INFINITY;
  bool unbound=true;
  for (int i=0; i< dk.size();i++)
  {
    if (dk(i)>0)
      unbound=unbound&&true;
    else
    {
      unbound=unbound&&false;
      if (dk(i)!=0)
      {
        double temp =alpha/(-1*dk(i));
        if (temp<ak)
          ak = temp;
      }
    }
  }
  if (unbound)
    return -INFINITY;
  //Eigen::VectorXd old_x = _x;
  _x = _x + ak*D*dk;
  cout<<"\n xk= "<<_x.transpose()<<endl;
  J = _c.transpose()*_x;
  cout<<"\n Jk= "<<J<<endl;
  //auto diff = (_c.transpose()*_x - _b.transpose()*wk);
  //diff.abs();
  
  return J;
}
