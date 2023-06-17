//
//  Fluid.hpp
//  Mable
//
//  Created by 須之内俊樹 on 2023/02/26.
//

#ifndef Fluid_hpp
#define Fluid_hpp

#include <iostream>
#include <vector>
#include <algorithm>
#include <string>
#include <sstream>
#include <iomanip>
#include <set>
#include <Eigen/Core>
#include <Eigen/IterativeLinearSolvers>
#include "Array3d.hpp"
#define Nx 16
#define Ny 16
#define Nz 16//グリッドの数

#define Tamb 25
#define g0 9.8
using ScalarType = double;
using IndexType = int64_t;
using Triplet = Eigen::Triplet<ScalarType,IndexType>;
using SparseMatrix = Eigen::SparseMatrix<ScalarType>;
typedef Eigen::SparseMatrix<double, Eigen::RowMajor, int64_t> SpMat;
struct Fluid{
    double dx;//セルの大きさ
    double dt;//時間の刻み幅
    //double rho;
    myArray3d u = myArray3d(Nx+1,Ny,Nz,0);//水平
    myArray3d old_u = myArray3d(Nx+1,Ny,Nz,0);//水平
    myArray3d v = myArray3d(Nx,Ny+1,Nz,0);//鉛直
    myArray3d old_v = myArray3d(Nx,Ny+1,Nz,0);//水平
    myArray3d w = myArray3d(Nx,Ny,Nz+1,0);
    myArray3d old_w = myArray3d(Nx,Ny,Nz+1,0);
    myArray3d p = myArray3d(Nx,Ny,Nz,0);//圧力
    myArray3d rho = myArray3d(Nx,Ny,Nz,0);
    myArray3d temp = myArray3d(Nx,Ny,Nz,Tamb);
    
    //double L;
    Fluid();
    //Fluid(double x,double t,double density);
    Fluid(double x,double t,myArray3d &density,myArray3d &temprature);
    //std::vector<int>DirichletBoundaryCondition(int i,int j,int k,myMap &map);
    void project();
//    void project(myMap &map);
    void faceAdvect();
    void centerAdvect(myArray3d &val);
    double TriLinearInterporation(double x,double y,double z,myArray3d &val);
    
};
#endif /* Fluid_hpp */
