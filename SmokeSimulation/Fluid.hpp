//
//  Fluid.hpp
//  Mable
//
//  Created by 須之内俊樹 on 2023/02/26.
// g++ main.cpp Fluid.cpp rendering.cpp Array3d.cpp output.cpp -I /Users/sunouchi/c++/library/eigen-3.4.0 -I /Users/sunouchi/c++/library/openmpi-4.1.5 -std=c++2a -fopenmp -O3 -I /Users/sunouchi/c++/library/MyLib

#ifndef Fluid_hpp
#define Fluid_hpp

#include <iostream>
#include <filesystem>
#include <vector>
#include <algorithm>
#include <string>
#include <sstream>
#include <iomanip>
#include <set>
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/IterativeLinearSolvers>
#include "Array3d.hpp"
#include "output.hpp"
#include "rendering.hpp"
#define Nx 32
#define Ny 32
#define Nz 64 //グリッドの数
#define Lx 1.0
#define Ly 1.0
#define Lz 2.0
#define range 1
#define Tamb 0
#define p0 0
#define Ramb 1
#define g0 9.8
#define beta 75
#define epcilon 1e-5
#define timestep 200
using ScalarType = double;
using IndexType = int64_t;
using Triplet = Eigen::Triplet<ScalarType,IndexType>;
using SparseMatrix = Eigen::SparseMatrix<ScalarType>;
typedef Eigen::SparseMatrix<double, Eigen::RowMajor, int64_t> SpMat;
struct Fluid{
    double dx;//セルの大きさ
    double dt;//時間の刻み幅
    TimeDisplayer TD;
    std::vector<double> times;
    //double rho;
    myArray3<double> u = myArray3<double>(Nx+1,Ny,Nz,0);//水平
    myArray3<double> old_u = myArray3<double>(Nx+1,Ny,Nz,0);//水平
    myArray3<double> v = myArray3<double>(Nx,Ny+1,Nz,0);//鉛直
    myArray3<double> old_v = myArray3<double>(Nx,Ny+1,Nz,0);//水平
    myArray3<double> w = myArray3<double>(Nx,Ny,Nz+1,0);
    myArray3<double> old_w = myArray3<double>(Nx,Ny,Nz+1,0);
    myArray3<double> p = myArray3<double>(Nx,Ny,Nz,p0);//圧力
    myArray3<double> rho_tgt = myArray3<double>(Nx,Ny,Nz,0);
    myArray3<double> rho_amb = myArray3<double>(Nx,Ny,Nz,Ramb);
    myArray3<double> temp = myArray3<double>(Nx,Ny,Nz,Tamb);
    myArray3<Eigen::Vector3d> centerRot = myArray3<Eigen::Vector3d>(Nx,Ny,Nz,Eigen::Vector3d::Zero());
    myArray3<Eigen::Vector3d> f = myArray3<Eigen::Vector3d>(Nx,Ny,Nz,Eigen::Vector3d::Zero());
    
    myArray3<double> vox_trans = myArray3<double>(Nx,Ny,Nz,0);
    //double L;
    Fluid();
    //Fluid(double x,double t,double density);
    Fluid(double x,double t,myArray3<double> &density,myArray3<double> &templature);
    //std::vector<int>DirichletBoundaryCondition(int i,int j,int k,myMap &map);
    void oneloop();
    void execute();
    void setDensity(int set_range);
    void setTemplature(int set_range);
    void setPressure(int set_range);
    void setV(int set_range);
    void project();
//    void project(myMap &map);
    void faceAdvect();
    void centerAdvect(myArray3<double> &val,double boundary_value);
    Eigen::Vector3d getBuoyanacy(int i,int j, int k);
    void setCenterRot();
    Eigen::Vector3d getConfinement(int i,int j,int k);
    double TriLinearInterporation(double x,double y,double z,myArray3<double> &val);
    void addForce();
    
};
#endif /* Fluid_hpp */
