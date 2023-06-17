//
//  Fluid.cpp












//力，密度，温度の初期化(Template Array3<double>)，温度と渦の力計算，計算のループ作成，密度を使ってテキトーな光の計算，平行投影のなんちゃって横断アルゴリズム




















#include "Fluid.hpp"
Fluid::Fluid(){}
//Fluid::Fluid(double x,double t,double density){
//    dx = x;
//    dt = t;
//    rho = density;
//    L = dx*Nx;
//    vfi.reset(f0.y());
//}
Fluid::Fluid(double x,double t,myArray3d &density,myArray3d &temprature){
    dx = x;
    dt = t;
    rho = density;
    temp = temprature;
}
double Fluid::TriLinearInterporation(double x,double y,double z,myArray3d &val){
    x = fmax(0.0, fmin(val.nx-1-1e-6,x/dx));
    y = fmax(0.0, fmin(val.ny-1-1e-6,y/dx));
    z = fmax(0.0, fmin(val.nz-1-1e-6,z/dx));
    int i = x;int j = y;int k = z;
    double s = x-i;double t = y-i;double u = z-i;
    Eigen::VectorXd f = {
        val.value[i][j][k],val.value[i+1][j][k],val.value[i][j+1][k],val.value[i][j][k+1],
        val.value[i+1][j+1][k],val.value[i][j+1][k+1],val.value[i+1][j][k+1],val.value[i+1][j+1][k+1]};
    Eigen::VectorXd c = {
        (1-s)*t*u,s*t*u,(1-s)*(1-t)*u,(1-s)*t*(1-u),
        s*(1-t)*u,(1-s)*(1-t)*(1-u),s*t*(1-u),s*(1-t)*(1-u)};
    return f.dot(c);
}
void Fluid::faceAdvect(){
    //x
    for(int i=1;i<u.nx;++i){
        for(int j=0;j<u.ny;++j){
            for(int k=0;j<u.nz;++k){
                double x = i*dx;double y = (j+0.5)*dx;double z = (k+0.5)*dx;
                double adv_x = x - dt*TriLinearInterporation(x, y-0.5*dx, z-0.5*dx, old_u);
                double adv_y = y - dt*TriLinearInterporation(x-0.5*dx, y, z-0.5*dx, old_v);
                double adv_z = z - dt*TriLinearInterporation(x-0.5*dx, y-0.5*dx, z, old_w);
                u.value[i][j][k] = TriLinearInterporation(adv_x, adv_y, adv_z, old_u);
            }
        }
    }
    for(int i=0;i<v.nx;++i){
        for(int j=1;j<v.ny;++j){
            for(int k=0;j<v.nz;++k){
                double x = (i+0.5)*dx;double y = j*dx;double z = (k+0.5)*dx;
                double adv_x = x - dt*TriLinearInterporation(x, y-0.5*dx, z-0.5*dx, old_u);
                double adv_y = y - dt*TriLinearInterporation(x-0.5*dx, y, z-0.5*dx, old_v);
                double adv_z = z - dt*TriLinearInterporation(x-0.5*dx, y-0.5*dx, z, old_w);
                v.value[i][j][k] = TriLinearInterporation(adv_x, adv_y, adv_z, old_v);
            }
        }
    }
    for(int i=0;i<w.nx;++i){
        for(int j=0;j<w.ny;++j){
            for(int k=1;j<w.nz;++k){
                double x = (i+0.5)*dx;double y = (j+0.5)*dx;double z = k*dx;
                double adv_x = x - dt*TriLinearInterporation(x, y-0.5*dx, z-0.5*dx, old_u);
                double adv_y = y - dt*TriLinearInterporation(x-0.5*dx, y, z-0.5*dx, old_v);
                double adv_z = z - dt*TriLinearInterporation(x-0.5*dx, y-0.5*dx, z, old_w);
                w.value[i][j][k] = TriLinearInterporation(adv_x, adv_y, adv_z, old_w);
            }
        }
    }
}
void Fluid::centerAdvect(myArray3d &val){
    myArray3d old_val = val;
    for(int i=0;i<val.nx;++i){
        for(int j=0;j<val.ny;++j){
            for(int k=0;j<val.nz;++k){
                double x = (i+0.5)*dx;double y = (j+0.5)*dx;double z = (k+0.5)*dx;
                double adv_x = x - dt*TriLinearInterporation(x, y-0.5*dx, z-0.5*dx, u);
                double adv_y = y - dt*TriLinearInterporation(x-0.5*dx, y, z-0.5*dx, v);
                double adv_z = z - dt*TriLinearInterporation(x-0.5*dx, y-0.5*dx, z, w);
                val.value[i][j][k] = TriLinearInterporation(adv_x, adv_y, adv_z, old_val);
            }
        }
    }
}

//std::vector<int>Fluid::DirichletBoundaryCondition(int i,int j,int k,myMap &map){
//    std::vector<int>ret(6,1);
//    std::vector<std::vector<int>>keys = {{i+1,j,k},{i,j+1,k},{i-1,j,k},{i,j-1,k},{i,j,k-1},{i,j,k+1}};
//    if(i == Nx-1)ret[0] = 0;
//    else if(!map.contains(keys[0]))ret[0] = 0;
//    if(j == Ny-1)ret[1] = 0;
//    else if(!map.contains(keys[1]))ret[1] = 0;
//
//    if(i == 0)ret[2] = 0;
//    else if(!map.contains(keys[2]))ret[2] = 0;
//    if(j == 0)ret[3] = 0;
//    else if(!map.contains(keys[3]))ret[3] = 0;
//
//    if(k == 0)ret[4] = 0;
//    else if(!map.contains(keys[4]))ret[4] = 0;
//    if(k == Nz-1)ret[5] = 0;
//    else if(!map.contains(keys[5]))ret[5] = 0;
//    return ret;
//}
//void Fluid::project(myMap &map){
void Fluid::project(){
    SparseMatrix A(Nx*Ny*Nz,Nx*Ny*Nz),B(Nx*Ny*Nz,Nx*Ny*Nz);
    Eigen::VectorXd b = Eigen::VectorXd::Zero(Nx*Ny*Nz);
    Eigen::VectorXd px;
//    std::set<int> DirichletKey;
//    std::vector<std::vector<int>>keys;
    //Tripletの計算
    std::vector<Triplet> triplets;
    for(int i=0;i<Nx;i++){
        for(int j=0;j<Ny;j++){
            for(int k=0;k<Nz;k++){
//                std::vector<int>key = {i,j,k};
//                if(!map.contains(key)){
//                    //前処理でAが変更されてしまうので，境界条件として別で無理矢理設定する．
//                    triplets.emplace_back(i+j*Nx+k*Nx*Ny,i+j*Nx+k*Nx*Ny,1);
//                    keys.push_back(key);
//                    DirichletKey.insert(i+j*Nx+k*Nx*Ny);
//                    continue;
//                }
                double scale = dt/(rho.value[i][j][k]*dx*dx);
                //std::cout << i << "," << j << std::endl;
                double D[6] = {1.0,1.0,-1.0,-1.0,-1.0,1.0};//周囲6方向に向かって働く、圧力の向き
                //double F[4] = {(double)(i<Nx-1),(double)(j<Ny-1),(double)(i>0),(double)(j>0)};//境界条件。壁なら0,流体なら1
                
                std::vector<int> F = {i<Nx-1,j<Ny-1,i>0,j>0,k>0,k<Nz-1};
                double U[6] = {
                    u.value[i+1][j][k],
                    v.value[i][j+1][k],
                    u.value[i][j][k],
                    v.value[i][j][k],
                    w.value[i][j][k],
                    w.value[i][j][k+1]};
                double sumP = 0;
                for(int n=0;n<6;n++){
                    sumP += -F[n]*scale;
                    //sumP += scale;
                    b(i+j*Nx+k*Nx*Ny) += D[n]*F[n]*U[n]/(dx);
                }
                //F = DirichletBoundaryCondition(i,j,k,map);
//                    for(int n=0;n<6;n++){
//                        std::cout << F_pri[n] << "," << F[n] << std::endl;
//                    }
                triplets.emplace_back(i+j*Nx+k*Nx*Ny,i+j*Nx+k*Nx*Ny, sumP);
                if(F[0])triplets.emplace_back(i+j*Nx+k*Nx*Ny,i+1+j*Nx+k*Nx*Ny, F[0]*scale);
                if(F[1])triplets.emplace_back(i+j*Nx+k*Nx*Ny,i+(j+1)*Nx+k*Nx*Ny, F[1]*scale);
                if(F[2])triplets.emplace_back(i+j*Nx+k*Nx*Ny,i-1+j*Nx+k*Nx*Ny, F[2]*scale);
                if(F[3])triplets.emplace_back(i+j*Nx+k*Nx*Ny,i+(j-1)*Nx+k*Nx*Ny, F[3]*scale);
                if(F[4])triplets.emplace_back(i+j*Nx+k*Nx*Ny,i+j*Nx+(k-1)*Nx*Ny, F[4]*scale);
                if(F[5])triplets.emplace_back(i+j*Nx+k*Nx*Ny,i+j*Nx+(k+1)*Nx*Ny, F[5]*scale);
            }
        }
    }
    A.setFromTriplets(triplets.begin(), triplets.end());
    Eigen::ConjugateGradient<SparseMatrix> solver;
//    for(int i=0;i<A.outerSize();++i){
//        for(SparseMatrix::InnerIterator it(A,i);it;++it){
//            if(it.row() == *DirichletKey.begin()){
//                it.valueRef() = 1;
//                DirichletKey.erase(DirichletKey.begin());
//            }
//        }
//    }
    solver.compute(A);
    px = solver.solve(b);
    for(int i=0;i<Nx;i++){
        for(int j=0;j<Ny;j++){
            for(int k=0;k<Nz;k++){
                p.value[i][j][k] = px(i+j*Nx+k*Nx*Ny);
            }
        }
    }
    for(int i=1; i<Nx;i++){
        for(int j=0;j<Ny;j++){
            for(int k=0;k<Nz;k++)u.value[i][j][k] = u.value[i][j][k] - dt/rho.value[i][j][k] * (p.value[i][j][k]-p.value[i-1][j][k])/dx;
        }
    }
    for(int i=0;i<Nx;i++){
        for(int j=1;j<Ny;j++){
            for(int k=0;k<Nz;k++)v.value[i][j][k] = v.value[i][j][k] - dt/rho.value[i][j][k] * (p.value[i][j][k]-p.value[i][j-1][k])/dx;
        }
    }
    for(int i=0;i<Nx;i++){
        for(int j=0;j<Ny;j++){
            for(int k=1;k<Nz;k++)w.value[i][j][k] = w.value[i][j][k] - dt/rho.value[i][j][k] * (p.value[i][j][k]-p.value[i][j][k-1])/dx;
        }
    }
}
void cal_buoyanacy(){
    Eigen::Vector3d dir_gravity = {0.0,1.0,0.0};
    for(int i=0;i<Nx;i++){
        for(int j=0;j<Ny;j++){
            for(int k=0;k<Nz;k++){
                
            }
        }
    }
}
