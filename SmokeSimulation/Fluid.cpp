//
//  Fluid.cpp












//力，密度，温度の初期化(Template Array3<double>)，温度と渦の力計算，計算のループ作成，密度を使ってテキトーな光の計算，平行投影のなんちゃって横断アルゴリズム




















#include "Fluid.hpp"
Fluid::Fluid(){
    dx = 1;
    dt = 1;
}
Fluid::Fluid(double x,double t,myArray3<double> &density,myArray3<double> &templature){
    dx = x;
    dt = t;
    rho = density;
    temp = templature;
}
void Fluid::setDensity(){
    for(int i=0;i<Nx;i++){
        for(int j=0;j<Ny;j++){
            for(int k=0;k<Nz;k++){
                if(Nx/3 < i && i < Nx/3*2 && k == 0 && Nz/3 < j && j < Nz/3*2)rho.value[i][j][k] = 2.0;
                else rho.value[i][j][k] = 1.0;
            }
        }
    }
}
void Fluid::setTemplature(){
    for(int i=0;i<Nx;i++){
        for(int j=0;j<Ny;j++){
            for(int k=0;k<Nz;k++){
                if(Nx/3 < i && i < Nx/3*2 && k == 0 && Nz/3 < j && j < Nz/3*2)temp.value[i][j][k] = 100;
            }
        }
    }
}
double Fluid::TriLinearInterporation(double x,double y,double z,myArray3<double> &val){
//    std::cout << "input_xyz = (" << x << "," << y << "," << z << ")" << std::endl;
//    std::cout << "x=fmax(0.0, fmin(" << val.nx-1-1e-6 << "," << x/dx << "));" << std::endl;
    x = fmax(0.0, fmin(val.nx-1-1e-6,x/dx));
//    std::cout << "x=fmax(0.0, fmin(" << val.ny-1-1e-6 << "," << y/dx << "));"<< std::endl;
    y = fmax(0.0, fmin(val.ny-1-1e-6,y/dx));
//    std::cout << "x=fmax(0.0, fmin(" << val.nz-1-1e-6 << "," << z/dx << "));"<< std::endl;
    z = fmax(0.0, fmin(val.nz-1-1e-6,z/dx));
//    std::cout << "fix_xyz = (" << x << "," << y << "," << z << ")" << std::endl;
    int i = x;int j = y;int k = z;
//    std::cout << "index_ijk = (" << i << "," << j << "," << k << ")" << std::endl;
    double s = x-i;double t = y-j;double u = z-k;
//    std::cout << "ratio_stu = (" << s << "," << t << "," << u << ")" << std::endl;
    Eigen::Vector<double,8> f = {
        val.value[i][j][k],val.value[i+1][j][k],val.value[i][j+1][k],val.value[i][j][k+1],
        val.value[i+1][j+1][k],val.value[i][j+1][k+1],val.value[i+1][j][k+1],val.value[i+1][j+1][k+1]};
    Eigen::Vector<double,8> c = {
        (1-s)*(1-t)*(1-u),s*(1-t)*(1-u),s*t*(1-u),(1-s)*t*(1-u),
        (1-s)*(1-t)*u,s*(1-t)*u,s*t*u,(1-s)*t*u};
//    std::cout << f[0] << "," << f[1] << "," << f[2] << "," << f[3] << std::endl;
//    std::cout << f[4] << "," << f[5] << "," << f[6] << "," << f[7] << std::endl;
//    std::cout << c[0] << "," << c[1] << "," << c[2] << "," << c[3] << std::endl;
//    std::cout << c[4] << "," << c[5] << "," << c[6] << "," << c[7] << std::endl;
//    std::cout << u <<"( ("<< t << "( (" << 1-s << ")"<<f[0] << "+" << s <<f[1] << ")+(" << 1-t << "( (" << 1-s << ")" << f[2] << "+" << s << f[4] << ") )" << std::endl;
//    for(int i=0;i<8;++i)std::cout << f[i]*c[i] << std::endl;
    return f.dot(c);
}
void Fluid::faceAdvect(){
    //x
    old_u = u;
    old_v = v;
    old_w = w;
    for(int i=1;i<u.nx;++i){
        for(int j=0;j<u.ny;++j){
            for(int k=0;k<u.nz;++k){
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
            for(int k=0;k<v.nz;++k){
                double x = (i+0.5)*dx;double y = j*dx;double z = (k+0.5)*dx;
                double adv_x = x - dt*TriLinearInterporation(x, y-0.5*dx, z-0.5*dx, old_u);
                double adv_y = y - dt*TriLinearInterporation(x-0.5*dx, y, z-0.5*dx, old_v);
                double adv_z = z - dt*TriLinearInterporation(x-0.5*dx, y-0.5*dx, z, old_w);
                v.value[i][j][k] = TriLinearInterporation(adv_x, adv_y, adv_z, old_v);
                //if(Nx/3 < i && i < Nx/3*2 && j < 3 && Nz/3 < k && k < Nz/3*2)
                //std::cout <<"(" << i << "," << j << "," << k << ")=(" << adv_x/dx << "," << adv_y/dx <<"," << adv_z/dx <<")"<< std::endl;
            }
        }
    }
    for(int i=0;i<w.nx;++i){
        for(int j=0;j<w.ny;++j){
            for(int k=1;k<w.nz;++k){
                double x = (i+0.5)*dx;double y = (j+0.5)*dx;double z = k*dx;
                double adv_x = x - dt*TriLinearInterporation(x, y-0.5*dx, z-0.5*dx, old_u);
                double adv_y = y - dt*TriLinearInterporation(x-0.5*dx, y, z-0.5*dx, old_v);
                double adv_z = z - dt*TriLinearInterporation(x-0.5*dx, y-0.5*dx, z, old_w);
                w.value[i][j][k] = TriLinearInterporation(adv_x, adv_y, adv_z, old_w);
            }
        }
    }
}
void Fluid::centerAdvect(myArray3<double> &val){
    myArray3<double> old_val = val;
    for(int i=0;i<val.nx;++i){
        for(int j=0;j<val.ny;++j){
            for(int k=0;k<val.nz;++k){
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
Eigen::Vector3d Fluid::getBuoyanacy(int i,int j, int k){
    Eigen::Vector3d dir_gravity = {0.0,1.0,0.0};
    return dx*(-g0*rho.value[i][j][k] + beta*(temp.value[i][j][k] - Tamb))*dir_gravity;
}

void Fluid::setCenterRot(){
    for(int i=1;i<Nx-1;i++){
        for(int j=1;j<Ny-1;j++){
            for(int k=1;k<Nz-1;k++){
                double u_P = (u.value[i-1][j][k+1] + u.value[i][j][k])/2;
                double u_M = (u.value[i-1][j][k] + u.value[i][j][k-1])/2;
                double v_P = (v.value[i+1][j-1][k] + v.value[i][j][k])/2;
                double v_M = (v.value[i][j-1][k] + v.value[i-1][j][k])/2;
                double w_P = (w.value[i][j+1][k-1] + w.value[i][j][k])/2;
                double w_M = (w.value[i][j][k-1] + w.value[i][j-1][k])/2;
                double x = (w_P - w_M - v_P+ v_M)/2*dx;
                double y = (u_P - u_M - w_P + w_M)/2*dx;
                double z = (v_P - v_M - u_P + u_M)/2*dx;
                centerRot.value[i][j][k] = {x,y,z};
            }
        }
    }
}
Eigen::Vector3d Fluid::getConfinement(int i,int j,int k){
    if(i==0||j==0||k==0)std::cout << "error_0" << std::endl;
    if(i==Nx-1||j==Ny-1||k==Nz-1)std::cout << "error_Nx" << std::endl;
    double u_P = centerRot.value[i][j][k+1].x();
    double u_M = centerRot.value[i][j][k-1].x();
    double v_P = centerRot.value[i+1][j][k].y();
    double v_M = centerRot.value[i-1][j][k].y();
    double w_P = centerRot.value[i][j+1][k].z();
    double w_M = centerRot.value[i][j-1][k].z();
    double x = (w_P - w_M - v_P+ v_M)/2*dx;
    double y = (u_P - u_M - w_P + w_M)/2*dx;
    double z = (v_P - v_M - u_P + u_M)/2*dx;
    double scale = epcilon*dx*centerRot.value[i][j][k].norm();
    return {x*scale,y*scale,z*scale};
}
void Fluid::addForce(){
    setCenterRot();
    for(int i=0;i<Nx;i++){
        for(int j=0;j<Ny;j++){
            for(int k=0;k<Nz;k++){
                f.value[i][j][k] = getBuoyanacy(i, j, k);
            }
        }
    }
    for(int i=1;i<Nx-1;i++){
        for(int j=1;j<Ny-1;j++){
            for(int k=1;k<Nz-1;k++){
                //f.value[i][j][k] += getConfinement(i, j, k);
            }
        }
    }
    //x
    for(int i=1;i<f.nx;i++){
        for(int j=0;j<f.ny;j++){
            for(int k=0;k<f.nz;k++){
                u.value[i][j][k] += (f.value[i-1][j][k].x() + f.value[i][j][k].x())/2;
            }
        }
    }
    //y
    for(int i=0;i<f.nx;i++){
        for(int j=1;j<f.ny;j++){
            for(int k=0;k<f.nz;k++){
                v.value[i][j][k] += (f.value[i][j-1][k].y() + f.value[i][j][k].x())/2;
            }
        }
    }
    //z
    for(int i=0;i<f.nx;i++){
        for(int j=0;j<f.ny;j++){
            for(int k=1;k<f.nz;k++){
                w.value[i][j][k] += (f.value[i][j][k-1].x() + f.value[i][j][k].x())/2;
            }
        }
    }
}

void Fluid::oneloop(){
    addForce();
    std::cout << "addForce" << std::endl;
    faceAdvect();
    //v.print();
    std::cout << "faceAdvect" << std::endl;
    project();
    std::cout << "project" << std::endl;
    centerAdvect(temp);
    std::cout << "centerAdvectTemp" << std::endl;
    centerAdvect(rho);
    std::cout << "centerAdvectRho" << std::endl;
}
void Fluid::execute(){
    dx = 0.01;
    dt = 0.001;
    setDensity();
    setTemplature();
    std::cout << "Initialize" << std::endl;
    std::string rootFolderName = "Result";
    std::string pressureFolderName = rootFolderName + "/pressure";
    std::string densityFolderName = rootFolderName + "/density";
    std::string templatureFolderName = rootFolderName + "/templature";
    std::filesystem::create_directories(rootFolderName);
    std::filesystem::create_directories(pressureFolderName);
    std::filesystem::create_directories(densityFolderName);
    std::filesystem::create_directories(templatureFolderName);
    for(int i=0;i<timestep;++i){
        p.print();
        oneloop();
        std::string OutputVTK_pre = pressureFolderName+  "/output"+std::to_string(i)+".vtk";
        std::string OutputVTK_den = densityFolderName+  "/output"+std::to_string(i)+".vtk";
        std::string OutputVTK_tem = templatureFolderName+  "/output"+std::to_string(i)+".vtk";
        outputVTK(OutputVTK_pre,p,dx);
        outputVTK(OutputVTK_den,rho,dx);
        outputVTK(OutputVTK_tem,temp,dx);
    }
}
