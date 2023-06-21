//
//  rendering.cpp
//  SmokeSimulation
//
//  Created by 須之内俊樹 on 2023/06/20.
//

#include "rendering.hpp"
void cal_voxTrans(myArray3<double> &rho,myArray3<double> &vox_Trans){
    for(int ix=0;ix<rho.nx;++ix){
        for(int iz=0;iz<rho.nz;++iz){
            double ray_trans = 1;
            for(int iy=0;iy<rho.ny;++iy){
                double extinction = 1.0/rho.value[ix][iy][iz];
                double vox_trans = exp(-extinction);
                vox_Trans.value[ix][iy][iz] = albedo*ambient*(1-vox_trans)*ray_trans;
                if(rho.value[ix][iy][iz] < threthold)vox_Trans.value[ix][iy][iz] *=2;
                ray_trans = ray_trans*vox_trans;
            }
        }
    }
}
bool generateImage(std::string &imageName,myArray3<double> &vox_Trance,myArray3<double> &rho){
    
    constexpr std::size_t width{ 256 }, height{ 256 }; //幅と高さ
    
    std::unique_ptr<RGBA[][width]> rgba(new(std::nothrow) RGBA[height][width]);
    if (!rgba) return false;

    //std::random_device rd;
    //std::mt19937 mt;
    //mt.seed(rd());
    
    std::uniform_int_distribution<> uid(0, 255);
    for (std::size_t row{}; row < height; ++row)
        for (std::size_t col{}; col < width; ++col) {
            int x = col - width/2;
            int z = row - height/5;
            rgba[row][col].r = background; //赤
            rgba[row][col].g = background;
            rgba[row][col].b = background;
            rgba[row][col].a = 255; //不透過
            if( 0<=x && x<vox_Trance.nx && 0<=z && z<vox_Trance.nz){
                double trans = 0;
                for(int iy = 0;iy < vox_Trance.ny;++iy){
                    //if(rho.value[x][iy][z] > threthold){
                        //rgba[row][col].a = 0;
                    //}
                    trans += 255*vox_Trance.value[x][iy][z];
                    //std::cout << "r+=" << 255*vox_Trance.value[x][iy][z] << ",r=" << red << std::endl;
                }
                rgba[row][col].r = 200; //赤
                rgba[row][col].g = 0;
                rgba[row][col].b = 0;
                rgba[row][col].a = (int)trans;
                //std::cout << std::endl;
            }
        }
    stbi_write_png(imageName.c_str(), static_cast<int>(width), static_cast<int>(height), static_cast<int>(sizeof(RGBA)), rgba.get(), 0);
    return true;
}
