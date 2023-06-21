//
//  rendering.hpp
//  SmokeSimulation
//
//  Created by 須之内俊樹 on 2023/06/20.
//
#ifndef rendering_hpp
#define rendering_hpp
#define albedo 0.5

#define ambient 0.9
#define threthold 1.0
#define background 255
#include "Array3d.hpp"
#include "Eigen/Core"
#include <cstddef>
#include <memory>
#include <new>
#include <random>
#include "stb_image_write.h"
struct RGBA {
    unsigned char r, g, b, a; //赤, 緑, 青, 透過
    RGBA() = default;
    constexpr RGBA(const unsigned char r_, const unsigned char g_, const unsigned char b_, const unsigned char a_) :r(r_), g(g_), b(b_), a(a_) {}
};
void rayshoot(myArray3<double> &rho);
void cal_voxTrans(myArray3<double> &rho,myArray3<double> &vox_Trans);
bool generateImage(std::string &imageName,myArray3<double> &vox_Trance,myArray3<double> &rho);
#endif /* rendering_hpp */
