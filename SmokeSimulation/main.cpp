//
//  main.cpp
//  SmokeSimulation
//
//  Created by 須之内俊樹 on 2023/06/13.
//

#include <iostream>
#include "Fluid.hpp"
#include "Array3d.hpp"
int main(int argc, const char * argv[]) {
//    double dx = 0.01;
//    double dt = 0.1;
//    myArray3<double> density;
//    myArray3<double> templature;
//    Fluid simulator = Fluid(dx,dt,density,templature);
    Fluid simulator = Fluid();
    simulator.execute();
    //myArray3<double> test = myArray3<double>(2,2,2,0);
}
