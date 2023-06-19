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
    Fluid simulator = Fluid();
    simulator.execute();
//    myArray3<double> test = myArray3<double>(2,2,3,0);
//    test.value[0][0][0] = 10;
//    test.value[1][0][0] = 10;
//    test.value[1][1][0] = 10;
//    test.value[0][1][0] = 10;
//    test.value[0][0][1] = 100;
//    test.value[1][0][1] = 100;
//    test.value[1][1][1] = 100;
//    test.value[0][1][1] = 100;
//    test.value[0][0][2] = 100;
//    test.value[1][0][2] = 100;
//    test.value[1][1][2] = 100;
//    test.value[0][1][2] = 100;
//
//    std::cout << simulator.TriLinearInterporation(0.6, 0.2, -1, test) << std::endl;
}
