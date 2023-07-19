//
//  output.hpp
//  SmokeSimulation
//
//  Created by 須之内俊樹 on 2023/06/18.
//

#ifndef output_hpp
#define output_hpp

#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <iomanip>
#include <Eigen/Core>
#include "Array3d.hpp"
#include <chrono>
struct TimeDisplayer{
    std::chrono::system_clock::time_point startTime;
    std::chrono::system_clock::time_point endTime;
    const char* str;
    void startTimer(const char* s);
    double endTimer();
};

void outputVTK(std::string OutputFileName,myArray3<double> &val,double dx);
#endif /* output_hpp */
