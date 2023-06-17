//
//  Array3d.hpp
//
//

#ifndef Array3d_hpp
#define Array3d_hpp

#include <vector>
#include <algorithm>
#include <string>
#include <sstream>
#include <Eigen/Core>
#include <iostream>
struct myArray3d{
    std::vector<std::vector<std::vector<double>>>value;
    int nx,ny,nz;
    int size;
    myArray3d();
    myArray3d(int size_x,int size_y,int size_z);
    myArray3d(int size_x,int size_y,int size_z,double val);
    void reset(double val);
    void print();
    std::vector<double> convert2Vector();
};
struct myMap{
    std::vector<std::vector<std::vector<std::vector<int>>>>value;
    int nx,ny,nz;
    myMap();
    myMap(int size_x,int size_y,int size_z);
    void reset();
    void print();
    bool contains(std::vector<int> &key);
    std::vector<int> at(std::vector<int> &key);
    std::vector<double> convert2Vector();
};
#endif /* Array3d_hpp */
