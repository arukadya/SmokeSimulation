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
template <typename T>
struct myArray3{
    std::vector<std::vector<std::vector<T>>>value;
    int nx,ny,nz;
    int size;
    myArray3();
//    myArray3(int size_x,int size_y,int size_z);
    myArray3(int size_x,int size_y,int size_z,T val);
    void reset(T val);
    void print();
    //std::vector<T> convert2Vector();
};
//struct myMap{
//    std::vector<std::vector<std::vector<std::vector<int>>>>value;
//    int nx,ny,nz;
//    myMap();
//    myMap(int size_x,int size_y,int size_z);
//    void reset();
//    void print();
//    bool contains(std::vector<int> &key);
//    std::vector<int> at(std::vector<int> &key);
//    std::vector<double> convert2Vector();
//};
template <typename T>
myArray3<T>::myArray3(){}
template <typename T>
myArray3<T>::myArray3(int size_x,int size_y,int size_z,T val){
    nx = size_x;
    ny = size_y;
    nz = size_z;
    size = nx*ny*nz;
    value.resize(nx);
    for(int i=0;i<nx;i++){
        value[i].resize(ny);
        for(int j=0;j<ny;j++){
            value[i][j].resize(nz);
            for(int k=0;k<nz;k++)value[i][j][k] = val;
        }
    }
}
template <typename T>
void myArray3<T>::reset(T val){
    for(int i=0;i<nx;i++){
        for(int j=0;j<ny;j++){
            for(int k=0;k<nz;k++)value[i][j][k] = val;
        }
    }
}
template <typename T>
void myArray3<T>::print(){
    for(int i=0;i<nx;i++){
        for(int j=0;j<ny;j++){
            for(int k=0;k<nz;k++)std::cout << value[i][j][k] << ",";
            std::cout << std::endl;
        }std::cout << std::endl;
    }std::cout << std::endl;
}
//template <>
//std::vector<double> myArray3<double>::convert2Vector(){
//    std::vector<double> data(nx*ny*nz);
//    for(int i=0;i<nx;++i){
//        for(int j=0;j<ny;++j){
//            for(int k=0;k<nz;++k){
//                data[i + j*nx + k*nx*ny] = value[i][j][k];
//            }
//        }
//    }
//    return data;
//}
//myMap::myMap(){}
//myMap::myMap(int size_x,int size_y,int size_z){
//    nx = size_x;
//    ny = size_y;
//    nz = size_z;
//    std::vector<int> zero = {};
//    value.resize(size_x);
//    //size.resize(size_x);
//    for(int i=0;i<size_x;++i){
//        value[i].resize(size_y);
//        //size[i].resize(size_y);
//        for(int j=0;j<size_y;j++){
//            value[i][j].resize(size_z);
//            //size[i][j].resize(size_z);
//            for(int k=0;k<size_z;k++){
//                value[i][j][k] = zero;
//                //size[i][j][k] = 0;
//            }
//        }
//    }
//}
//void myMap::reset(){
//    for(int i=0;i<nx;i++){
//        for(int j=0;j<ny;j++){
//            for(int k=0;k<nz;k++)value[i][j][k] = std::vector<int>();
//        }
//    }
//}
//void myMap::print(){
//    for(int i=0;i<nx;i++){
//        for(int j=0;j<ny;j++){
//            for(int k=0;k<nz;k++){
//                if(!value[i][j][k].empty())std::cout << value[i][j][k].size() << " ";
//                else std::cout << 0 << " ";
//            }
//            std::cout << std::endl;
//        }
//        std::cout << std::endl;
//    }
//    std::cout << std::endl;
//}
//bool myMap::contains(std::vector<int> &key){
//    if(value[key[0]][key[1]][key[2]].empty())return false;
//    else return true;
//}
//std::vector<int> myMap::at(std::vector<int> &key){
//    return value[key[0]][key[1]][key[2]];
//}

#endif /* Array3d_hpp */
