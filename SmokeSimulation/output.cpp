//
//  output.cpp
//  SmokeSimulation
//
//  Created by 須之内俊樹 on 2023/06/18.
//

#include "output.hpp"
void outputVTK(std::string OutputFileName,myArray3<double> &val,double dx){
    std::vector<float> origin = {0,0,0};
    int nx = val.nx;
    int ny = val.ny;
    int nz = val.nz;
    std::ofstream writing_file;
    std::cout << OutputFileName << std::endl;
    writing_file.open(OutputFileName, std::ios::out);
    if(writing_file.is_open())std::cout << "success!!" << std::endl;
    std::string writing_text ="# vtk DataFile Version 2.0\nIsosurface\nASCII\nDATASET STRUCTURED_POINTS\n";
    writing_file << writing_text << std::endl;
    writing_file << "DIMENSIONS " << nx <<" "<< ny <<" "<< nz << std::endl;
    writing_file << "ORIGIN " << origin[0] <<" "<< origin[1] <<" "<< origin[2] << std::endl;
    writing_file << "SPACING " << dx <<" "<< dx <<" "<< dx << std::endl;
    writing_file << "POINT_DATA " << val.size << std::endl;
    writing_file << "SCALARS value float 1" << std::endl;
    writing_file << "LOOKUP_TABLE default" << std::endl;
//    for(int i=0;i<nx;i++){
//        for(int j=0;j<ny;j++){
//            for(int k=0;k<nz;k++){
//                writing_file << val.value[i][j][k] << std::endl;
//            }
//        }
//    }
    for(int k=0;k<nz;k++){
        for(int j=0;j<ny;j++){
            for(int i=0;i<nx;i++){
                writing_file << val.value[i][j][k] << std::endl;
            }
        }
    }
    writing_file.close();
}
