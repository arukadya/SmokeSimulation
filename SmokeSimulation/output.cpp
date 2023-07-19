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
    std::string writing_text ="# vtk DataFile Version 2.0\nIsosurface\nASCII\nDATASET STRUCTURED_POINTS\n";
    writing_file << writing_text << std::endl;
    writing_file << "DIMENSIONS " << nx <<" "<< ny <<" "<< nz << std::endl;
    writing_file << "ORIGIN " << origin[0] <<" "<< origin[1] <<" "<< origin[2] << std::endl;
    writing_file << "SPACING " << dx <<" "<< dx <<" "<< dx << std::endl;
    writing_file << "POINT_DATA " << val.size << std::endl;
    writing_file << "SCALARS value float 1" << std::endl;
    writing_file << "LOOKUP_TABLE default" << std::endl;
    for(int k=0;k<nz;k++){
        for(int j=0;j<ny;j++){
            for(int i=0;i<nx;i++){
                if(abs(val.value[i][j][k]) < 1e-10) writing_file << 0 << std::endl;
                else writing_file << val.value[i][j][k] << std::endl;
            }
        }
    }
    writing_file.close();
}
void TimeDisplayer::startTimer(const char* s){
    startTime = std::chrono::system_clock::now();
    str = s;
}
double TimeDisplayer::endTimer(){
    endTime = std::chrono::system_clock::now();
    double time = static_cast<double>(std::chrono::duration_cast<std::chrono::milliseconds>(endTime - startTime).count());
    std::cout << str << ":" << time << "ms" << std::endl;
    return time;
    //std::cout << std::endl;
}
//void outputPLT(std::string OutputFileName ,std::string pltFileName,std::vector<double> &data){
//    //vertex
//    std::ofstream writing_file;
//    std::filesystem::create_directories("outputData");
//    writing_file.open("outputData/" + OutputFileName + "_v.dat",std::ios::out);
//    writing_file.close();
//    writing_file.open(pltFileName , std::ios::out);
//    std::filesystem::create_directories("protImage");
//    writing_file << "plot " << "'outputData/" + OutputFileName + "_v.dat'"<<" w linespoint pointsize 2 pointtype 4" << std::endl;
//    writing_file << "set terminal jpeg" << std::endl;
//    writing_file << "set output " << "'outputImage/" + OutputFileName +".jpeg'" << std::endl;
//    writing_file << "replot" << std::endl;
//    writing_file.close();
//
//}
