#include <iostream>
#include <iomanip>
#include <complex>
#include <vector>
using namespace std;



bool mbseries(complex<double> c){
    
    std::complex<double> z = 0;
    const int max_iter = 100

    for(int i = 0; i < max_iter; i++){
        z = z*z + c;
        if(std::abs(z) >= 2){
            return false;
        } 
    }

    if(std::abs(z) <= 2){
            return true;
    } 

}


bool * mandlebrot_set(int xpts,int ypts ){


    const short int xmin = -2;
    const short int xmax = 2;
    const short int ymin = -1;
    const short int ymax = 1;
    const int mesh_size = xpts * ypts;
    double xinc = (xmax - xmin)/(xpts - 1)
    double yinc = (ymax - ymin)/(ypts - 1)

    mb_mesh.resize(mesh_size);

    complex<double> c;

    for(int iy = 0; iy < ypts; iy++){
        for(int ix = 0; ix < xpts; ix++){

            mbseries(c)
        }
    }
}

int main{

    vector<bool> mb_mesh;
    return 0;
}


