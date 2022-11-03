#include <iostream>
#include <fstream>
#include <complex>
#include <vector>
#include <cstdio>
using namespace std;



uint32_t mbseries(complex<double> c){
    
    std::complex<double> z(0,0);
    const int max_iter = 4000;

    for(int i = 0; i < max_iter; i++){
    
        z = z*z + c;

        if(abs(z) > 2){
            return i;
        } 

    }

    return 0;
}


void bit_mesh(int xpts,int ypts, vector<uint32_t> &mesh ){


    const short int xmin = -2;
    const short int xmax = 2;
    const short int ymin = -1;
    const short int ymax = 1;

    const int mesh_size = xpts * ypts;
    mesh.reserve(mesh_size);
    
    double xinc = 1.0*(xmax - xmin)/(xpts - 1.0);
    double yinc = 1.0*(ymax - ymin)/(ypts - 1.0);

    double cur_y = -1.0;
    double cur_x = -2.0;
    complex<double> c(cur_x,cur_y);

    for(int iy = 0; iy < ypts; ++iy){

        for(int ix = 0; ix < xpts; ++ix){

            mesh.push_back(mbseries(c));
            c = complex<double>(c.real() + xinc, c.imag());
        }
        c = complex<double>(xmin, c.imag()+yinc);
    }
}

int main(){

    vector<uint32_t> ultra_res_mesh;
    //vector<uint32_t> ultra_res_mesh;

    // High resolution 10k x 10k
    bit_mesh(10000, 10000, ultra_res_mesh);
    ofstream ultra_fout("ultra_res_mesh.txt");
    for(const auto &i : ultra_res_mesh){
        ultra_fout << i << "\n";
    } 
    cout << "\nUltra high resolution mesh size  = ";
    cout << ultra_res_mesh.size() << endl;
   
    return 0;
}


