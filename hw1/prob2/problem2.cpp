#include <iostream>
#include <fstream>
#include <complex>
#include <vector>
#include <cstdio>
using namespace std;



uint32_t mbseries(complex<double> c){
    
    // declaring max iterations and the complex number of interest
    std::complex<double> z(0,0);
    const int max_iter = 10000;

    for(int i = 0; i < max_iter; i++){
    
        z = z*z + c;

        if(abs(z) > 2){
            return i;
        } 

    }

    return 0;
}


void bit_mesh(int xpts,int ypts, vector<uint32_t> &mesh ){

    //declare system bounds
    const short int xmin = -2;
    const short int xmax = 2;
    const short int ymin = -1;
    const short int ymax = 1;

    //expand vector size to fit the complete dataset
    const int mesh_size = xpts * ypts;
    mesh.reserve(mesh_size);
    
    //determine increment size and declare starting point
    double xinc = 1.0*(xmax - xmin)/(xpts - 1.0);
    double yinc = 1.0*(ymax - ymin)/(ypts - 1.0);

    double cur_y = -1.0;
    double cur_x = -2.0;
    complex<double> c(cur_x,cur_y);

    //loop through the mesh populating line by line in the x-direction
    for(int iy = 0; iy < ypts; ++iy){

        for(int ix = 0; ix < xpts; ++ix){

            mesh.push_back(mbseries(c));
            c = complex<double>(c.real() + xinc, c.imag());
        }
        c = complex<double>(xmin, c.imag()+yinc);
    }
}


void area_calc(vector<uint32_t> &mesh ){

    const int plot_range_area = 8;
    uint32_t pix_tot = mesh.size();
    uint32_t pix_converged = 0;
    double area_pix;
    double mb_area;

    for(const auto &i : mesh){
         
         if(i==0){
            pix_converged+=1;
        }
    } 

    area_pix = (1.0 * pix_converged)/ pix_tot;
    mb_area = 1.0 * plot_range_area * area_pix;

    printf("Percentage of Pixels that converged = %.3f\n", area_pix);
    printf("Area of Mandelbrot set = %.3f\n", mb_area);

}


int main(){

    // instantiate vectors that will be used
    vector<uint32_t> low_res_mesh;
    vector<uint32_t> mid_res_mesh;
    vector<uint32_t> high_res_mesh;


    // Call functions create mesh and output mesh into a textfile
    // Low resolution 10 x 10
    bit_mesh(10,10, low_res_mesh);
    ofstream low_fout("LOW_res.txt");
    for(const auto &i : low_res_mesh){
        low_fout << i << "\n";
    } 
    cout << "\nlow resolution mesh size = ";
    cout << low_res_mesh.size() << endl;

    area_calc(low_res_mesh);

    // Mid resolution 100 x 100
    bit_mesh(100,100, mid_res_mesh);
    ofstream mid_fout("MID_res.txt");
    for(const auto &i : mid_res_mesh){
        mid_fout << i << "\n";
    } 
    cout << "\nmid resolution mesh size = ";
    cout << mid_res_mesh.size() << endl;
    
    area_calc(mid_res_mesh);

    // High resolution 1000 x 1000
    bit_mesh(1000,1000, high_res_mesh);
    ofstream high_fout("HIGH_res.txt");
    for(const auto &i : high_res_mesh){
        high_fout << i << "\n";
    } 
    cout << "\nhigh resolution mesh size = ";
    cout << high_res_mesh.size() << endl;

    area_calc(high_res_mesh);
   
    return 0;
}


