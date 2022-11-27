#include <iostream>
#include <fstream>
#include <complex>
#include <vector>
#include <cstdio>
#include "mpi.h"
#include <chrono>



uint32_t mbseries(std::complex<double> c){
    
    // declaring max iterations and the complex number of interest
    std::complex<double> z(0,0);
    const int max_iter = 1000;

    for(int i = 0; i < max_iter; i++){
    
        z = z*z + c;

        if(abs(z) > 2){
            return i+1;
        } 

    }

    return 0;
}


void bit_mesh(int xpts_global,int ypts_global, std::vector<uint32_t> &total_mesh ){

    int rank;
    int worldSize;

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &worldSize);

    //declare global system bounds
    const short int xmin_global = -2;
    const short int xmax_global = 2;
    const short int ymin_global = -1;
    const short int ymax_global = 1;
        
    //determine increment size and declare starting point
    double xinc = 1.0*(xmax_global - xmin_global)/(xpts_global - 1.0);
    double yinc = 1.0*(ymax_global - ymin_global)/(ypts_global- 1.0);

    int xpts_loc;
    double xmin_loc;

    int root_xpts = xpts_global/worldSize + xpts_global%worldSize;
    int non_root_xpts = xpts_global/worldSize;

    if(rank == 0 ){
        xpts_loc = root_xpts;
        xmin_loc = 1.0*xmin_global ;

        //int ypts_loc = ypts_global/worldSize + ypts_global%worldSize;
        
    }
    else{
        xpts_loc = non_root_xpts;
        xmin_loc = 1.0*xmin_global + xinc*(1.0*xpts_loc*rank 
                                      + 1.0*(xpts_global%worldSize));

        //int ypts_loc = ypts_global/worldSize;
    }

    //expand vector size to fit the complete dataset
    const int total_mesh_size = xpts_global * ypts_global;
    total_mesh.reserve(total_mesh_size);

    int mesh_size_loc = xpts_loc*ypts_global;
    int root_mesh_size = root_xpts*ypts_global;
    int non_root_mesh_size = non_root_xpts*ypts_global;

    std::vector<uint32_t> holder;
    holder.reserve(mesh_size_loc);

    double cur_y = -1.0;
    double cur_x = xmin_loc;
    std::complex<double> c(cur_x,cur_y);

    //loop through the mesh populating line by line in the x-direction
    for(int iy = 0; iy < ypts_global; ++iy){

        for(int ix = 0; ix < xpts_loc; ++ix){

            holder.insert(holder.begin() + ix,mbseries(c));
            //holder.push_back(mbseries(c));
            c = std::complex<double>(c.real() + xinc, c.imag());
        }
        std::printf("Rank %d:\tholder.size()=%ld\n",rank, holder.size());
    MPI_Barrier(MPI_COMM_WORLD);

    if(rank==0){

        total_mesh.resize(total_mesh_size);
        int counts[worldSize];
        std::fill_n(counts,worldSize, non_root_mesh_size);
        counts[0] = root_mesh_size;

        int displacements[worldSize];
        displacements[0] = 0;

        switch (worldSize)
        {
        case 1:
            break;
        case 2:
            displacements[1] = root_mesh_size;
        default:            
            displacements[1] = root_mesh_size;
            for(int i = 2; i < worldSize; i++){
                displacements[i] = ypts_global * ( root_xpts + non_root_xpts*(i-1) );
            }
        }
        //xpts_loc*rank + (xpts_global%worldSize)
    
        MPI_Gatherv(holder.data(), mesh_size_loc, MPI_UINT32_T, total_mesh.data(), counts, displacements, MPI_UINT32_T, 0, MPI_COMM_WORLD);    
    }
    else{
        MPI_Gatherv(holder.data(), mesh_size_loc, MPI_UINT32_T, NULL, NULL,NULL, MPI_UINT32_T, 0,MPI_COMM_WORLD);  
    }

        holder.clear();
        holder.reserve(xpts_loc);    
        c = std::complex<double>(xmin_loc, c.imag()+yinc);
    }

    
}


void area_calc(std::vector<uint32_t> &mesh ){

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


int main(int argc, char *argv[] ){

    MPI_Init(&argc, &argv);

    int rank;
    int worldSize;

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &worldSize);
/*
    double t_init, t_final, t_elapsed;

    if(rank == 0 ){
    t_init = MPI_Wtime();
    }
*/
    // instantiate vectors that will be used
    std::vector<uint32_t> rank_mesh, total_mesh;

    // Call functions create mesh and output mesh into a textfile

    // High resolution 1000 x 1000
    bit_mesh(10,10, total_mesh);

    if(rank == 0 ){
        std::ofstream high_fout("HIGH_res.txt");
        for(const auto &i : total_mesh){
            high_fout << i << "\n";
        } 
        std::cout << "\nhigh resolution mesh size = ";
        std::cout << total_mesh.size() << std::endl;

        //area_calc(rank_mesh);
    }
    MPI_Finalize();
    return 0;
}


