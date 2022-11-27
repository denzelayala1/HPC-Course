#include <iostream>
#include <fstream>
#include <complex>
#include <vector>
#include <cstdio>
#include "mpi.h"
#include <chrono>
#include <stdlib.h>

uint32_t mbseries(std::complex<double> c);
void area_calc(std::vector<uint32_t> &mesh );
void mesh_unifier(std::vector<uint32_t> &output_mesh,std::vector<uint32_t> &holder);
void plot_maker(int xpts_global,int ypts_global, std::vector<uint32_t> &output_mesh );
void export_mesh(std::vector<uint32_t> &mesh);

int main(int argc, char *argv[] ){

    MPI_Init(&argc, &argv);

    int rank;
    int worldSize;

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &worldSize);

    double t_init, t_final, t_elapsed;

    //if(rank == 0 ){
    t_init = MPI_Wtime();
    

    // instantiate vectors that will be used
    std::vector<uint32_t> output_mesh;

    // Call functions create mesh and output mesh into a textfile

    // High resolution 1000 x 1000
    plot_maker(1000,1000, output_mesh);

    t_final = MPI_Wtime();
    t_elapsed = t_final - t_init;
    printf("Rank %d\tElapsed Time: %.10f seconds\n", rank, t_elapsed);

    if(rank == 0 ){
        export_mesh(output_mesh);
        area_calc(output_mesh);
    }

    MPI_Finalize();
    return 0;
}


/*

 Populates a vector with the Mandelbrot set in the domain -2<x<2 & -1<y<1 

 The function will find the step size from x and y resolution. From there it calculates values starting from the bottom left corner. It first calculates an entire row (x-values) then increments once in the y-direction and calculates the next row. The output vector is populated every row.


 @param xpts_global the resolution in the x-direction
 @param ypts_global the resolution in the y-direction
 @param output_mesh the vector that is being populated
*/
void plot_maker(int xpts_global,int ypts_global, std::vector<uint32_t> &output_mesh ){

    //declare global system bounds

    const short int xmin_global = -2;
    const short int xmax_global = 2;
    const short int ymin_global = -1;
    const short int ymax_global = 1;

    // x increment size. same for all ranks
    double xinc = 1.0*(xmax_global - xmin_global)/(xpts_global - 1.0);
    // y increment size. same for all ranks
    double yinc = 1.0*(ymax_global - ymin_global)/(ypts_global- 1.0);

    //MPI variables

    int rank;    //  MPI rank of this core
    int worldSize;    //  MPI number of cores used

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &worldSize);

    
    const int non_root_xpts = xpts_global/worldSize; // the number of points all non-root cores will calculate
    const int remainder = xpts_global%worldSize; // The left over points from floor division
    const int root_xpts = non_root_xpts + remainder;// Non-root points plus the remainder  

    int xpts_loc;     // Number of x-direction points THIS core will calculate
    double xmin_loc;    // The "starting" x-point for THIS core
    if(rank == 0 ){
        xpts_loc = root_xpts;
        xmin_loc = (double) xmin_global ;
    }
    else{
        xpts_loc = non_root_xpts;
        xmin_loc = (double)xmin_global + xinc*(xpts_loc*rank + remainder);
    }

    std::vector<uint32_t> holder; // A local singular row  of calculated points that gets rewritten every increment in the y-direction
    std::vector<uint32_t> local_mesh; // local_mesh will contain the "entire" image calculated on a core

    //the max number of pixels for the output
    const int output_mesh_size = xpts_global * ypts_global;
    // The total number of pixels THIS core will calculate 
    const int mesh_size_loc = xpts_loc*ypts_global;

    /*
    Reserving the appropriate amount of space for each type of mesh
    */
    output_mesh.reserve(output_mesh_size);
    local_mesh.reserve(mesh_size_loc);
    holder.reserve(xpts_loc);

    double cur_y = -1.0;    //current y position
    double cur_x = xmin_loc;   //current x position 
    std::complex<double> c(cur_x,cur_y);    // complex number given by the current x and y

/*
CORE OF THE FUNCTION

loop through the mesh populating line by line in the x-direction

*/
    for(int iy = 0; iy < ypts_global; ++iy){

        for(int ix = 0; ix < xpts_loc; ++ix){

            //holder.insert(holder.begin() + ix,mbseries(c));
            holder.push_back(mbseries(c));
            c = std::complex<double>(c.real() + xinc, c.imag());
        }


        // append this row to the local image
        local_mesh.insert(local_mesh.end(), holder.begin(), holder.end());
        
        //populate output mesh with this row
        mesh_unifier(output_mesh, holder);

        //prepare for the next row of the image
        holder.clear();
        holder.reserve(xpts_loc);
        c = std::complex<double>(xmin_loc, c.imag()+yinc);
    }

    export_mesh(local_mesh);
}


void mesh_unifier(std::vector<uint32_t> &output_mesh,std::vector<uint32_t> &holder){

    int rank;   //  MPI rank of this core
    int worldSize;  //  MPI number of cores used
    int root;
    MPI_Comm_size(MPI_COMM_WORLD, &worldSize);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    root = 0; //worldSize - 1;

    int tag_size = 1; // confirm this data is the size of the vector
    int tag_row = 2; // confirm this data is the contents of the vector

    if(rank == root){
        

        // ROOT EXCLUSIVE: stores one row of x-points from rank 'n' core that will be appended to the final output mesh.
        std::vector<uint32_t> temp_mesh;

        output_mesh.insert(output_mesh.end(), holder.begin(), holder.end());

        for(int i_rank = 1; i_rank < worldSize; i_rank++){
            
            int row_length;
            MPI_Recv(&row_length, 1, MPI_INT, i_rank, i_rank,MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            temp_mesh.resize(row_length);

            MPI_Recv(temp_mesh.data(), row_length, MPI_UINT32_T, i_rank, i_rank, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            output_mesh.insert(output_mesh.end(), temp_mesh.begin(), temp_mesh.end());
            
            temp_mesh.clear();
        }       
    }
    else{

        int row_length = holder.size();
        MPI_Send(&row_length, 1, MPI_INT, root, rank, MPI_COMM_WORLD);
        MPI_Send(holder.data(), row_length,MPI_UINT32_T, root, rank, MPI_COMM_WORLD);
    }

}


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
    printf("Area of Mandelbrot set = %.8f\n", mb_area);

}


void export_mesh(std::vector<uint32_t> &mesh){

    int rank;
    int worldSize;

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &worldSize);

    char fname [100];

    std::sprintf(fname, "Cores-%d__res-%ld__Rank-%d__Mandelbrot_set.txt", worldSize, mesh.size(), rank );
    
    std::ofstream fout(fname);
    for(const auto &i : mesh){
        fout << i << "\n";
    }
    
    char move[150];
    std::sprintf(move, "mv %s outputs", fname);
    std::system(move);
}