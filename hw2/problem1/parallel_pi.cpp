#include <iostream>
#include <math.h>
#include "mpi.h"
#include <chrono>


void calc_Pi(int bins, int rank, int world_size){

    double t_init, t_final, t_elapsed;
    if(rank == 0 ){
    t_init = MPI_Wtime();
    }

    double step_size = 1.0/((double)bins );
    double x_input = step_size*( (double)rank+ 0.5);
    double local_area = 0;

    for(int i = rank; i < bins  ; i += world_size ){

            local_area += step_size/(1 + x_input*x_input);
            x_input += step_size* world_size ;
    }

    double tot_area;
    MPI_Reduce(&local_area, &tot_area, 1,MPI_DOUBLE, MPI_SUM, 0,MPI_COMM_WORLD);

    double pi = 4*tot_area;


    MPI_Barrier(MPI_COMM_WORLD);
    t_final = MPI_Wtime();
    t_elapsed = t_final - t_init;
    if(rank == 0 ){
    printf("CORES: %d   BINS: %d \nElapsed Time: %.10f seconds (%f milliseconds) \n", world_size, bins, t_elapsed, 1000* t_elapsed);
       std::printf("Estimated value of Pi = %.16f\n\n", pi);
    }

}

int main(int argc, char *argv[] ){



    int rank;
    int worldSize;
    int length;


    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &worldSize);

    if(argc > 1){

        for(int i=1; i < argc; i++){
            int bins = atoi(argv[i]);

            if(rank == 0 ){
            std::printf("You have selected %d bins\n\n", bins);
            }
            calc_Pi(bins, rank, worldSize);
        }
    }
    else{
        if(rank == 0 ){
        std::printf("Default run with 10,000 bins\n\n");
        }
        calc_Pi(10000, rank, worldSize);
    }

    MPI_Finalize();

    return 0;
}

