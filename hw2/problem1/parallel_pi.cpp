#include <iostream>
#include <math.h>
#include "mpi.h"


double calc_Pi(int bins, int ranks){

    double stepsize = 1.0/((double)bins -1 );
    double x_input = stepsize/2;
    double fn_out = 0;

    for(int i = 0; i < bins - 1; i++ ){

        fn_out += stepsize/(1+ x_input*x_input);
        x_input += stepsize;
    }

    return 4 * fn_out;

}

int main(int argc, char *argv[] ){

    int rank;
    int worldSize;
    char name[80];
    int length;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &worldSize);

    std::printf("the world size = %d", worldSize);



    return 0;
}