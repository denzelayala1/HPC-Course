#include <iostream>
#include <math.h>
#include <chrono>


double serial_Pi(int bins){
    double stepsize = 1.0/((double)bins);
    double x_input = stepsize/2;
    double fn_out = 0;

    for(int i = 0; i < bins; i++ ){

        fn_out += stepsize/(1+ x_input*x_input);
        x_input += stepsize;
    }

    return 4 * fn_out;

}

int main(int argc, char *argv[] ){

    auto start = std::chrono::high_resolution_clock::now();

    double pi_10_bins = serial_Pi(10);
    double pi_100_bins = serial_Pi(100);
    double pi_1000_bins = serial_Pi(1000);

    std::printf("serial pi_10_bins = %.8f\nserial pi_100_bins = %.8f\nserial pi_1000_bins = %.8f\n", pi_10_bins, pi_100_bins, pi_1000_bins);

    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);


    return 0;
}