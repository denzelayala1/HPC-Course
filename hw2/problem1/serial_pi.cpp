#include <iostream>
#include <math.h>


double serial_Pi(int bins){
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

    double pi_10_bins = serial_Pi(10);
    double pi_100_bins = serial_Pi(100);
    double pi_1000_bins = serial_Pi(1000);

    std::printf("serial pi_10_bins = %.8f\nserial pi_100_bins = %.8f\nserial pi_1000_bins = %.8f\n ", pi_10_bins, pi_100_bins, pi_1000_bins);

    return 0;
}