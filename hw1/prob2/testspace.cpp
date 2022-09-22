#include <iostream>
#include <iomanip>
#include <complex>
#include <vector>
using namespace std;

int main(){
    complex<int> comps(0, 0);
    int inc = 1;
    cout << comps << "\t" ;


    for(int y = 0; y < 3 ; ++y){
        
        for(int ix = 0; ix < 3; ++ix){

                comps = complex<int>(comps.real()+inc, comps.imag());
                cout << comps << "\t" ;
        }
        comps = complex<int>(0, comps.imag() + inc);
    }
    return 0;
}