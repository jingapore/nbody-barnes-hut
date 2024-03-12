
#include "model.h"
#include <cmath>

#include <iostream>

using std::array;

array<double, DIM_SIZE> eval_force(const double * r1, double m1, const double * r2, double m2, double g, double rlimit){
    double M = m1 * m2;
    array<double, DIM_SIZE> f;
    
    double norm = 0;
    for(int c = 0; c < DIM_SIZE; c++){
        norm += pow(r1[c] - r2[c], 2);
    }

    norm = pow(norm, 0.5);
    double norm_cubed = std::max(pow(rlimit, 3), pow(norm, 3));
    
    for(int c = 0; c < DIM_SIZE; c++){
        f[c] = M * (r1[c] - r2[c]) * g;
        f[c] /= norm_cubed; // 3.0 over here is not dim_size, but is the universal eq to calculate force projection
    }
    return f;
}
