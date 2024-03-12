
#ifndef MODEL_H_N_BODY
#define MODEL_H_N_BODY

#include <array>
#include "../misc/body.h"

using std::array;

array<double, DIM_SIZE> eval_force(const double * r1, double m1, const double * r2, double m2, double g, double rlimit);

#endif // MODEL_H_N_BODY
