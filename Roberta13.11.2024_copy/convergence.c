#include "comum.h"
//---check convergence---

void convergence(int itc, double error, double residual_p, double residual_u, double residual_v){
    itc++;
    error = fmax(residual_p, residual_u);
    error = fmax(residual_v, error);
}   