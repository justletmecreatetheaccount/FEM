#pragma once

#include "project/fem.h"


#define _FACTOR (1e5)
#define _DISPLACEMENT (2.0 * 1e-6)
#define _E (170.0 * 1e9) // https://www.makeitfrom.com/material-properties/SAE-ASTM-Grade-G3500-F10007-Grey-Cast-Iron and https://en.wikipedia.org/wiki/Young%27s_modulus#Temperature_dependence and https://journals.sagepub.com/doi/full/10.1177/1687814018819563
#define _NU (0.29)
#define _RHO (7.5 * 1e3)
#define _G (9.81)

typedef struct Plot {
    femProblem* theProblem;
    femGeo* theGeometry;
    double* theForces;
    double* theSoluce;
} Plot;


Plot elasticity_solve(const char *meshfile, const char *outfile);
void elasticity_solve_csr(const char *meshfile, const char *outfile);




