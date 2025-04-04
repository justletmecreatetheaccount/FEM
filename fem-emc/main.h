#pragma once

#include "project/fem.h"

#ifdef __cplusplus
extern "C" {
#endif

#define _FACTOR (1e5)
#define _DISPLACEMENT (2.0 * 10e-7)
#define _E (170.0 * 10e9) // https://www.makeitfrom.com/material-properties/SAE-ASTM-Grade-G3500-F10007-Grey-Cast-Iron and https://en.wikipedia.org/wiki/Young%27s_modulus#Temperature_dependence
#define _NU (0.29)
#define _RHO (7.5 * 10e3)
#define _G (9.81)

typedef struct Plot {
    femProblem* theProblem;
    femGeo* theGeometry;
    double* theForces;
    double* theSoluce;
} Plot;


Plot elasticity_solve(const char *meshfile, const char *outfile);



#ifdef __cplusplus
}
#endif // C++
