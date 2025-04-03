#pragma once

#include "project/fem.h"

#ifdef __cplusplus
extern "C" {
#endif

#define _DISPLACEMENT (2.0 * 10e-7)
#define _E (211.0 * 10e9)
#define _NU (0.3)
#define _RHO (7.85 * 10e3)
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
