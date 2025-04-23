#include "../benchmark.h"
#include "cuda_utils.cuh"
#include "fem.hpp"
#include <bits/stdc++.h>
#include <cstdlib>
#include <cstring>

void elasticity_solve(const char *meshfile, const char *outfile, double E,
                      double nu, double rho, double g) {
  // Read the mesh and create the problem
  fem::Problem problem = fem::Problem(meshfile, E, nu, rho, g);
  // Boundary conditions are Dirichlet X and Y

  // Assemble and solve
  // femElasticityPrint(theProblem);
  //

  cuda::test();
  // Write out the solution
}
