#include "../benchmark.h"
#include "fem.hpp"
#include "file_io.hpp"
#include <stdio.h>

void elasticity_solve(const char *meshfile, const char *outfile, double E,
                      double nu, double rho, double g) {
  // Read the mesh and the problem
  fem::Problem problem = fem::read_mesh_file(meshfile);
  // Boundary conditions are Dirichlet X and Y

  // Assemble and solve
  // femElasticityPrint(theProblem);

  // Write out the solution
}
