#include "../benchmark.h"
#include <stdio.h>
#include "fem.h"

void elasticity_solve(const char *meshfile, const char *outfile, double E, double nu, double rho, double g) {
  // Read the mesh and the problem
  femGeo *theGeometry = geoMeshRead(meshfile);
  femProblem *theProblem = femElasticityCreate(theGeometry, E, nu, rho, g, PLANAR_STRAIN);

  // Boundary conditions are Dirichlet X and Y
  femElasticityAddBoundaryCondition(theProblem, "Plate left", DIRICHLET_X, 21e3);
  femElasticityAddBoundaryCondition(theProblem, "Plate right", DIRICHLET_X, 21e3);    

  // Assemble and solve
  // femElasticityPrint(theProblem);
  double *theSoluce = femElasticitySolve(theProblem);
  femElasticityForces(theProblem);
  
  // Write out the solution
  int nNodes = theGeometry->theNodes->nNodes;
  femSolutionWrite(nNodes, 2, theSoluce, outfile);

  // free the allocated ressources
  femElasticityFree(theProblem);
  geoFinalize(theGeometry);
}
