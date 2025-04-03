#include "main.h"
#include "fem.h"
#include <stdio.h>

Plot elasticity_solve(const char *meshfile, const char *outfile)
{
    // Read the mesh and the problem
    femGeo *theGeometry = geoMeshRead(meshfile);
    femProblem *theProblem = femElasticityCreate(theGeometry, _E, _NU, _RHO, _G, PLANAR_STRAIN);

  if (!strcmp(meshfile, "../data/mesh_ref.txt")) {

    // Boundary conditions are Dirichlet X and Y
    femElasticityAddBoundaryCondition(theProblem, "Base", DIRICHLET_X, 0.0);
    femElasticityAddBoundaryCondition(theProblem, "Base", DIRICHLET_Y, 0.0);
    femElasticityAddBoundaryCondition(theProblem, "Symmetry", DIRICHLET_X, 0.0);
    femElasticityAddBoundaryCondition(theProblem, "Symmetry", DIRICHLET_Y, 0.0);
  } else {

    // Boundary conditions are Dirichlet X and Y
    femElasticityAddBoundaryCondition(theProblem, "Plate left", DIRICHLET_X, _DISPLACEMENT);
    femElasticityAddBoundaryCondition(theProblem, "Plate right", DIRICHLET_X,
                                      -_DISPLACEMENT);
    femElasticityAddBoundaryCondition(theProblem, "Plate bottom", DIRICHLET_Y,
                                      0.0);
    femElasticityAddBoundaryCondition(theProblem, "Plate bottom", DIRICHLET_X, 0.0);
  }

    // Assemble and solve
    // femElasticityPrint(theProblem);
    double *theSoluce = femElasticitySolve(theProblem);
    double *theForces = femElasticityForces(theProblem);
    Plot res = {theProblem, theGeometry, theForces, theSoluce};

    int nNodes = theGeometry->theNodes->nNodes;
    femSolutionWrite(nNodes, 2, theSoluce, outfile);
    return res;

}
