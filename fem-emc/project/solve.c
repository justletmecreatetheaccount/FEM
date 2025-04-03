#include "../main.h"
#include "fem.h"
#include <stdio.h>

void elasticity_solve(const char *meshfile, const char *outfile)
{
    // Read the mesh and the problem
    femGeo *theGeometry = geoMeshRead(meshfile);
    femProblem *theProblem = femElasticityCreate(theGeometry, _E, _NU, _RHO, _G, PLANAR_STRAIN);

    // Boundary conditions
    femElasticityAddBoundaryCondition(theProblem, "Plate left", DIRICHLET_X, 
                                       _DISPLACEMENT);
    femElasticityAddBoundaryCondition(theProblem, "Plate top-left", DIRICHLET_X,
                                       _DISPLACEMENT);
    femElasticityAddBoundaryCondition(theProblem, "Plate right", DIRICHLET_X,
                                      -_DISPLACEMENT);
    femElasticityAddBoundaryCondition(theProblem, "Plate top-right", DIRICHLET_X,
                                      -_DISPLACEMENT);
    femElasticityAddBoundaryCondition(theProblem, "Plate bottom", DIRICHLET_Y,
                                      0.0);
    femElasticityAddBoundaryCondition(theProblem, "Plate bottom", DIRICHLET_X,
                                      0.0);

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
