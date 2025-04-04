#include "main.h"
#include "fem.h"
#include <stdio.h>


void elasticity_solve_csr(const char *meshfile, const char *outfile, double E,
                      double nu, double rho, double g) {
  // Read the mesh and the problem
  femGeo *theGeometry = geoMeshRead(meshfile);

  femProblemCsr *theProblem =
      femElasticityCreateCsr(theGeometry, E, nu, rho, g, PLANAR_STRAIN);

  if (!strcmp(meshfile, "../data/mesh_ref.txt")) {

    // Boundary conditions are Dirichlet X and Y
    femElasticityAddBoundaryConditionCsr(theProblem, "Base", DIRICHLET_X, 0.0);
    femElasticityAddBoundaryConditionCsr(theProblem, "Base", DIRICHLET_Y, 0.0);
    femElasticityAddBoundaryConditionCsr(theProblem, "Symmetry", DIRICHLET_X, 0.0);
    femElasticityAddBoundaryConditionCsr(theProblem, "Symmetry", DIRICHLET_Y, 0.0);
  } else {

    // Boundary conditions are Dirichlet X and Y
    femElasticityAddBoundaryConditionCsr(theProblem, "Plate left", DIRICHLET_X, 1e-6);
    femElasticityAddBoundaryConditionCsr(theProblem, "Plate right", DIRICHLET_X,
                                      -1e-6);
    femElasticityAddBoundaryConditionCsr(theProblem, "Plate bottom", DIRICHLET_Y,
                                      0.0);
    femElasticityAddBoundaryConditionCsr(theProblem, "Plate bottom", DIRICHLET_X, 0.0);
  }

  // Assemble and solve

  double *theSoluce = femElasticitySolveCsr(theProblem);
  //double* theForces = femElasticityForces(theProblem);
  double *normDisplacement = malloc(theGeometry->theNodes->nNodes * sizeof(double));
  //double *forcesX = malloc(theGeometry->theNodes->nNodes * sizeof(double));
  //double *forcesY = malloc(theGeometry->theNodes->nNodes * sizeof(double));
  
  // for (int i=0; i<theGeometry->theNodes->nNodes; i++){
  //     normDisplacement[i] = sqrt(theSoluce[2*i+0]*theSoluce[2*i+0] + 
  //                                 theSoluce[2*i+1]*theSoluce[2*i+1]);
  //     forcesX[i] = theForces[2*i+0];
  //     forcesY[i] = theForces[2*i+1]; }

  double hMin = femMin(normDisplacement,theGeometry->theNodes->nNodes);  
  double hMax = femMax(normDisplacement,theGeometry->theNodes->nNodes);  
  printf(" ==== Minimum displacement          : %14.7e [m] \n",hMin);
  printf(" ==== Maximum displacement          : %14.7e [m] \n",hMax);

//
//  -5- Calcul de la force globaleresultante
//

  double theGlobalForce[2] = {0, 0};
  // for (int i=0; i<theProblem->geometry->theNodes->nNodes; i++) {
  //     theGlobalForce[0] += theForces[2*i+0];
  //     theGlobalForce[1] += theForces[2*i+1]; }
  printf(" ==== Global horizontal force       : %14.7e [N] \n",theGlobalForce[0]);
  printf(" ==== Global vertical force         : %14.7e [N] \n",theGlobalForce[1]);

  // Write out the solution
  int nNodes = theGeometry->theNodes->nNodes;
  femSolutionWrite(nNodes, 2, theSoluce, outfile);

  // free the allocated ressources
  femElasticityFreeCsr(theProblem);
  geoFinalize(theGeometry);
}

Plot elasticity_solve(const char *meshfile, const char *outfile) {
  // Read the mesh and the problem
  femGeo *theGeometry = geoMeshRead(meshfile);
  femProblem *theProblem =
      femElasticityCreate(theGeometry, _E, _NU, _RHO, _G, PLANAR_STRAIN);

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
