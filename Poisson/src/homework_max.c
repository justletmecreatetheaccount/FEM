#include "fem.h"
#include <stdio.h>

#ifndef NOPOISSONCREATE

femPoissonProblem *femPoissonCreate(const char *filename) {
  femGeo *theGeometry = geoMeshCreate(filename);
  femPoissonProblem *theProblem = malloc(sizeof(femPoissonProblem));
  theProblem->geo = theGeometry;
  femMesh *theMesh = theGeometry->theElements;
  if (theMesh->nLocalNode == 4) {
    theProblem->space = femDiscreteCreate(4, FEM_QUAD);
    theProblem->rule = femIntegrationCreate(4, FEM_QUAD);
  } else if (theMesh->nLocalNode == 3) {
    theProblem->space = femDiscreteCreate(3, FEM_TRIANGLE);
    theProblem->rule = femIntegrationCreate(3, FEM_TRIANGLE);
  }
  theProblem->system = femFullSystemCreate(theMesh->nodes->nNodes);
  return theProblem;
}

#endif
#ifndef NOPOISSONBOUNDARY

void femPoissonFindBoundaryNodes(femPoissonProblem *theProblem) {
  femGeo *theGeometry = theProblem->geo;
  // Mesh of boundary edges
  femMesh *theEdges = theGeometry->theEdges;
  int nBoundary = 0;

  int *count = malloc(theEdges->nodes->nNodes * sizeof(int));
  for (int i = 0; i < theEdges->nodes->nNodes; i++) {
    count[i] = 0;
  }
  for (int i = 0; i < theEdges->nElem; i++) {
    for (int j = 0; j < theEdges->nLocalNode; j++) {
      count[theEdges->elem[i * theEdges->nLocalNode + j]]++;
    }
  }
  for (int i = 0; i < theEdges->nodes->nNodes; i++) {
    if (count[i] > 0) {
      nBoundary++;
    }
  }

  femDomain *theBoundary = malloc(sizeof(femDomain));
  theGeometry->nDomains++;
  theGeometry->theDomains = realloc(
      theGeometry->theDomains, theGeometry->nDomains * sizeof(femDomain *));
  theGeometry->theDomains[theGeometry->nDomains - 1] = theBoundary;
  theBoundary->nElem = nBoundary;
  theBoundary->elem = malloc(nBoundary * sizeof(int));
  theBoundary->mesh = NULL;
  sprintf(theBoundary->name, "Boundary");

  int k = 0;
  for (int i = 0; i < theEdges->nodes->nNodes; i++) {
    if (count[i] > 0) {
      theBoundary->elem[k] = i;
      k++;
    }
  }

  free(count);
}

#endif
#ifndef NOPOISSONFREE

void femPoissonFree(femPoissonProblem *theProblem) {
  geoMeshFree(theProblem->geo);
  femDiscreteFree(theProblem->space);
  femIntegrationFree(theProblem->rule);
  femFullSystemFree(theProblem->system);
  free(theProblem);

}

#endif
#ifndef NOPOISSONLOCAL

void femPoissonLocal(femPoissonProblem *theProblem, const int iElem, int *map,
                     double *x, double *y) {
  femMesh *theMesh = theProblem->geo->theElements;

  femNodes *theNodes = theMesh->nodes;
  int nLocalNode = theMesh->nLocalNode;

  for (int j = 0; j < nLocalNode; j++) {
    int nodeIndex = theMesh->elem[iElem * nLocalNode + j];
    map[j] = nodeIndex;
    x[j] = theNodes->X[nodeIndex];
    y[j] = theNodes->Y[nodeIndex];
  }
}

#endif
#ifndef NOPOISSONSOLVE

void femPoissonSolve(femPoissonProblem *theProblem) {

  femMesh *theMesh = theProblem->geo->theElements;
  femDomain *theBoundary = geoGetDomain(theProblem->geo, "Boundary");
  femFullSystem *theSystem = theProblem->system;
  femIntegration *theRule = theProblem->rule;
  femDiscrete *theSpace = theProblem->space;

  if (theSpace->n > 4)
    Error("Unexpected discrete space size !");
  double x[4], y[4], phi[4], dphidxsi[4], dphideta[4], dphidx[4], dphidy[4];
  int iElem, iInteg, iEdge, i, j, map[4];
  int nLocal = theMesh->nLocalNode;

  double xsi, eta, weight;
  double dxdxsi, dxdeta, dydxsi, dydeta;

  // On parcourt nos elements
  for (int e = 0; e < theMesh->nElem; e++) {
    femPoissonLocal(theProblem, e, map, x, y);

    // On parcourt les points d'integration
    for (int k = 0; k < theRule->n; k++) {
      eta = theRule->eta[k];
      xsi = theRule->xsi[k];
      weight = theRule->weight[k];

      // On calcule la valeur des fonctions de forme au
      // point d'integration et leur dérivées.
      femDiscretePhi2(theSpace, xsi, eta, phi);
      femDiscreteDphi2(theSpace, xsi, eta, dphidxsi, dphideta);

      dxdxsi = 0; dxdeta = 0; dydxsi = 0; dydeta = 0;

      // On calcule les dérivées des coordonnées x et y
      // par rapport à xsi et eta au point d'integration.
      // Ceci ce fait en itérant sur les points locaux de
      // l'élément et en regardant la valeur de la dérivée
      // de la fonction de forme aux eta et xsi actuels.
      for (int i = 0; i < theSpace->n; i++) {
        dxdxsi += x[i] * dphidxsi[i];
        dxdeta += x[i] * dphideta[i];
        dydxsi += y[i] * dphidxsi[i];
        dydeta += y[i] * dphideta[i];
      }

      double J_e = fabs(dxdxsi * dydeta - dxdeta * dydxsi);

      // On calcule les dérivées de chacune des fonctions de forme
      // sur l'élément !!!sans diviser par la jacobienne!!! pour
      // le moment pour opti un peu.
      for (int i = 0; i < theSpace->n; i++) {
          dphidx[i] = (dphidxsi[i] * dydeta - dphideta[i] * dydxsi);
          dphidy[i] = (dphideta[i] * dxdxsi - dphidxsi[i] * dxdeta);
      }

      // Remplissage de la matrice A et du vecteur B
      // en passant par map pour avoir les indices globaux
      // à partir des indices locaux. On divise par la
      // jacobienne.
      for (int i = 0; i < theSpace->n; i++) {
        for (int j = 0; j < theSpace->n; j++) {
          theSystem->A[map[i]][map[j]] += (dphidx[i] * dphidx[j] + dphidy[i] * dphidy[j]) * weight / J_e;
        }
        theSystem->B[map[i]] += phi[i] * weight * J_e;
      }

    }
  }

  for (i = 0; i < theBoundary->nElem; i++) {
    femFullSystemConstrain(theSystem, theBoundary->elem[i], 0.0);
  }

  femFullSystemEliminate(theSystem);
}

#endif
