#include "fem.h"

# ifndef NOPOISSONCREATE

femPoissonProblem *femPoissonCreate(const char *filename)
{
    femGeo* theGeometry = geoMeshCreate(filename);
    femPoissonProblem *theProblem = malloc(sizeof(femPoissonProblem));
    theProblem->geo  = theGeometry;
    femMesh *theMesh = theGeometry->theElements;
    if (theMesh->nLocalNode == 4) {
        theProblem->space = femDiscreteCreate(4,FEM_QUAD);
        theProblem->rule = femIntegrationCreate(4,FEM_QUAD); }
    else if (theMesh->nLocalNode == 3) {
        theProblem->space = femDiscreteCreate(3,FEM_TRIANGLE);
        theProblem->rule = femIntegrationCreate(3,FEM_TRIANGLE); }
    theProblem->system = femFullSystemCreate(theMesh->nodes->nNodes);
    return theProblem;
}

# endif
# ifndef NOPOISSONBOUNDARY


void femPoissonFindBoundaryNodes(femPoissonProblem *theProblem)
{
    femGeo* theGeometry = theProblem->geo;  
    femMesh* theEdges = theGeometry->theEdges;
    printf("Number of edges %d \n", theEdges->nElem); 
    

    femDomain *theBoundary = malloc(sizeof(femDomain));
    theGeometry->nDomains++;
    theGeometry->theDomains = realloc(theGeometry->theDomains,theGeometry->nDomains*sizeof(femDomain*));
    theGeometry->theDomains[theGeometry->nDomains-1] = theBoundary;
    theBoundary->elem = malloc(theEdges->nElem * theEdges->nLocalNode*sizeof(int));

    int nBoundary = 0;
    int j = 0;

    for (int i = 0; i < theEdges->nElem * theEdges->nLocalNode; i++) {
        int alreadyIn = 0;
        for (int k = 0; k < j; k++) {
            if (theBoundary->elem[k] == theEdges->elem[i]) {
                alreadyIn = 1;
                break;
            }
        }
        
        if (alreadyIn == 0) {
            theBoundary->elem[j] = theEdges->elem[i];
            j++;
            nBoundary++;
        }
    }
    theBoundary->elem = realloc(theBoundary->elem,nBoundary*sizeof(int));
    theBoundary->nElem = nBoundary;
    
    
    
    theBoundary->mesh = NULL;
    sprintf(theBoundary->name,"Boundary");
}    
# endif
# ifndef NOPOISSONFREE

void femPoissonFree(femPoissonProblem *theProblem)
{

    geoMeshFree(theProblem->geo);
    femDiscreteFree(theProblem->space);
    femIntegrationFree(theProblem->rule);
    femFullSystemFree(theProblem->system);
    free(theProblem);
}
    
# endif
# ifndef NOPOISSONLOCAL

void femPoissonLocal(femPoissonProblem *theProblem, const int iElem, int *map, double *x, double *y)
{
    femMesh *theMesh = theProblem->geo->theElements;
    femNodes *Nodes = theMesh->nodes;
    for (int i = 0; i < theMesh->nLocalNode; i++){
        map[i] = theMesh->elem[iElem * theMesh->nLocalNode + i];
        x[i] = Nodes->X[map[i]];
        y[i] = Nodes->Y[map[i]];
    }

}

# endif
# ifndef NOPOISSONSOLVE

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

# endif



