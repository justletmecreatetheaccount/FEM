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
    

    femDomain *theBoundary = malloc(sizeof(femDomain));
    theGeometry->nDomains++;
    theGeometry->theDomains = realloc(theGeometry->theDomains,theGeometry->nDomains*sizeof(femDomain*));
    theGeometry->theDomains[theGeometry->nDomains-1] = theBoundary;
    theBoundary->elem = malloc(theEdges->nElem * theEdges->nLocalNode*sizeof(int));

    int nBoundary = 0;
    int alreadyIn = 0;
    int j = 0;

    for (int i = 0; i < theEdges->nElem * theEdges->nLocalNode; i++) {
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
    
    // A completer :-)
}
    
# endif
# ifndef NOPOISSONLOCAL

void femPoissonLocal(femPoissonProblem *theProblem, const int iElem, int *map, double *x, double *y)
{
    femMesh *theMesh = theProblem->geo->theElements;
    femNodes *theNodes = theMesh->nodes;

    for (int i = 0; i < theMesh->nLocalNode; i++) {
        map[i] = theMesh->elem[iElem*theMesh->nLocalNode + i];
        x[i] = theNodes->X[map[i]];
        y[i] = theNodes->Y[map[i]];
    }

}

# endif
# ifndef NOPOISSONSOLVE

void femPoissonSolve(femPoissonProblem *theProblem)
{

    femMesh *theMesh = theProblem->geo->theElements;
    femDomain *theBoundary = geoGetDomain(theProblem->geo,"Boundary");
    femFullSystem *theSystem = theProblem->system;
    femIntegration *theRule = theProblem->rule;
    femDiscrete *theSpace = theProblem->space;
 
    if (theSpace->n > 4) Error("Unexpected discrete space size !");  
    double x[4],y[4],phi[4],dphidxsi[4],dphideta[4],dphidx[4],dphidy[4];
    int iElem,iInteg,iEdge,i,j,map[4];
    int nLocal = theMesh->nLocalNode;

    for (iElem = 0; iElem < theMesh->nElem; iElem++) {
        femPoissonLocal(theProblem,iElem,map,x,y);
        for (i = 0; i < theRule->n; i++) {
            femDiscretePhi2(theSpace,theRule->xsi[i],theRule->eta[i],phi);
            femDiscreteDphi2(theSpace,theRule->xsi[i],theRule->eta[i],dphidxsi,dphideta);
            double detJ = 0;
            for (j = 0; j < nLocal; j++) {
                detJ += dphidxsi[j]*x[j] + dphideta[j]*y[j]; }
            detJ = fabs(detJ);
            for (i = 0; i < nLocal; i++) {
                for (j = 0; j < nLocal; j++) {
                    double Aij = 0;
                    for (iInteg = 0; iInteg < theRule->n; iInteg++) {
                        femDiscreteDphi2(theSpace,theRule->xsi[iInteg],theRule->eta[iInteg],dphidx,dphidy);
                        Aij += (dphidx[i]*dphidx[j] + dphidy[i]*dphidy[j])*detJ*theRule->weight[iInteg]; }
                    theSystem->A[map[i]][map[j]] += Aij; }
                theSystem->B[map[i]] += 0; }
        }
    }
    for (i = 0; i < theBoundary->nElem; i++) {
        for (j = 0; j < nLocal; j++) {
            theSystem->B[theBoundary->elem[i]] = 0; }
    }
    femFullSystemConstrain(theSystem,0,0);
    femFullSystemEliminate(theSystem);
}

# endif



