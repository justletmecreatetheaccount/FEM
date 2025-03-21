#include "fem.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>


void femElasticityAssembleElements(femProblem *theProblem){
    femFullSystem  *theSystem = theProblem->system;
    femIntegration *theRule = theProblem->rule;
    femDiscrete    *theSpace = theProblem->space;
    femGeo         *theGeometry = theProblem->geometry;
    femNodes       *theNodes = theGeometry->theNodes;
    femMesh        *theMesh = theGeometry->theElements;
    femMesh        *theEdges = theGeometry->theEdges;
    double x[4],y[4],phi[4],dphidxsi[4],dphideta[4],dphidx[4],dphidy[4];
    int iElem,iInteg,iEdge,i,j,d,map[4],mapX[4],mapY[4];
    int nLocal = theMesh->nLocalNode;
    double *soluce = theProblem->soluce;
    double a   = theProblem->A;
    double b   = theProblem->B;
    double c   = theProblem->C;      
    double rho = theProblem->rho;
    double g   = theProblem->g;
    double **A = theSystem->A;
    double *B  = theSystem->B;
    
    
    for (iElem = 0; iElem < theMesh->nElem; iElem++) {
        for (j=0; j < nLocal; j++) {
            map[j]  = theMesh->elem[iElem*nLocal+j];
            mapX[j] = 2*map[j];
            mapY[j] = 2*map[j] + 1;
            x[j]    = theNodes->X[map[j]];
            y[j]    = theNodes->Y[map[j]];} 
        
        for (iInteg=0; iInteg < theRule->n; iInteg++) {    
            double xsi    = theRule->xsi[iInteg];
            double eta    = theRule->eta[iInteg];
            double weight = theRule->weight[iInteg];  
            femDiscretePhi2(theSpace,xsi,eta,phi);
            femDiscreteDphi2(theSpace,xsi,eta,dphidxsi,dphideta);
            
            double dxdxsi = 0.0;
            double dxdeta = 0.0;
            double dydxsi = 0.0; 
            double dydeta = 0.0;
            for (i = 0; i < theSpace->n; i++) {  
                dxdxsi += x[i]*dphidxsi[i];       
                dxdeta += x[i]*dphideta[i];   
                dydxsi += y[i]*dphidxsi[i];   
                dydeta += y[i]*dphideta[i]; }
            double jac = fabs(dxdxsi * dydeta - dxdeta * dydxsi);
            
            for (i = 0; i < theSpace->n; i++) {    
                dphidx[i] = (dphidxsi[i] * dydeta - dphideta[i] * dydxsi) / jac;       
                dphidy[i] = (dphideta[i] * dxdxsi - dphidxsi[i] * dxdeta) / jac; }            
            for (i = 0; i < theSpace->n; i++) { 
                for(j = 0; j < theSpace->n; j++) {
                    A[mapX[i]][mapX[j]] += (dphidx[i] * a * dphidx[j] + 
                                            dphidy[i] * c * dphidy[j]) * jac * weight;                                                                                            
                    A[mapX[i]][mapY[j]] += (dphidx[i] * b * dphidy[j] + 
                                            dphidy[i] * c * dphidx[j]) * jac * weight;                                                                                           
                    A[mapY[i]][mapX[j]] += (dphidy[i] * b * dphidx[j] + 
                                            dphidx[i] * c * dphidy[j]) * jac * weight;                                                                                            
                    A[mapY[i]][mapY[j]] += (dphidy[i] * a * dphidy[j] + 
                                            dphidx[i] * c * dphidx[j]) * jac * weight; }}
             for (i = 0; i < theSpace->n; i++) {
                B[mapY[i]] -= phi[i] * g * rho * jac * weight; }}} 
}

void femElasticityAssembleNeumann(femProblem *theProblem){
    femFullSystem  *theSystem = theProblem->system;
    femIntegration *theRule = theProblem->ruleEdge;
    femDiscrete    *theSpace = theProblem->spaceEdge;
    femGeo         *theGeometry = theProblem->geometry;
    femNodes       *theNodes = theGeometry->theNodes;
    femMesh        *theEdges = theGeometry->theEdges;
    double x[2],y[2],phi[2], dphidxsi[2];
    int iBnd,iElem,iInteg,iEdge,i,j,d,map[2],mapU[2];
    int nLocal = 2;
    double *B  = theSystem->B;


    for(iBnd=0; iBnd < theProblem->nBoundaryConditions; iBnd++){
        femBoundaryCondition *theCondition = theProblem->conditions[iBnd];
        femDomain *theDomain = theCondition->domain;
        femBoundaryType type = theCondition->type;
        double value = theCondition->value;

        if (type == NEUMANN_X) {
            for (iElem = 0; iElem < theDomain->nElem; iElem++) {
                for (j=0; j < nLocal; j++) {
                    map[j] = theDomain->elem[iElem*nLocal+j];
                    mapU[j] = 2*map[j];
                    x[j] = theNodes->X[map[j]];
                    y[j] = theNodes->Y[map[j]];}
                for (iInteg=0; iInteg < theRule->n; iInteg++) {    
                    double xsi = theRule->xsi[iInteg];
                    double weight = theRule->weight[iInteg];  
                    femDiscretePhi(theSpace,xsi,phi);
                    femDiscreteDphi(theSpace, xsi, dphidxsi);

                    // Compute dx/dξ and dy/dξ
                    double dxdxsi = 0.0;
                    double dydxsi = 0.0;
                    for (int k = 0; k < nLocal; k++) {
                        dxdxsi += x[k] * dphidxsi[k];  // dx/dξ
                        dydxsi += y[k] * dphidxsi[k];  // dy/dξ
                    }

                    // Compute Jacobian determinant
                    double jac = sqrt(dxdxsi * dxdxsi + dydxsi * dydxsi);

                    for (i = 0; i < theSpace->n; i++) {
                        B[mapU[i]] += phi[i] * value * weight * jac; }}}

        } else if (type == NEUMANN_Y) {
            for (iElem = 0; iElem < theDomain->nElem; iElem++) {
                printf("iElem = %d\n", iElem);
                for (j=0; j < nLocal; j++) {
                    map[j] = theDomain->elem[iElem*nLocal+j];
                    printf("map[j] = %d\n", map[j]);
                    mapU[j] = 2*map[j] + 1;
                    x[j] = theNodes->X[map[j]];
                    y[j] = theNodes->Y[map[j]];}

                for (iInteg=0; iInteg < theRule->n; iInteg++) {
                    double xsi  = theRule->xsi[iInteg];
                    double weight = theRule->weight[iInteg];  
                    femDiscretePhi(theSpace,xsi,phi);
                    femDiscreteDphi(theSpace, xsi, dphidxsi);

                    // Compute dx/dξ and dy/dξ
                    double dxdxsi = 0.0;
                    double dydxsi = 0.0;
                    for (int k = 0; k < nLocal; k++) {
                        dxdxsi += x[k] * dphidxsi[k];  // dx/dξ
                        dydxsi += y[k] * dphidxsi[k];  // dy/dξ
                    }

                    // Compute Jacobian determinant
                    double jac = sqrt(dxdxsi * dxdxsi + dydxsi * dydxsi);
                    printf("jac = %f\n", jac);

                    for (i = 0; i < nLocal; i++) {
                        printf("B[mapU[%i]] = %f\n", mapU[i], B[mapU[i]]);
                        B[mapU[i]] += phi[i] * value * jac * weight; }}}
        }
    }
}

double** Aint;
double* Bint;

double *femElasticitySolve(femProblem *theProblem){
    femFullSystem *theSystem = theProblem->system;

    for (int i = 0; i < theSystem->size; i++) {
        theSystem->B[i] = 0;
    }

    femElasticityAssembleNeumann(theProblem);
    femElasticityAssembleElements(theProblem);

    Aint = malloc(theSystem->size * sizeof(double*));
    Bint = malloc(theSystem->size * sizeof(double));
    for (int i = 0; i < theSystem->size; i++) {
        Aint[i] = malloc(theSystem->size*sizeof(double));
        Bint[i] = theSystem->B[i];
        for (int j = 0; j < theSystem->size; j++) {
            Aint[i][j] = theSystem->A[i][j];
        }
    }

    int* theConstrainedNodes = theProblem->constrainedNodes;
    for (int i = 0; i < theSystem->size; i++) {
        if (theConstrainedNodes[i] != -1) {
            double value = theProblem->conditions[theConstrainedNodes[i]]->value;
            femFullSystemConstrain(theSystem, i, value);
        }
    }

    femFullSystemEliminate(theSystem);

    for (int i = 0; i < theSystem->size; i++) {
        theProblem->soluce[i] = theSystem->B[i];
    }

    return theProblem->soluce;
}

double * femElasticityForces(femProblem *theProblem){        
    femFullSystem *theSystem = theProblem->system;
    double *soluce = theProblem->soluce;
    double *residuals = theProblem->residuals;
    int size = theSystem->size;
    double **A = Aint;
    double *B = Bint;
    int i,j;
    for (i = 0; i < size; i++) {
        residuals[i] = -B[i];
        for (j = 0; j < size; j++) residuals[i] += A[i][j] * soluce[j]; }
    for (i = 0; i < size; i++) {
        free(A[i]);
    }
    free(A);
    free(B);
    return residuals;
}

//  ==== Minimum displacement          :  0.0000000e+00 [m] 
//  ==== Maximum displacement          :  4.5464210e-06 [m] 
//  ==== Global horizontal force       :  2.6302062e-10 [N] 
//  ==== Global vertical force         : -7.0221752e+04 [N] 
//  ==== Weight                        :  5.0221752e+04 [N] 