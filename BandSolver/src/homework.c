

#include"fem.h"


#ifndef NORENUMBER 

double* assignableArrayPos;

int compareNodePos(const void* firstNode, const void* secondNode) {
    int* iOne = (int *) firstNode;
    int* iTwo = (int *) secondNode;
    double diff = assignableArrayPos[*iOne] - assignableArrayPos[*iTwo];
    return (diff < 0) - (diff > 0);
}

void femMeshRenumber(femMesh *theMesh, femRenumType renumType)
{
    int i;
    int* toso = (int*) malloc(sizeof(int) * theMesh->nodes->nNodes);
    for (i = 0; i < theMesh->nodes->nNodes; i++) {
        toso[i] = i;
    }
    
    switch (renumType) {
        case FEM_NO :
            for (i = 0; i < theMesh->nodes->nNodes; i++) 
                theMesh->nodes->number[i] = i;
            break;
// 
// A modifier :-)
// debut
//
        case FEM_XNUM : 
        assignableArrayPos = theMesh->nodes->X;
        qsort(toso, theMesh->nodes->nNodes, sizeof(int), compareNodePos);
        for (i = 0; i < theMesh->nodes->nNodes; i++) 
                theMesh->nodes->number[toso[i]] = i;
        break;
        case FEM_YNUM : 
        assignableArrayPos = theMesh->nodes->Y;
        qsort(toso, theMesh->nodes->nNodes, sizeof(int), compareNodePos);
        for (i = 0; i < theMesh->nodes->nNodes; i++) 
                theMesh->nodes->number[toso[i]] = i;
        break;            
// 
// end
//

        default : Error("Unexpected renumbering option"); 
        }
    free(toso);
}


#endif
#ifndef NOBAND 

int femMeshComputeBand(femMesh *theMesh)
{
    int bandwidth = 0;
    int nElem = theMesh->nElem;
    int nLocalNode = theMesh->nLocalNode;
    int *elem = theMesh->elem;

    // Loop over all elements
    for (int e = 0; e < nElem; e++) {
        // Loop over all nodes in the element
        for (int i = 0; i < nLocalNode; i++) {
            int node_i = theMesh->nodes->number[elem[e * nLocalNode + i]]; // Global index of node i
            for (int j = 0; j < nLocalNode; j++) {
                int node_j = theMesh->nodes->number[elem[e * nLocalNode + j]]; // Global index of node j
                int diff = abs(node_i - node_j); // Difference in global indices
                if (diff > bandwidth) {
                    bandwidth = diff; // Update bandwidth if this difference is larger
                }
            }
        }
    }

    return bandwidth + 1;
}


#endif
#ifndef NOBANDASSEMBLE


void femBandSystemAssemble(femBandSystem* myBandSystem, double *Aloc, double *Bloc, int *map, int nLoc)
{
    for (int i = 0; i < nLoc; i++) {
        int row = map[i];

        for (int j = 0; j < nLoc; j++) {
            int col = map[j];

            if (row <= col) {
                myBandSystem->A[row][col] += Aloc[i * nLoc + j];
            }
        }
        myBandSystem->B[row] += Bloc[i];
    }
}


#endif
#ifndef NOBANDELIMINATE


double  *femBandSystemEliminate(femBandSystem *myBand)
{
    double  **A, *B, factor;
    int     i, j, k, jend, size, band;
    A    = myBand->A;
    B    = myBand->B;
    size = myBand->size;
    band = myBand->band;


    /* Gauss elimination */

    for (k=0; k < size; k++) {
        if ( fabs(A[k][k]) <= 1e-16 ) {
            printf("Pivot index %d  ",k);
            printf("Pivot value %e  ",A[k][k]);
            Error("Cannot eliminate with such a pivot"); }
        jend = fmin(size,k + band);
        for (i = k + 1 ; i <  jend; i++) {
            factor = A[k][i] / A[k][k];
            for (j = i ; j < jend; j++) 
                A[i][j] = A[i][j] - A[k][j] * factor;
            B[i] = B[i] - B[k] * factor; 
        }
    } //code de fem.c devoir 3

    /* Back-substitution */

    for (i = size - 1; i >= 0; i--) {
        factor = 0;

        for (j = i + 1; j < fmin(i + band, size); j++)
            factor += A[i][j] * B[j];
        B[i] = (B[i] - factor) / A[i][i];
    }

    return(myBand->B);//code de fem.c devoir 3
}


#endif

