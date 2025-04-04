/*
 *  fem.c
 *  Library for LEPL1110 : Finite Elements for dummies
 *
 *  Copyright (C) 2021 UCL-IMMC : Vincent Legat
 *  All rights reserved.
 *
 */

#include "fem.h"
#include <bits/pthreadtypes.h>
#include <stdio.h>
#include <stdlib.h>

femGeo theGeometry;

femGeo *geoGetGeometry() { return &theGeometry; }


void geoFinalize() {

  if (theGeometry.theNodes) {
    free(theGeometry.theNodes->X);
    free(theGeometry.theNodes->Y);
    free(theGeometry.theNodes);
  }
  if (theGeometry.theElements) {
    free(theGeometry.theElements->elem);
    free(theGeometry.theElements);
  }
  if (theGeometry.theEdges) {
    free(theGeometry.theEdges->elem);
    free(theGeometry.theEdges);
  }
  for (int i = 0; i < theGeometry.nDomains; i++) {
    free(theGeometry.theDomains[i]->elem);
    free(theGeometry.theDomains[i]);
  }
  free(theGeometry.theDomains);
}

femGeo *geoMeshRead(const char *filename) {
  FILE *file = fopen(filename, "r");

  int trash, *elem;

  femNodes *theNodes = malloc(sizeof(femNodes));
  theGeometry.theNodes = theNodes;
  ErrorScan(fscanf(file, "Number of nodes %d \n", &theNodes->nNodes));
  theNodes->X = malloc(sizeof(double) * (theNodes->nNodes));
  theNodes->Y = malloc(sizeof(double) * (theNodes->nNodes));
  for (int i = 0; i < theNodes->nNodes; i++) {
    ErrorScan(fscanf(file, "%d : %le %le \n", &trash, &theNodes->X[i],
                     &theNodes->Y[i]));
  }

  femMesh *theEdges = malloc(sizeof(femMesh));
  theGeometry.theEdges = theEdges;
  theEdges->nLocalNode = 2;
  theEdges->nodes = theNodes;
  ErrorScan(fscanf(file, "Number of edges %d \n", &theEdges->nElem));
  theEdges->elem = malloc(sizeof(int) * theEdges->nLocalNode * theEdges->nElem);
  for (int i = 0; i < theEdges->nElem; ++i) {
    elem = theEdges->elem;
    ErrorScan(fscanf(file, "%6d : %6d %6d \n", &trash, &elem[2 * i],
                     &elem[2 * i + 1]));
  }

  femMesh *theElements = malloc(sizeof(femMesh));
  theGeometry.theElements = theElements;
  theElements->nLocalNode = 0;
  theElements->nodes = theNodes;
  char elementType[MAXNAME];
  ErrorScan(
      fscanf(file, "Number of %s %d \n", elementType, &theElements->nElem));
  if (strncasecmp(elementType, "triangles", MAXNAME) == 0) {
    theElements->nLocalNode = 3;
    theElements->elem =
        malloc(sizeof(int) * theElements->nLocalNode * theElements->nElem);
    for (int i = 0; i < theElements->nElem; ++i) {
      elem = theElements->elem;
      ErrorScan(fscanf(file, "%6d : %6d %6d %6d \n", &trash, &elem[3 * i],
                       &elem[3 * i + 1], &elem[3 * i + 2]));
    }
  }
  if (strncasecmp(elementType, "quads", MAXNAME) == 0) {
    theElements->nLocalNode = 4;
    theElements->elem =
        malloc(sizeof(int) * theElements->nLocalNode * theElements->nElem);
    for (int i = 0; i < theElements->nElem; ++i) {
      elem = theElements->elem;
      ErrorScan(fscanf(file, "%6d : %6d %6d %6d %6d \n", &trash, &elem[4 * i],
                       &elem[4 * i + 1], &elem[4 * i + 2], &elem[4 * i + 3]));
    }
  }

  ErrorScan(fscanf(file, "Number of domains %d\n", &theGeometry.nDomains));
  int nDomains = theGeometry.nDomains;
  theGeometry.theDomains = malloc(sizeof(femDomain *) * nDomains);
  for (int iDomain = 0; iDomain < nDomains; iDomain++) {
    femDomain *theDomain = malloc(sizeof(femDomain));
    theGeometry.theDomains[iDomain] = theDomain;
    theDomain->mesh = theEdges;
    ErrorScan(fscanf(file, "  Domain : %6d \n", &trash));
    ErrorScan(fscanf(file, "  Name : %[^\n]s \n", (char *)&theDomain->name));
    ErrorScan(fscanf(file, "  Number of elements : %6d\n", &theDomain->nElem));
    theDomain->elem = malloc(sizeof(int) * 2 * theDomain->nElem);
    for (int i = 0; i < theDomain->nElem; i++) {
      ErrorScan(fscanf(file, "%6d", &theDomain->elem[i]));
      if ((i + 1) != theDomain->nElem && (i + 1) % 10 == 0)
        ErrorScan(fscanf(file, "\n"));
    }
  }

  fclose(file);
  return &theGeometry;
}

void femSolutionWrite(int nNodes, int nfields, double *data,
                      const char *filename) {
  FILE *file = fopen(filename, "w");
  if (!file) {
    printf("Error at %s:%d\nUnable to open file %s\n", __FILE__, __LINE__,
           filename);
    exit(-1);
  }
  fprintf(file, "Size %d,%d\n", nNodes, nfields);
  for (int i = 0; i < nNodes; i++) {
    for (int j = 0; j < nfields - 1; j++) {
      fprintf(file, "%.18le,", data[i * nfields + j]);
    }
    fprintf(file, "%.18le", data[i * nfields + nfields - 1]);
    fprintf(file, "\n");
  }
  fclose(file);
}


int geoGetDomain(char *name) {
  int theIndex = -1;
  int nDomains = theGeometry.nDomains;
  for (int iDomain = 0; iDomain < nDomains; iDomain++) {
    femDomain *theDomain = theGeometry.theDomains[iDomain];
    if (strncasecmp(name, theDomain->name, MAXNAME) == 0)
      theIndex = iDomain;
  }
  return theIndex;
}

static const double _gaussQuad4Xsi[4] = {-0.577350269189626, -0.577350269189626,
                                         0.577350269189626, 0.577350269189626};
static const double _gaussQuad4Eta[4] = {0.577350269189626, -0.577350269189626,
                                         -0.577350269189626, 0.577350269189626};
static const double _gaussQuad4Weight[4] = {
    1.000000000000000, 1.000000000000000, 1.000000000000000, 1.000000000000000};
static const double _gaussTri3Xsi[3] = {0.166666666666667, 0.666666666666667,
                                        0.166666666666667};
static const double _gaussTri3Eta[3] = {0.166666666666667, 0.166666666666667,
                                        0.666666666666667};
static const double _gaussTri3Weight[3] = {0.166666666666667, 0.166666666666667,
                                           0.166666666666667};
static const double _gaussEdge2Xsi[2] = {0.577350269189626, -0.577350269189626};
static const double _gaussEdge2Weight[2] = {1.000000000000000,
                                            1.000000000000000};

femIntegration *femIntegrationCreate(int n, femElementType type) {
  femIntegration *theRule = malloc(sizeof(femIntegration));
  if (type == FEM_QUAD && n == 4) {
    theRule->n = 4;
    theRule->xsi = _gaussQuad4Xsi;
    theRule->eta = _gaussQuad4Eta;
    theRule->weight = _gaussQuad4Weight;
  } else if (type == FEM_TRIANGLE && n == 3) {
    theRule->n = 3;
    theRule->xsi = _gaussTri3Xsi;
    theRule->eta = _gaussTri3Eta;
    theRule->weight = _gaussTri3Weight;
  } else if (type == FEM_EDGE && n == 2) {
    theRule->n = 2;
    theRule->xsi = _gaussEdge2Xsi;
    theRule->eta = NULL;
    theRule->weight = _gaussEdge2Weight;
  } else
    Error("Cannot create such an integration rule !");
  return theRule;
}

void femIntegrationFree(femIntegration *theRule) { free(theRule); }

void _q1c0_x(double *xsi, double *eta) {
  xsi[0] = 1.0;
  eta[0] = 1.0;
  xsi[1] = -1.0;
  eta[1] = 1.0;
  xsi[2] = -1.0;
  eta[2] = -1.0;
  xsi[3] = 1.0;
  eta[3] = -1.0;
}

void _q1c0_phi(double xsi, double eta, double *phi) {
  phi[0] = (1.0 + xsi) * (1.0 + eta) / 4.0;
  phi[1] = (1.0 - xsi) * (1.0 + eta) / 4.0;
  phi[2] = (1.0 - xsi) * (1.0 - eta) / 4.0;
  phi[3] = (1.0 + xsi) * (1.0 - eta) / 4.0;
}

void _q1c0_dphidx(double xsi, double eta, double *dphidxsi, double *dphideta) {
  dphidxsi[0] = (1.0 + eta) / 4.0;
  dphidxsi[1] = -(1.0 + eta) / 4.0;
  dphidxsi[2] = -(1.0 - eta) / 4.0;
  dphidxsi[3] = (1.0 - eta) / 4.0;
  dphideta[0] = (1.0 + xsi) / 4.0;
  dphideta[1] = (1.0 - xsi) / 4.0;
  dphideta[2] = -(1.0 - xsi) / 4.0;
  dphideta[3] = -(1.0 + xsi) / 4.0;
}

void _p1c0_x(double *xsi, double *eta) {
  xsi[0] = 0.0;
  eta[0] = 0.0;
  xsi[1] = 1.0;
  eta[1] = 0.0;
  xsi[2] = 0.0;
  eta[2] = 1.0;
}

void _p1c0_phi(double xsi, double eta, double *phi) {
  phi[0] = 1 - xsi - eta;
  phi[1] = xsi;
  phi[2] = eta;
}

void _p1c0_dphidx(double xsi, double eta, double *dphidxsi, double *dphideta) {
  dphidxsi[0] = -1.0;
  dphidxsi[1] = 1.0;
  dphidxsi[2] = 0.0;
  dphideta[0] = -1.0;
  dphideta[1] = 0.0;
  dphideta[2] = 1.0;
}

void _e1c0_x(double *xsi) {
  xsi[0] = -1.0;
  xsi[1] = 1.0;
}

void _e1c0_phi(double xsi, double *phi) {
  phi[0] = (1 - xsi) / 2.0;
  phi[1] = (1 + xsi) / 2.0;
}

void _e1c0_dphidx(double xsi, double *dphidxsi) {
  dphidxsi[0] = -0.5;
  dphidxsi[1] = 0.5;
}

femDiscrete *femDiscreteCreate(int n, femElementType type) {
  femDiscrete *theSpace = malloc(sizeof(femDiscrete));
  theSpace->type = type;
  theSpace->n = 0;
  theSpace->x = NULL;
  theSpace->phi = NULL;
  theSpace->dphidx = NULL;
  theSpace->x2 = NULL;
  theSpace->phi2 = NULL;
  theSpace->dphi2dx = NULL;

  if (type == FEM_QUAD && n == 4) {
    theSpace->n = 4;
    theSpace->x2 = _q1c0_x;
    theSpace->phi2 = _q1c0_phi;
    theSpace->dphi2dx = _q1c0_dphidx;
  } else if (type == FEM_TRIANGLE && n == 3) {
    theSpace->n = 3;
    theSpace->x2 = _p1c0_x;
    theSpace->phi2 = _p1c0_phi;
    theSpace->dphi2dx = _p1c0_dphidx;
  } else if (type == FEM_EDGE && n == 2) {
    theSpace->n = 2;
    theSpace->x = _e1c0_x;
    theSpace->phi = _e1c0_phi;
    theSpace->dphidx = _e1c0_dphidx;
  } else
    Error("Cannot create such a discrete space !");
  return theSpace;
}

void femDiscreteFree(femDiscrete *theSpace) { free(theSpace); }

void femDiscreteXsi2(femDiscrete *mySpace, double *xsi, double *eta) {
  mySpace->x2(xsi, eta);
}

void femDiscretePhi2(femDiscrete *mySpace, double xsi, double eta,
                     double *phi) {
  mySpace->phi2(xsi, eta, phi);
}

void femDiscreteDphi2(femDiscrete *mySpace, double xsi, double eta,
                      double *dphidxsi, double *dphideta) {
  mySpace->dphi2dx(xsi, eta, dphidxsi, dphideta);
}

void femDiscreteXsi(femDiscrete *mySpace, double *xsi) { mySpace->x(xsi); }

void femDiscretePhi(femDiscrete *mySpace, double xsi, double *phi) {
  mySpace->phi(xsi, phi);
}

void femDiscreteDphi(femDiscrete *mySpace, double xsi, double *dphidxsi) {
  mySpace->dphidx(xsi, dphidxsi);
}

void femDiscretePrint(femDiscrete *mySpace) {
  int i, j;
  int n = mySpace->n;
  double xsi[4], eta[4], phi[4], dphidxsi[4], dphideta[4];

  if (mySpace->type == FEM_EDGE) {
    femDiscreteXsi(mySpace, xsi);
    for (i = 0; i < n; i++) {
      femDiscretePhi(mySpace, xsi[i], phi);
      femDiscreteDphi(mySpace, xsi[i], dphidxsi);
      for (j = 0; j < n; j++) {
        printf("(xsi=%+.1f) : ", xsi[i]);
        printf(" phi(%d)=%+.1f", j, phi[j]);
        printf("   dphidxsi(%d)=%+.1f \n", j, dphidxsi[j]);
      }
      printf(" \n");
    }
  }

  if (mySpace->type == FEM_QUAD || mySpace->type == FEM_TRIANGLE) {
    femDiscreteXsi2(mySpace, xsi, eta);
    for (i = 0; i < n; i++) {
      femDiscretePhi2(mySpace, xsi[i], eta[i], phi);
      femDiscreteDphi2(mySpace, xsi[i], eta[i], dphidxsi, dphideta);
      for (j = 0; j < n; j++) {
        printf("(xsi=%+.1f,eta=%+.1f) : ", xsi[i], eta[i]);
        printf(" phi(%d)=%+.1f", j, phi[j]);
        printf("   dphidxsi(%d)=%+.1f", j, dphidxsi[j]);
        printf("   dphideta(%d)=%+.1f \n", j, dphideta[j]);
      }
      printf(" \n");
    }
  }
}

femFullSystem *femFullSystemCreate(int size) {
  femFullSystem *theSystem = malloc(sizeof(femFullSystem));
  femFullSystemAlloc(theSystem, size);
  femFullSystemInit(theSystem);

  return theSystem;
}

void femFullSystemFree(femFullSystem *theSystem) {
  free(theSystem->A);
  free(theSystem->B);
  free(theSystem);
}

void femFullSystemAlloc(femFullSystem *mySystem, int size) {
  int i;
  double *elem = malloc(sizeof(double) * size * (size + 1));
  mySystem->A = malloc(sizeof(double *) * size);
  mySystem->B = elem;
  mySystem->A[0] = elem + size;
  mySystem->size = size;
  for (i = 1; i < size; i++)
    mySystem->A[i] = mySystem->A[i - 1] + size;
}

void femFullSystemInit(femFullSystem *mySystem) {
  int i, size = mySystem->size;
  for (i = 0; i < size * (size + 1); i++)
    mySystem->B[i] = 0;
}

void femFullSystemPrint(femFullSystem *mySystem) {
  double **A, *B;
  int i, j, size;

  A = mySystem->A;
  B = mySystem->B;
  size = mySystem->size;

  for (i = 0; i < size; i++) {
    for (j = 0; j < size; j++)
      if (A[i][j] == 0)
        printf("         ");
      else
        printf(" %+.1e", A[i][j]);
    printf(" :  %+.1e \n", B[i]);
  }
}

double dotVectors(int size, double *A, double *B) {
  int i;
  double result = 0;
  for (i = 0; i < size; i++)
    result += A[i] * B[i];
  return result;
}

void dotScalarVector(int size, double A, double *B, double *C) {
  int i;
  for (i = 0; i < size; i++)
    C[i] = A * B[i];
}

void dotMatrixVector(int size, double **A, double *B, double *C) {
  int i, j;
  for (i = 0; i < size; i++) {
    C[i] = 0;
    for (j = 0; j < size; j++)
      C[i] += A[i][j] * B[j];
  }
}
void substracVector(int size, double *A, double *B, double *C) {
  int i;
  for (i = 0; i < size; i++)
    C[i] = A[i] - B[i];
}

void addVector(int size, double *A, double *B, double *C) {
  int i;
  for (i = 0; i < size; i++)
    C[i] = A[i] + B[i];
}

double vectorNorm(int size, double *A) {
  int i;
  double norm = 0;
  for (i = 0; i < size; i++)
    norm += A[i] * A[i];
  norm = sqrt(norm);
  return norm;
}

void incompleteCholesky(femFullSystem *mySystem) {
  int size = mySystem->size;
  double **A = mySystem->A;
  //double *B = mySystem->B;
  double **L = malloc(sizeof(double *) * size);
  for (int i = 0; i < size; i++) {
    L[i] = malloc(sizeof(double) * size);
    for (int j = 0; j < size; j++) {
      L[i][j] = 0;
    }
  }

  for (int i = 0; i < size; i++) {
    double sum = 0;
    for (int k = 0; k < i; k++) {
      sum += L[i][k] * L[i][k];
    }
    L[i][i] = sqrt(A[i][i] - sum);

    for (int j = i + 1; j < size; j++) {
      sum = 0;
      for (int k = 0; k < i; k++) {
        sum += L[i][k] * L[j][k];
      }
      L[j][i] = (A[j][i] - sum) / L[i][i];
    }
  }

  mySystem->L = L;
}

void applyPreconditionerJacob(femFullSystem *mySystem, double *r, double *z) {
  // M is the diagonal of A
  for (int i = 0; i < mySystem->size; i++) {
    z[i] = r[i] / mySystem->A[i][i];
  }
}

void applyPreconditionerCholesky(femFullSystem *mySystem, double *r,
                                 double *z) {
  // First perform forward substitution (L solve)
  for (int i = 0; i < mySystem->size; i++) {
    z[i] = r[i];
    for (int j = 0; j < i; j++) {
      z[i] -= mySystem->L[i][j] * z[j];
    }
    z[i] /= mySystem->L[i][i];
  }

  // Then perform backward substitution (Lᵀ solve)
  for (int i = mySystem->size - 1; i >= 0; i--) {
    for (int j = i + 1; j < mySystem->size; j++) {
      z[i] -= mySystem->L[j][i] * z[j]; // Note L[j][i] for Lᵀ
    }
    z[i] /= mySystem->L[i][i];
  }
}

double *femConjugateGradientJacob(femFullSystem *mySystem) {
  double **A, *B;

  A = mySystem->A;
  B = mySystem->B;
  int size = mySystem->size;
  double *Residual = malloc(sizeof(double) * size);
  double *PreconditionedResidual = malloc(sizeof(double) * size);
  double *search_direction = malloc(sizeof(double) * size);
  double *A_search_direction = malloc(sizeof(double) * size);
  double *utility = malloc(sizeof(double) * size);

  // Initialize residual vector
  dotMatrixVector(size, A, utility, Residual);
  substracVector(size, B, Residual, Residual);

  // Skip all if utility right guess
  if (vectorNorm(size, Residual) < TOL) {
    memcpy(B, utility, sizeof(double) * size);
    free(Residual);
    free(search_direction);
    free(A_search_direction);
    free(utility);
    return B;
  }

  // Apply preconditioner: z = M⁻¹r
  applyPreconditionerJacob(mySystem, Residual, PreconditionedResidual);

  // Initialize solution vector
  memcpy(B, utility, sizeof(double) * size);
  // Initialize search direction vector
  memcpy(search_direction, PreconditionedResidual, sizeof(double) * size);

  // Iterate until convergence
  while (1) {

    dotMatrixVector(size, A, search_direction, A_search_direction);
    double old_rz_product = dotVectors(size, Residual, PreconditionedResidual);
    double step_size =
        old_rz_product / dotVectors(size, search_direction, A_search_direction);

    // Update solution
    dotScalarVector(size, step_size, search_direction, utility);
    addVector(size, B, utility, B);

    // Update residual
    dotScalarVector(size, step_size, A_search_direction, utility);
    substracVector(size, Residual, utility, Residual);

    if (vectorNorm(size, Residual) < TOL) {
      break;
    }

    // Apply preconditioner: z = M⁻¹r
    applyPreconditionerJacob(mySystem, Residual, PreconditionedResidual);

    // Update search direction
    double beta =
        dotVectors(size, Residual, PreconditionedResidual) / old_rz_product;
    dotScalarVector(size, beta, search_direction, utility);
    addVector(size, PreconditionedResidual, utility, search_direction);
  }
  free(Residual);
  free(PreconditionedResidual);
  free(search_direction);
  free(A_search_direction);
  free(utility);
  return B;
}

double *femConjugateGradientCholesky(femFullSystem *mySystem) {
  incompleteCholesky(mySystem);
  double **A, *B;

  A = mySystem->A;
  B = mySystem->B;
  int size = mySystem->size;
  double *Residual = malloc(sizeof(double) * size);
  double *PreconditionedResidual = malloc(sizeof(double) * size);
  double *search_direction = malloc(sizeof(double) * size);
  double *A_search_direction = malloc(sizeof(double) * size);
  double *utility = malloc(sizeof(double) * size);

  // Initialize residual vector
  dotMatrixVector(size, A, utility, Residual);
  substracVector(size, B, Residual, Residual);

  // Skip all if utility right guess
  if (vectorNorm(size, Residual) < TOL) {
    memcpy(B, utility, sizeof(double) * size);
    free(Residual);
    free(search_direction);
    free(A_search_direction);
    free(utility);
    return B;
  }

  // Apply preconditioner: z = M⁻¹r
  applyPreconditionerCholesky(mySystem, Residual, PreconditionedResidual);

  // Initialize solution vector
  memcpy(B, utility, sizeof(double) * size);
  // Initialize search direction vector
  memcpy(search_direction, PreconditionedResidual, sizeof(double) * size);

  // Iterate until convergence
  while (1) {

    dotMatrixVector(size, A, search_direction, A_search_direction);
    double old_rz_product = dotVectors(size, Residual, PreconditionedResidual);
    double step_size =
        old_rz_product / dotVectors(size, search_direction, A_search_direction);

    // Update solution
    dotScalarVector(size, step_size, search_direction, utility);
    addVector(size, B, utility, B);

    // Update residual
    dotScalarVector(size, step_size, A_search_direction, utility);
    substracVector(size, Residual, utility, Residual);

    if (vectorNorm(size, Residual) < TOL) {
      break;
    }

    // Apply preconditioner: z = M⁻¹r
    applyPreconditionerCholesky(mySystem, Residual, PreconditionedResidual);

    // Update search direction
    double beta =
        dotVectors(size, Residual, PreconditionedResidual) / old_rz_product;
    dotScalarVector(size, beta, search_direction, utility);
    addVector(size, PreconditionedResidual, utility, search_direction);
  }
  free(Residual);
  free(PreconditionedResidual);
  free(search_direction);
  free(A_search_direction);
  free(utility);
  return B;
}

double *femConjugateGradient(femFullSystem *mySystem) {
  double **A, *B;

  A = mySystem->A;
  B = mySystem->B;
  int size = mySystem->size;
  double *Residual = malloc(sizeof(double) * size);
  double *search_direction = malloc(sizeof(double) * size);
  double *A_search_direction = malloc(sizeof(double) * size);
  double *utility = malloc(sizeof(double) * size);

  // Initialize residual vector
  dotMatrixVector(size, A, utility, Residual);
  substracVector(size, B, Residual, Residual);

  // Skip all if utility right guess
  if (vectorNorm(size, Residual) < TOL) {
    memcpy(B, utility, sizeof(double) * size);
    free(Residual);
    free(search_direction);
    free(A_search_direction);
    free(utility);
    return B;
  }
  // Initialize solution vector
  memcpy(B, utility, sizeof(double) * size);
  // Initialize search direction vector
  memcpy(search_direction, Residual, sizeof(double) * size);

  // Iterate until convergence
  while (1) {

    dotMatrixVector(size, A, search_direction, A_search_direction);
    double old_rz_product = dotVectors(size, Residual, Residual);
    double step_size =
        old_rz_product / dotVectors(size, search_direction, A_search_direction);

    // Update solution
    dotScalarVector(size, step_size, search_direction, utility);
    addVector(size, B, utility, B);

    // Update residual
    dotScalarVector(size, step_size, A_search_direction, utility);
    substracVector(size, Residual, utility, Residual);

    if (vectorNorm(size, Residual) < TOL) {
      break;
    }

    // Update search direction
    double beta = dotVectors(size, Residual, Residual) / old_rz_product;
    dotScalarVector(size, beta, search_direction, utility);
    addVector(size, Residual, utility, search_direction);
  }
  free(Residual);
  free(search_direction);
  free(A_search_direction);
  free(utility);
  return B;
}

double *femFullSystemEliminate(femFullSystem *mySystem) {
  double **A, *B, factor;
  int i, j, k, size;

  A = mySystem->A;
  B = mySystem->B;
  size = mySystem->size;

  /* Gauss elimination */

  for (k = 0; k < size; k++) {
    if (fabs(A[k][k]) <= 1e-16) {
      printf("Pivot index %d  ", k);
      printf("Pivot value %e  ", A[k][k]);
      Error("Cannot eliminate with such a pivot");
    }
    for (i = k + 1; i < size; i++) {
      factor = A[i][k] / A[k][k];
      for (j = k + 1; j < size; j++)
        A[i][j] = A[i][j] - A[k][j] * factor;
      B[i] = B[i] - B[k] * factor;
    }
  }

  /* Back-substitution */

  for (i = size - 1; i >= 0; i--) {
    factor = 0;
    for (j = i + 1; j < size; j++)
      factor += A[i][j] * B[j];
    B[i] = (B[i] - factor) / A[i][i];
  }

  return (mySystem->B);
}

void femFullSystemConstrain(femFullSystem *mySystem, int myNode,
                            double myValue) {
  double **A, *B;
  int i, size;

  A = mySystem->A;
  B = mySystem->B;
  size = mySystem->size;

  for (i = 0; i < size; i++) {
    B[i] -= myValue * A[i][myNode];
    A[i][myNode] = 0;
  }

  for (i = 0; i < size; i++)
    A[myNode][i] = 0;

  A[myNode][myNode] = 1;
  B[myNode] = myValue;
}

femProblem *femElasticityCreate(femGeo *theGeo, double E, double nu, double rho,
                                double g, femElasticCase iCase) {
  femProblem *theProblem = malloc(sizeof(femProblem));
  theProblem->E = E;
  theProblem->nu = nu;
  theProblem->g = g;
  theProblem->rho = rho;

  if (iCase == PLANAR_STRESS) {
    theProblem->A = E / (1 - nu * nu);
    theProblem->B = E * nu / (1 - nu * nu);
    theProblem->C = E / (2 * (1 + nu));
  } else if (iCase == PLANAR_STRAIN || iCase == AXISYM) {
    theProblem->A = E * (1 - nu) / ((1 + nu) * (1 - 2 * nu));
    theProblem->B = E * nu / ((1 + nu) * (1 - 2 * nu));
    theProblem->C = E / (2 * (1 + nu));
  }

  theProblem->planarStrainStress = iCase;
  theProblem->nBoundaryConditions = 0;
  theProblem->conditions = NULL;

  int size = 2 * theGeo->theNodes->nNodes;
  theProblem->constrainedNodes = malloc(size * sizeof(int));
  theProblem->soluce = malloc(size * sizeof(double));
  theProblem->residuals = malloc(size * sizeof(double));
  for (int i = 0; i < size; i++) {
    theProblem->constrainedNodes[i] = -1;
    theProblem->soluce[i] = 0.0;
    theProblem->residuals[i] = 0.0;
  }

  theProblem->geometry = theGeo;
  if (theGeo->theElements->nLocalNode == 3) {
    theProblem->space = femDiscreteCreate(3, FEM_TRIANGLE);
    theProblem->rule = femIntegrationCreate(3, FEM_TRIANGLE);
  }
  if (theGeo->theElements->nLocalNode == 4) {
    theProblem->space = femDiscreteCreate(4, FEM_QUAD);
    theProblem->rule = femIntegrationCreate(4, FEM_QUAD);
  }
  theProblem->spaceEdge = femDiscreteCreate(2, FEM_EDGE);
  theProblem->ruleEdge = femIntegrationCreate(2, FEM_EDGE);
  theProblem->system = femFullSystemCreate(size);

  // femDiscretePrint(theProblem->space);
  // femDiscretePrint(theProblem->spaceEdge);

  return theProblem;
}

void femElasticityFree(femProblem *theProblem) {
  femFullSystemFree(theProblem->system);
  femIntegrationFree(theProblem->rule);
  femDiscreteFree(theProblem->space);
  femIntegrationFree(theProblem->ruleEdge);
  femDiscreteFree(theProblem->spaceEdge);
  free(theProblem->conditions);
  free(theProblem->constrainedNodes);
  free(theProblem->soluce);
  free(theProblem->residuals);
  free(theProblem);
}

void femElasticityAddBoundaryCondition(femProblem *theProblem, char *nameDomain,
                                       femBoundaryType type, double value) {
  int iDomain = geoGetDomain(nameDomain);
  if (iDomain == -1)
    Error("Undefined domain :-(");

  femBoundaryCondition *theBoundary = malloc(sizeof(femBoundaryCondition));
  theBoundary->domain = theProblem->geometry->theDomains[iDomain];
  theBoundary->value = value;
  theBoundary->type = type;
  theProblem->nBoundaryConditions++;
  int size = theProblem->nBoundaryConditions;

  if (theProblem->conditions == NULL)
    theProblem->conditions = malloc(size * sizeof(femBoundaryCondition *));
  else
    theProblem->conditions =
        realloc(theProblem->conditions, size * sizeof(femBoundaryCondition *));
  theProblem->conditions[size - 1] = theBoundary;

  int shift = -1;
  if (type == DIRICHLET_X)
    shift = 0;
  if (type == DIRICHLET_Y)
    shift = 1;
  if (shift == -1)
    return;
  int *elem = theBoundary->domain->elem;
  int nElem = theBoundary->domain->nElem;
  for (int e = 0; e < nElem; e++) {
    for (int i = 0; i < 2; i++) {
      int node = theBoundary->domain->mesh->elem[2 * elem[e] + i];
      theProblem->constrainedNodes[2 * node + shift] = size - 1;
    }
  }
}


double femElasticityIntegrate(femProblem *theProblem,
                              double (*f)(double x, double y)) {
  femIntegration *theRule = theProblem->rule;
  femGeo *theGeo = theProblem->geometry;
  femNodes *theNodes = theGeo->theNodes;
  femMesh *theMesh = theGeo->theElements;
  femDiscrete *theSpace = theProblem->space;

  double x[4], y[4], phi[4], dphidxsi[4], dphideta[4];
  int iElem, iInteg, i, map[4];
  int nLocal = theMesh->nLocalNode;
  double value = 0.0;
  for (iElem = 0; iElem < theMesh->nElem; iElem++) {
    for (i = 0; i < nLocal; i++) {
      map[i] = theMesh->elem[iElem * nLocal + i];
      x[i] = theNodes->X[map[i]];
      y[i] = theNodes->Y[map[i]];
    }
    for (iInteg = 0; iInteg < theRule->n; iInteg++) {
      double xsi = theRule->xsi[iInteg];
      double eta = theRule->eta[iInteg];
      double weight = theRule->weight[iInteg];
      femDiscretePhi2(theProblem->space, xsi, eta, phi);
      femDiscreteDphi2(theSpace, xsi, eta, dphidxsi, dphideta);
      double dxdxsi = 0.0;
      double dxdeta = 0.0;
      double dydxsi = 0.0;
      double dydeta = 0.0;
      for (i = 0; i < theSpace->n; i++) {
        dxdxsi += x[i] * dphidxsi[i];
        dxdeta += x[i] * dphideta[i];
        dydxsi += y[i] * dphidxsi[i];
        dydeta += y[i] * dphideta[i];
      }
      double jac = fabs(dxdxsi * dydeta - dxdeta * dydxsi);
      for (i = 0; i < theProblem->space->n; i++) {
        value += phi[i] * f(x[i], y[i]) * jac * weight;
      }
    }
  }
  return value;
}

double femMin(double *x, int n) {
  double myMin = x[0];
  int i;
  for (i = 1; i < n; i++)
    myMin = fmin(myMin, x[i]);
  return myMin;
}

double femMax(double *x, int n) {
  double myMax = x[0];
  int i;
  for (i = 1; i < n; i++)
    myMax = fmax(myMax, x[i]);
  return myMax;
}

void femError(char *text, int line, char *file) {
  printf("\n-------------------------------------------------------------------"
         "------------- ");
  printf("\n  Error in %s at line %d : \n  %s\n", file, line, text);
  printf("---------------------------------------------------------------------"
         " Yek Yek !! \n\n");
  exit(69);
}

void femErrorGmsh(int ierr, int line, char *file) {
  if (ierr == 0)
    return;
  printf("\n-------------------------------------------------------------------"
         "------------- ");
  printf("\n  Error in %s at line %d : \n  error code returned by gmsh %d\n",
         file, line, ierr);
  printf("---------------------------------------------------------------------"
         " Yek Yek !! \n\n");
  exit(69);
}

void femErrorScan(int test, int line, char *file) {
  if (test >= 0)
    return;

  printf("\n-------------------------------------------------------------------"
         "------------- ");
  printf("\n  Error in fscanf or fgets in %s at line %d : \n", file, line);
  printf("---------------------------------------------------------------------"
         " Yek Yek !! \n\n");
  exit(69);
}

void femWarning(char *text, int line, char *file) {
  printf("\n-------------------------------------------------------------------"
         "------------- ");
  printf("\n  Warning in %s at line %d : \n  %s\n", file, line, text);
  printf("---------------------------------------------------------------------"
         " Yek Yek !! \n\n");
}


// CSR PART OF THE CODE

void femCsrSystemFree(femCsrSystem *theSystem) {
  free(theSystem->data);
  free(theSystem->column);
  free(theSystem->row_pointer);
  free(theSystem->B);
  free(theSystem);
}

femCsrSystem *femCsrSystemCreate(int size, femGeo *theGeo) {
  femCsrSystem *theSystem = malloc(sizeof(femCsrSystem));
  femCsrSystemAlloc(theSystem, theGeo,size);
  femCsrSystemInit(theSystem, theGeo);

  return theSystem;
}

void femCsrSystemAlloc(femCsrSystem *mySystem, femGeo *theGeo, int size) {
  int number_internal_links = (theGeo->theElements->nElem * theGeo->theElements->nLocalNode - theGeo->theEdges->nElem) / 2;
  int links_per_element = (theGeo->theElements->nLocalNode - 1) * theGeo->theElements->nLocalNode / 2;
  mySystem->data = malloc(sizeof(double) * (3 * size / 2 + 4 * (links_per_element * theGeo->theElements->nElem - number_internal_links))); 
  //size for the diagonal, 4 bc eatch entry is a 2x2 matrix
  mySystem->column = malloc(sizeof(int) * (3 * size / 2 + 4 * (links_per_element * theGeo->theElements->nElem - number_internal_links)));
  mySystem->row_pointer = malloc(sizeof(int) * (size + 1));
  //mySystem->row_pointer[0] = 0;
  //mySystem->row_pointer[size] = size + 4 * number_individual_links; // init here to avoid redoing the math for th number of individual links
  mySystem->B = malloc(sizeof(double) * size);
  mySystem->size = size;
}

void femCsrSystemInit(femCsrSystem *mySystem, femGeo *theGeo) {
  int size = mySystem->size;
  int nElem = theGeo->theElements->nElem;
  int nlocal = theGeo->theElements->nLocalNode;
  int buffer_size = size / 2; // size chosen arbitrarily
  int* row_pointer = mySystem->row_pointer;
  int* column = mySystem->column;
  int* elem = theGeo->theElements->elem;
  // will store all the other nodes already linked to the node
  int* buffer = malloc( buffer_size * sizeof(int)); 

  // initialize the row pointer
  row_pointer[0] = 0;


  int found_matches = 0;
  for (int i = 0; i < size / 2; i++) {
    for (int j = 0; j < nElem; j++) {
      for (int k = 0; k < nlocal; k++) {
        int node = elem[j * nlocal + k];
        if (node == i) {
          // found the node
          for (int l = 0; l < nlocal; l++) {
            int other_node = elem[j * nlocal + l];
            // check if the other node is not the same as the node or link already exists
            if (other_node != node && other_node > i) {
              // check if the node is already in the buffer 
              int already_found = 0;
              for (int m = 0; m < found_matches; m++) {
                if (buffer[m] == other_node) {
                  already_found = 1;
                  break;
                }
              }
              if (!already_found) {
                buffer[found_matches] = other_node;
                found_matches++;
              }
            }
          }
        }
      }
    }
    row_pointer[i * 2 + 1] = found_matches * 2 + row_pointer[2 * i] + 2;
    row_pointer[i * 2 + 2] = found_matches * 2 + row_pointer[2 * i + 1] + 1;
    // itself
    column[row_pointer[2 * i]] = 2 * i;
    column[row_pointer[2 * i] + 1] = 2 * i + 1;
    column[row_pointer[2 * i + 1]] = 2 * i + 1;
    for (int j = 0; j < found_matches; j++) {
      column[2 + row_pointer[2 * i] + 2 * j] = 2 * buffer[j];
      column[2 + row_pointer[2 * i] + 2 * j + 1] = 2 * buffer[j] + 1;
      column[1 + row_pointer[2 * i + 1] + 2 * j] = 2 * buffer[j];
      column[1 + row_pointer[2 * i + 1] + 2 * j + 1] = 2 * buffer[j] + 1;
    }
    for (int j = 0; j < buffer_size; j++) { // not sorted but idcare
      buffer[j] = -1;
    }
    found_matches = 0;
  }

  for (int i = 0; i < size; i++)
    mySystem->B[i] = 0;

  free(buffer);
  int number_internal_links = (theGeo->theElements->nElem * theGeo->theElements->nLocalNode - theGeo->theEdges->nElem) / 2;
  int links_per_element = (theGeo->theElements->nLocalNode - 1) * theGeo->theElements->nLocalNode / 2;
  if (mySystem->row_pointer[size] != (3 * size / 2 + 4 * (links_per_element * theGeo->theElements->nElem - number_internal_links))) {
    printf("scream\n"); // fucking debugger
    printf("row_pointer[%d] = %d\n", size, mySystem->row_pointer[size]);
    printf("alocated size = %d\n", (3 * size / 2 + 4 * (links_per_element * theGeo->theElements->nElem - number_internal_links)));
  }
}


void femElasticityAddBoundaryConditionCsr(femProblemCsr *theProblem, char *nameDomain,
                                       femBoundaryType type, double value) {
  int iDomain = geoGetDomain(nameDomain);
  if (iDomain == -1)
    Error("Undefined domain :-(");

  femBoundaryCondition *theBoundary = malloc(sizeof(femBoundaryCondition));
  theBoundary->domain = theProblem->geometry->theDomains[iDomain];
  theBoundary->value = value;
  theBoundary->type = type;
  theProblem->nBoundaryConditions++;
  int size = theProblem->nBoundaryConditions;

  if (theProblem->conditions == NULL)
    theProblem->conditions = malloc(size * sizeof(femBoundaryCondition *));
  else
    theProblem->conditions =
        realloc(theProblem->conditions, size * sizeof(femBoundaryCondition *));
  theProblem->conditions[size - 1] = theBoundary;

  int shift = -1;
  if (type == DIRICHLET_X)
    shift = 0;
  if (type == DIRICHLET_Y)
    shift = 1;
  if (shift == -1)
    return;
  int *elem = theBoundary->domain->elem;
  int nElem = theBoundary->domain->nElem;
  for (int e = 0; e < nElem; e++) {
    for (int i = 0; i < 2; i++) {
      int node = theBoundary->domain->mesh->elem[2 * elem[e] + i];
      theProblem->constrainedNodes[2 * node + shift] = size - 1;
    }
  }
}

void femCsrSystemConstrain(femCsrSystem *mySystem, int myNode,
                            double myValue) {
  double *data, *B;
  int i, size, *column, *row_pointer;

  data = mySystem->data;
  column = mySystem->column;
  row_pointer = mySystem->row_pointer;
  B = mySystem->B;
  size = mySystem->size;

  // upper triangle
  for (i = 0; i < myNode; i++) {
    for (int j = row_pointer[i]; j < row_pointer[i + 1]; j++) {
      if (column[j] == myNode) {
        // Found the column index for the current row
        B[i] -= myValue * data[j];
        data[j] = 0;
        break;;
      }
    }
  }

  // lower triangle
  for (int j = row_pointer[myNode]; j < row_pointer[myNode + 1]; j++) {
    B[column[j]] -= myValue * data[j];
  }

  for (i = row_pointer[myNode]; i < row_pointer[myNode + 1]; i++) {
    data[i] = 0;
    if (column[i] == myNode) {
      // Found the column index for the current row
      data[i] = 1;
      continue;
    }
  }

  B[myNode] = myValue;
}

femProblemCsr *femElasticityCreateCsr(femGeo *theGeo, double E, double nu, double rho,
                                double g, femElasticCase iCase) {
  femProblemCsr *theProblem = malloc(sizeof(femProblemCsr));
  theProblem->E = E;
  theProblem->nu = nu;
  theProblem->g = g;
  theProblem->rho = rho;

  if (iCase == PLANAR_STRESS) {
    theProblem->A = E / (1 - nu * nu);
    theProblem->B = E * nu / (1 - nu * nu);
    theProblem->C = E / (2 * (1 + nu));
  } else if (iCase == PLANAR_STRAIN || iCase == AXISYM) {
    theProblem->A = E * (1 - nu) / ((1 + nu) * (1 - 2 * nu));
    theProblem->B = E * nu / ((1 + nu) * (1 - 2 * nu));
    theProblem->C = E / (2 * (1 + nu));
  }

  theProblem->planarStrainStress = iCase;
  theProblem->nBoundaryConditions = 0;
  theProblem->conditions = NULL;

  int size = 2 * theGeo->theNodes->nNodes;
  theProblem->constrainedNodes = malloc(size * sizeof(int));
  theProblem->soluce = malloc(size * sizeof(double));
  theProblem->residuals = malloc(size * sizeof(double));
  for (int i = 0; i < size; i++) {
    theProblem->constrainedNodes[i] = -1;
    theProblem->soluce[i] = 0.0;
    theProblem->residuals[i] = 0.0;
  }

  theProblem->geometry = theGeo;
  if (theGeo->theElements->nLocalNode == 3) {
    theProblem->space = femDiscreteCreate(3, FEM_TRIANGLE);
    theProblem->rule = femIntegrationCreate(3, FEM_TRIANGLE);
  }
  if (theGeo->theElements->nLocalNode == 4) {
    theProblem->space = femDiscreteCreate(4, FEM_QUAD);
    theProblem->rule = femIntegrationCreate(4, FEM_QUAD);
  }
  theProblem->spaceEdge = femDiscreteCreate(2, FEM_EDGE);
  theProblem->ruleEdge = femIntegrationCreate(2, FEM_EDGE);
  theProblem->systemCsr = femCsrSystemCreate(size, theGeo);

  // femDiscretePrint(theProblem->space);
  // femDiscretePrint(theProblem->spaceEdge);

  return theProblem;
}

void femElasticityFreeCsr(femProblemCsr *theProblem) {
  femCsrSystemFree(theProblem->systemCsr);
  femIntegrationFree(theProblem->rule);
  femDiscreteFree(theProblem->space);
  femIntegrationFree(theProblem->ruleEdge);
  femDiscreteFree(theProblem->spaceEdge);
  free(theProblem->conditions);
  free(theProblem->constrainedNodes);
  free(theProblem->soluce);
  free(theProblem->residuals);
  free(theProblem);
}

void dotMatrixVectorCsr(int size, double *data, int* column, int* row_pointer, double *B, double *result) {
  for (int i = 0; i < size; i++) {
    result[i] = 0;
  }

  for (int i = 0; i < size; i++) {
    for (int j = row_pointer[i] + 1; j < row_pointer[i + 1]; j++) {
      result[i] += data[j] * B[column[j]];
      result[column[j]] += data[j] * B[i];
    }
  }

  // diagonal
  for (int i = 0; i < size; i++) {
    result[i] += data[row_pointer[i]] * B[i];
  }
}

double *femConjugateGradientCsr(femCsrSystem *mySystem) {
  double *B;
  double *data = mySystem->data;
  int *column = mySystem->column;
  int *row_pointer = mySystem->row_pointer;
  B = mySystem->B;
  int size = mySystem->size;
  double *Residual = malloc(sizeof(double) * size);
  double *search_direction = malloc(sizeof(double) * size);
  double *A_search_direction = malloc(sizeof(double) * size);
  double *utility = malloc(sizeof(double) * size);

  double *testing_array = malloc(sizeof(double) * size);

  // Initialize residual vector
  dotMatrixVectorCsr(size, data, column, row_pointer, utility, Residual);
  substracVector(size, B, Residual, Residual);

  // Skip all if utility right guess
  if (vectorNorm(size, Residual) < TOL) {
    memcpy(B, utility, sizeof(double) * size);
    free(Residual);
    free(search_direction);
    free(A_search_direction);
    free(utility);
    return B;
  }
  // Initialize solution vector
  memcpy(B, utility, sizeof(double) * size);
  // Initialize search direction vector
  memcpy(search_direction, Residual, sizeof(double) * size);

  // Iterate until convergence
  while (1) {
    dotMatrixVectorCsr(size, data, column, row_pointer, search_direction, A_search_direction);

    double old_rz_product = dotVectors(size, Residual, Residual);
    double step_size =
        old_rz_product / dotVectors(size, search_direction, A_search_direction);

    // Update solution
    dotScalarVector(size, step_size, search_direction, utility);
    addVector(size, B, utility, B);

    // Update residual
    dotScalarVector(size, step_size, A_search_direction, utility);
    substracVector(size, Residual, utility, Residual);

    if (vectorNorm(size, Residual) < TOL) {
      break;
    }

    // Update search direction
    double beta = dotVectors(size, Residual, Residual) / old_rz_product;
    dotScalarVector(size, beta, search_direction, utility);
    addVector(size, Residual, utility, search_direction);
  }
  free(Residual);
  free(search_direction);
  free(A_search_direction);
  free(utility);
  free(testing_array);
  return B;
}