
/*
 *  fem.h
 *  Library for LEPL1110 : Finite Elements for dummies
 *
 *  Copyright (C) 2024 UCL-IMMC : Vincent Legat
 *  All rights reserved.
 *
 */

#ifndef _FEM_H_
#define _FEM_H_

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#define ErrorScan(a)   femErrorScan(a,__LINE__,__FILE__)
#define Error(a)       femError(a,__LINE__,__FILE__)
#define Warning(a)     femWarning(a,  __LINE__, __FILE__)
#define FALSE 0 
#define TRUE  1
#define MAXNAME 256

typedef enum {FEM_TRIANGLE,FEM_QUAD} femElementType;

// A structure representing an array of nodes
// X and Y are the coordinates of the nodes
// The number of nodes is stored in nNodes
typedef struct {
    int nNodes;
    double *X;
    double *Y;
} femNodes;

// A structure containing the mesh
// nLocalNode is the number of nodes per element
// nElem is the number of elements
// elem
// nodes is a pointer to a femNodes a structure representing an array of nodes
typedef struct {
    int nLocalNode;
    int nElem;
    int *elem;
    femNodes *nodes;
} femMesh;

// Boh jsp
typedef struct {
    femMesh *mesh;
    int nElem;
    int *elem;
    char name[MAXNAME];
} femDomain;

// A structure representing the geometry of the problem
typedef struct {
    femElementType elementType;
    femNodes *theNodes;
    femMesh  *theElements;
    femMesh  *theEdges;
    int nDomains;
    femDomain **theDomains;
} femGeo;


// A structure containting functions to compute 
typedef struct {
    int n;
    void (*x2)(double *xsi, double *eta);
    void (*phi2)(double xsi, double eta, double *phi);
    void (*dphi2dx)(double xsi, double eta, double *dphidxsi, double *dphideta);
} femDiscrete;

// A structure containing the integration rule
// n is the number of integration points
// xsi, eta and weight are the coordinates and the weights of the integration points
typedef struct {
    int n;
    const double *xsi;
    const double *eta;
    const double *weight;
} femIntegration;

// A structure containing the matrix A and the vector B of the linear system
typedef struct {
    double *B;
    double **A;
    int size;
} femFullSystem;

// A structure containing the geometry, the discrete space, the integration rule and the linear system
typedef struct {
    femGeo *geo;
    femDiscrete *space;
    femIntegration *rule;
    femFullSystem *system;
} femPoissonProblem;


femGeo*              geoMeshCreate(const char *filename);
void                 geoMeshFree(femGeo* theGeometry);
void                 geoMeshPrint(femGeo* theGeometry);
femDomain*           geoGetDomain(femGeo* theGeometry, char *name);

femIntegration      *femIntegrationCreate(int n, femElementType type);
void                 femIntegrationFree(femIntegration *theRule);

femDiscrete*         femDiscreteCreate(int n, femElementType type);
void                 femDiscreteFree(femDiscrete* mySpace);
void                 femDiscretePrint(femDiscrete* mySpace);
void                 femDiscreteXsi2(femDiscrete* mySpace, double *xsi, double *eta);
void                 femDiscretePhi2(femDiscrete* mySpace, double xsi, double eta, double *phi);
void                 femDiscreteDphi2(femDiscrete* mySpace, double xsi, double eta, double *dphidxsi, double *dphideta);

femFullSystem*       femFullSystemCreate(int size);
void                 femFullSystemFree(femFullSystem* mySystem);
void                 femFullSystemPrint(femFullSystem* mySystem);
void                 femFullSystemInit(femFullSystem* mySystem);
void                 femFullSystemAlloc(femFullSystem* mySystem, int size);
double*              femFullSystemEliminate(femFullSystem* mySystem);
void                 femFullSystemConstrain(femFullSystem* mySystem, int myNode, double value);

femPoissonProblem   *femPoissonCreate(const char *filename);
void                 femPoissonFindBoundaryNodes(femPoissonProblem *theProblem);
void                 femPoissonLocal(femPoissonProblem *theProblem, const int i, int *map, double *x, double *y);
void                 femPoissonFree(femPoissonProblem *theProblem);
void                 femPoissonSolve(femPoissonProblem *theProblem);

double               femMin(double *x, int n);
double               femMax(double *x, int n);
void                 femError(char *text, int line, char *file);
void                 femErrorScan(int test, int line, char *file);
void                 femWarning(char *text, int line, char *file);


#endif