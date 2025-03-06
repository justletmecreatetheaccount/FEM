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

/**
 * @brief Structure containing the nodes of the global mesh
 * @param nNodes Number of nodes
 * @param X x-coordinates of the nodes
 * @param Y y-coordinates of the nodes
 */
typedef struct {
    int nNodes;
    double *X;
    double *Y;
} femNodes;


/**
 * @brief Structure containing the elements of a mesh. Typically there are 2 meshes in a project: 
 * one for the edges (segments) and one for the elements inside (e.g. triangles or quads)
 * @param nLocalNode Number of nodes per element, e.g. 2 for a segment (e.g. elements are edges), 3 for a triangle, 4 for a quad
 * @param nElem Number of elements 
 * @param elem List of elements' nodes indices (size: `nElem * nLocalNode`). 
 * It is arranged as follows: `elem[i*nLocalNode + j]` is the index of the j-th node of the i-th element.
 * For example, for a triangle, `elem[3*i]`, `elem[3*i+1]` and `elem[3*i+2]` are the indices of the 3 nodes of the i-th triangle
 * @param nodes Pointer to the nodes' structure
 */
typedef struct {
    int nLocalNode;
    int nElem;
    int *elem;
    femNodes *nodes;
} femMesh;


/**
 * @brief Structure containing one domain (a boundary) of the general geometry (e.g. the left side of a rectangle)
 * @param mesh Pointer to the mesh structure: contains elements of 2 nodes (edges)
 * @param nElem Number of elements in the domain
 * @param elem List of domain's elements' indices
 * @param name Name of the domain
 */
typedef struct {
    femMesh *mesh;
    int nElem;
    int *elem;
    char name[MAXNAME];
} femDomain;

/**
 * @brief Structure containing the general geometry of the problem
 * @param elementType Type of elements used in the mesh (e.g. triangles or quads)
 * @param theNodes Pointer to the nodes structure containing the nodes of the global mesh
 * @param theElements Pointer to the elements structure
 * @param theEdges Pointer to the edges structure
 * @param nDomains Number of domains in the geometry
 * @param theDomains List of domains
 */
typedef struct {
    femElementType elementType;
    femNodes *theNodes;
    femMesh  *theElements;
    femMesh  *theEdges;
    int nDomains;
    femDomain **theDomains;
} femGeo;


/**
 * @brief Structure containing the discrete integration rule for a given element type (e.g. triangles or quads)
 * @param n Number of integration points (e.g. 3 for a triangle)
 * @param x2 Function that computes the coordinates of the integration points
 * @param phi2 Function that computes the shape functions at a given integration point
 * @param dphi2dx Function that computes the derivatives of the shape functions at a given integration point
 */
typedef struct {
    int n;
    /**
     * @brief Function that computes the coordinates of the integration points

     * @param xsi array of size n to store the xsi coordinates of the integration points
     * @param eta array of size n to store the eta coordinates of the integration points
     */
    void (*x2)(double *xsi, double *eta);

    /**
     * @brief Function that computes the shape functions at a given integration point
     * @param xsi xsi coordinate of the integration point
     * @param eta eta coordinate of the integration point
     * @param phi array of size n to store the n shape functions evaluated at the integration point
     */
    void (*phi2)(double xsi, double eta, double *phi);

    /**
     * @brief Function that computes the derivatives of the shape functions at a given integration point
     * @param xsi xsi coordinate of the integration point
     * @param eta eta coordinate of the integration point
     * @param dphidxsi array of size n to store the n derivatives of the shape functions with respect to xsi
     * @param dphideta array of size n to store the n derivatives of the shape functions with respect to eta
     */
    void (*dphi2dx)(double xsi, double eta, double *dphidxsi, double *dphideta);
} femDiscrete;
    
/**
 * @brief Structure containing the integration rule for a given element type (e.g. triangles or quads)
 * @param n Number of integration points (e.g. 4 for a quad)
 * @param xsi Array of size n containing the xsi coordinates of the integration points
 * @param eta Array of size n containing the eta coordinates of the integration points
 * @param weight Array of size n containing the weights of the integration points
 */
typedef struct {
    int n;
    const double *xsi;
    const double *eta;
    const double *weight;
} femIntegration;

/**
 * @brief Structure containing the full system of equations
 * @param B Right-hand side of the system
 * @param A Stiffness matrix of the system
 * @param size Size of the system
 */
typedef struct {
    double *B;
    double **A;
    int size;
} femFullSystem;

/**
 * @brief Structure containing the problem to solve
 * @param geo Pointer to the geometry structure
 * @param space Pointer to the discrete structure
 * @param rule Pointer to the integration rule structure
 * @param system Pointer to the full system structure
 */
typedef struct {
    femGeo *geo;
    femDiscrete *space;
    femIntegration *rule;
    femFullSystem *system;
} femPoissonProblem;


#pragma region Geometry
/**
 * @brief Read the mesh of the geometry from a file
 * @param filename Name of the file to read the mesh from
 */
femGeo*              geoMeshCreate(const char *filename);

/**
 * @brief Free the geometry. Internal function.
 */
void                 geoMeshFree(femGeo* theGeometry);

/**
 * @brief Print the geometry's mesh's information
 */
void                 geoMeshPrint(femGeo* theGeometry);

/**
 * @brief Get the domain index from its name
 * @param name Name of the domain
 * @return pointer to the domain structure if found, NULL otherwise
 */
femDomain*           geoGetDomain(femGeo* theGeometry, char *name);
#pragma endregion

#pragma region Integration

/**
 * @brief Initialize the integration rule for a given element type
 * @param n Number of integration points
 * @param type Type of element (e.g. triangles or quads)
 * @return Pointer to the integration rule structure
 */
femIntegration      *femIntegrationCreate(int n, femElementType type);

/**
 * @brief Free the integration rule
 * @param theRule Pointer to the integration rule structure
 */
void                 femIntegrationFree(femIntegration *theRule);
#pragma endregion

#pragma region Discrete

/**
 * @brief Create a new discrete structure for a given element type
 * @param n Number of nodes per element
 * @param type Type of element (e.g. triangles or quads)
 * @return Pointer to the discrete structure
 */
femDiscrete*         femDiscreteCreate(int n, femElementType type);

/**
 * @brief Free the discrete structure
 * @param mySpace Pointer to the discrete structure
 */
void                 femDiscreteFree(femDiscrete* mySpace);

/**
 * @brief Print the discrete structure
 * @param mySpace Pointer to the discrete structure
 */
void                 femDiscretePrint(femDiscrete* mySpace);

/**
 * @brief Compute the coordinates of the integration points
 * @param mySpace Pointer to the discrete structure
 * @param xsi Array of size n to store the xsi coordinates of the integration points
 * @param eta Array of size n to store the eta coordinates of the integration points
 */
void                 femDiscreteXsi2(femDiscrete* mySpace, double *xsi, double *eta);

/**
 * @brief Compute the shape functions at a given integration point
 * @param mySpace Pointer to the discrete structure
 * @param xsi xsi coordinate of the integration point
 * @param eta eta coordinate of the integration point
 * @param phi Array of size n to store the n shape functions evaluated at the integration point
 */
void                 femDiscretePhi2(femDiscrete* mySpace, double xsi, double eta, double *phi);

/**
 * @brief Compute the derivatives of the shape functions at a given integration point
 * @param mySpace Pointer to the discrete structure
 * @param xsi xsi coordinate of the integration point
 * @param eta eta coordinate of the integration point
 * @param dphidxsi Array of size n to store the n derivatives of the shape functions with respect to xsi
 * @param dphideta Array of size n to store the n derivatives of the shape functions with respect to eta
 */
void                 femDiscreteDphi2(femDiscrete* mySpace, double xsi, double eta, double *dphidxsi, double *dphideta);
#pragma endregion

#pragma region FullSystem

/**
 * @brief Create a new system of equations
 * @param size Size of the system
 * @return Pointer to the full system structure
 */
femFullSystem*       femFullSystemCreate(int size);

/**
 * @brief Free the full system of equations
 * @param mySystem Pointer to the full system structure
 */
void                 femFullSystemFree(femFullSystem* mySystem);

/**
 * @brief Print the full system of equations
 * @param mySystem Pointer to the full system structure
 */
void                 femFullSystemPrint(femFullSystem* mySystem);

/**
 * @brief Initialize the full system of equations to 0 (A = 0, B = 0)
 * @param mySystem Pointer to the full system structure
 */
void                 femFullSystemInit(femFullSystem* mySystem);

/**
 * @brief Allocate memory for the full system of equations
 * @param mySystem Pointer to the full system structure
 * @param size Size of the system
 */
void                 femFullSystemAlloc(femFullSystem* mySystem, int size);

/**
 * @brief Eliminate the full system of equations using Gauss elimination
 * @param mySystem Pointer to the full system structure
 * @return Pointer to the solution
 */
double*              femFullSystemEliminate(femFullSystem* mySystem);

/**
 * @brief Constrain a node in the full system of equations
 * @param mySystem Pointer to the full system structure
 * @param myNode Index of the node to constrain
 * @param value Value to constrain the node to
 */
void                 femFullSystemConstrain(femFullSystem* mySystem, int myNode, double value);
#pragma endregion

#pragma region Poisson

/**
 * @brief Read the geometry from a file and create a new problem structure
 * @param filename Name of the file to read the problem from
 * @return Pointer to the problem structure
 */
femPoissonProblem   *femPoissonCreate(const char *filename);

/**
 * @brief Create a new domain structure in the problem whose elements are the nodes on the boundary (stored as elements of 1 node)
 * See an [example](https://hokkaydo.be/notes/objective_poisson.png) of what result is expected
 * @param theProblem Pointer to the problem structure in which to create the domain
 */
void                 femPoissonFindBoundaryNodes(femPoissonProblem *theProblem);

/**
 * @brief Copies the local indices and coordinates of the i-th element's nodes in the `map`, `x` and `y` arrays
 * 
 * @param theProblem Pointer to the problem structure
 * @param i Indices of the nodes
 * @param map Array of size nLocalNode to store the local indices of the nodes
 * @param x Array of size nLocalNode to store the x-coordinates of the nodes
 * @param y  Array of size nLocalNode to store the y-coordinates of the nodes
 */
void                 femPoissonLocal(femPoissonProblem *theProblem, const int i, int *map, double *x, double *y);

/**
 * @brief Free all the memory allocated for the structures in the problem
 * @param theProblem Pointer to the problem structure
 */
void                 femPoissonFree(femPoissonProblem *theProblem);

/**
 * @brief Assemble the stiffness matrix and the right-hand side of the system of equations then constrain the nodes on the boundary before solving the system
 * 
 * @param theProblem Pointer to the problem structure
 */
void                 femPoissonSolve(femPoissonProblem *theProblem);
#pragma endregion

#pragma region Utils
/**
 * @brief Compute the minimum of an array
 * @param x Array of size n
 * @param n Size of the array
 * @return Minimum value of the array
 */
double              femMin(double *x, int n);

/**
 * @brief Compute the maximum of an array
 * @param x Array of size n
 * @param n Size of the array
 * @return Maximum value of the array
 */
double              femMax(double *x, int n);

/**
 * @brief Write formatted global error message to the standard output
 * @param text Text to write
 * @param line Line number
 * @param file File name
 */
void                femError(char *text, int line, char *file);

/**
 * @brief Write formatted file reading error message to the standard output and exit the program
 * @param test Test to perform
 * @param line Line number
 * @param file File name
 */
void                femErrorScan(int test, int line, char *file);

/**
 * @brief Write formatted GMSH warning message to the standard output
 * @param text Text to write
 * @param line Line number
 * @param file File name
 */
void                femErrorGmsh(int test, int line, char *file);

/**
 * @brief Write formatted global warning message to the standard output
 * @param text Text to write
 * @param line Line number
 * @param file File name
 */
void                femWarning(char *text, int line, char *file);
#pragma endregion

#endif
