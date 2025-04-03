/*
 *  fem.c
 *  Library for LEPL1110 : Finite Elements for dummies
 *
 *  Copyright (C) 2021 UCL-IMMC : Vincent Legat
 *  All rights reserved.
 *
 */

#include "fem.h"

femGeo theGeometry;

femGeo *geoGetGeometry() { return &theGeometry; }

double geoSizeDefault(double x, double y) { return theGeometry.h; }

double geoGmshSize(int dim, int tag, double x, double y, double z, double lc, void *data)
{
    return theGeometry.geoSize(x, y);
}
void geoInitialize()
{
    int ierr;
    theGeometry.geoSize = geoSizeDefault;
    gmshInitialize(0, NULL, 1, 0, &ierr);
    ErrorGmsh(ierr);
    gmshModelAdd("MyGeometry", &ierr);
    ErrorGmsh(ierr);
    gmshModelMeshSetSizeCallback(geoGmshSize, NULL, &ierr);
    ErrorGmsh(ierr);
    theGeometry.theNodes = NULL;
    theGeometry.theElements = NULL;
    theGeometry.theEdges = NULL;
    theGeometry.nDomains = 0;
    theGeometry.theDomains = NULL;
}

void geoFinalize()
{
    int ierr;

    if (theGeometry.theNodes)
    {
        free(theGeometry.theNodes->X);
        free(theGeometry.theNodes->Y);
        free(theGeometry.theNodes);
    }
    if (theGeometry.theElements)
    {
        free(theGeometry.theElements->elem);
        free(theGeometry.theElements);
    }
    if (theGeometry.theEdges)
    {
        free(theGeometry.theEdges->elem);
        free(theGeometry.theEdges);
    }
    for (int i = 0; i < theGeometry.nDomains; i++)
    {
        free(theGeometry.theDomains[i]->elem);
        free(theGeometry.theDomains[i]);
    }
    free(theGeometry.theDomains);
    gmshFinalize(&ierr);
    ErrorGmsh(ierr);
}

void geoSetSizeCallback(double (*geoSize)(double x, double y))
{
    theGeometry.geoSize = geoSize;
}

void geoMeshImport()
{
    int ierr;

    /* Importing nodes */

    size_t nNode, n, m, *node;
    double *xyz, *trash;
    gmshModelMeshGetNodes(&node, &nNode, &xyz, &n,
                          &trash, &m, -1, -1, 0, 0, &ierr);
    ErrorGmsh(ierr);
    femNodes *theNodes = malloc(sizeof(femNodes));
    theNodes->nNodes = nNode;
    theNodes->X = malloc(sizeof(double) * (theNodes->nNodes));
    theNodes->Y = malloc(sizeof(double) * (theNodes->nNodes));
    for (int i = 0; i < theNodes->nNodes; i++)
    {
        theNodes->X[i] = xyz[3 * node[i] - 3];
        theNodes->Y[i] = xyz[3 * node[i] - 2];
    }
    theGeometry.theNodes = theNodes;
    gmshFree(node);
    gmshFree(xyz);
    gmshFree(trash);


    /* Importing elements */
    /* Pas super joli : a ameliorer pour eviter la triple copie */

    size_t nElem, *elem;
    gmshModelMeshGetElementsByType(1, &elem, &nElem,
                                   &node, &nNode, -1, 0, 1, &ierr);
    ErrorGmsh(ierr);
    femMesh *theEdges = malloc(sizeof(femMesh));
    theEdges->nLocalNode = 2;
    theEdges->nodes = theNodes;
    theEdges->nElem = nElem;
    theEdges->elem = malloc(sizeof(int) * 2 * theEdges->nElem);
    for (int i = 0; i < theEdges->nElem; i++)
        for (int j = 0; j < theEdges->nLocalNode; j++)
            theEdges->elem[2 * i + j] = node[2 * i + j] - 1;
    theGeometry.theEdges = theEdges;
    int shiftEdges = elem[0];
    gmshFree(node);
    gmshFree(elem);

    gmshModelMeshGetElementsByType(2, &elem, &nElem,
                                   &node, &nNode, -1, 0, 1, &ierr);
    ErrorGmsh(ierr);
    if (nElem != 0)
    {
        femMesh *theElements = malloc(sizeof(femMesh));
        theElements->nLocalNode = 3;
        theElements->nodes = theNodes;
        theElements->nElem = nElem;
        theElements->elem = malloc(sizeof(int) * 3 * theElements->nElem);
        for (int i = 0; i < theElements->nElem; i++)
            for (int j = 0; j < theElements->nLocalNode; j++)
                theElements->elem[3 * i + j] = node[3 * i + j] - 1;
        theGeometry.theElements = theElements;
        gmshFree(node);
        gmshFree(elem);
    }

    int nElemTriangles = nElem;
    gmshModelMeshGetElementsByType(3, &elem, &nElem,
                                   &node, &nNode, -1, 0, 1, &ierr);
    ErrorGmsh(ierr);
    if (nElem != 0 && nElemTriangles != 0)
        Error("Cannot consider hybrid geometry with triangles and quads :-(");

    if (nElem != 0)
    {
        femMesh *theElements = malloc(sizeof(femMesh));
        theElements->nLocalNode = 4;
        theElements->nodes = theNodes;
        theElements->nElem = nElem;
        theElements->elem = malloc(sizeof(int) * 4 * theElements->nElem);
        for (int i = 0; i < theElements->nElem; i++)
            for (int j = 0; j < theElements->nLocalNode; j++)
                theElements->elem[4 * i + j] = node[4 * i + j] - 1;
        theGeometry.theElements = theElements;
        gmshFree(node);
        gmshFree(elem);
    }

    /* Importing 1D entities */

    int *dimTags;
    gmshModelGetEntities(&dimTags, &n, 1, &ierr);
    ErrorGmsh(ierr);
    theGeometry.nDomains = n / 2;
    theGeometry.theDomains = malloc(sizeof(femDomain *) * n / 2);

    for (int i = 0; i < n / 2; i++)
    {
        int dim = dimTags[2 * i + 0];
        int tag = dimTags[2 * i + 1];
        femDomain *theDomain = malloc(sizeof(femDomain));
        theGeometry.theDomains[i] = theDomain;
        theDomain->mesh = theEdges;
        sprintf(theDomain->name, "Entity %d ", tag - 1);

        int *elementType;
        size_t nElementType, **elementTags, *nElementTags, nnElementTags, **nodesTags, *nNodesTags, nnNodesTags;
        gmshModelMeshGetElements(&elementType, &nElementType, &elementTags, &nElementTags, &nnElementTags, &nodesTags, &nNodesTags, &nnNodesTags, dim, tag, &ierr);
        theDomain->nElem = nElementTags[0];
        theDomain->elem = malloc(sizeof(int) * 2 * theDomain->nElem);
        for (int j = 0; j < theDomain->nElem; j++)
        {
            theDomain->elem[j] = elementTags[0][j] - shiftEdges;
        }
        gmshFree(nElementTags);
        gmshFree(nNodesTags);
        gmshFree(elementTags);
        gmshFree(nodesTags);
        gmshFree(elementType);
    }
    gmshFree(dimTags);

    return;
}

void geoMeshPrint()
{
    femNodes *theNodes = theGeometry.theNodes;
    if (theNodes != NULL)
    {
        printf("Number of nodes %d \n", theNodes->nNodes);
        for (int i = 0; i < theNodes->nNodes; i++)
        {
            printf("%6d : %14.7e %14.7e \n", i, theNodes->X[i], theNodes->Y[i]);
        }
    }
    femMesh *theEdges = theGeometry.theEdges;
    if (theEdges != NULL)
    {
        printf("Number of edges %d \n", theEdges->nElem);
        int *elem = theEdges->elem;
        for (int i = 0; i < theEdges->nElem; i++)
        {
            printf("%6d : %6d %6d \n", i, elem[2 * i], elem[2 * i + 1]);
        }
    }
    femMesh *theElements = theGeometry.theElements;
    if (theElements != NULL)
    {
        if (theElements->nLocalNode == 3)
        {
            printf("Number of triangles %d \n", theElements->nElem);
            int *elem = theElements->elem;
            for (int i = 0; i < theElements->nElem; i++)
            {
                printf("%6d : %6d %6d %6d\n", i, elem[3 * i], elem[3 * i + 1], elem[3 * i + 2]);
            }
        }
        if (theElements->nLocalNode == 4)
        {
            printf("Number of quads %d \n", theElements->nElem);
            int *elem = theElements->elem;
            for (int i = 0; i < theElements->nElem; i++)
            {
                printf("%6d : %6d %6d %6d %6d\n", i, elem[4 * i], elem[4 * i + 1], elem[4 * i + 2], elem[4 * i + 3]);
            }
        }
    }
    int nDomains = theGeometry.nDomains;
    printf("Number of domains %d\n", nDomains);
    for (int iDomain = 0; iDomain < nDomains; iDomain++)
    {
        femDomain *theDomain = theGeometry.theDomains[iDomain];
        printf("  Domain : %6d \n", iDomain);
        printf("  Name : %s\n", theDomain->name);
        printf("  Number of elements : %6d\n", theDomain->nElem);
        for (int i = 0; i < theDomain->nElem; i++)
        {
            //         if (i != theDomain->nElem  && (i % 10) != 0)  printf(" - ");
            printf("%6d", theDomain->elem[i]);
            if ((i + 1) != theDomain->nElem && (i + 1) % 10 == 0)
                printf("\n");
        }
        printf("\n");
    }
}

void geoMeshWrite(const char *filename)
{
    FILE *file = fopen(filename, "w");

    femNodes *theNodes = theGeometry.theNodes;
    fprintf(file, "Number of nodes %d \n", theNodes->nNodes);
    for (int i = 0; i < theNodes->nNodes; i++)
    {
        fprintf(file, "%6d : %14.7e %14.7e \n", i, theNodes->X[i], theNodes->Y[i]);
    }

    femMesh *theEdges = theGeometry.theEdges;
    fprintf(file, "Number of edges %d \n", theEdges->nElem);
    int *elem = theEdges->elem;
    for (int i = 0; i < theEdges->nElem; i++)
    {
        fprintf(file, "%6d : %6d %6d \n", i, elem[2 * i], elem[2 * i + 1]);
    }

    femMesh *theElements = theGeometry.theElements;
    if (theElements->nLocalNode == 3)
    {
        fprintf(file, "Number of triangles %d \n", theElements->nElem);
        elem = theElements->elem;
        for (int i = 0; i < theElements->nElem; i++)
        {
            fprintf(file, "%6d : %6d %6d %6d\n", i, elem[3 * i], elem[3 * i + 1], elem[3 * i + 2]);
        }
    }
    if (theElements->nLocalNode == 4)
    {
        fprintf(file, "Number of quads %d \n", theElements->nElem);
        elem = theElements->elem;
        for (int i = 0; i < theElements->nElem; i++)
        {
            fprintf(file, "%6d : %6d %6d %6d %6d\n", i, elem[4 * i], elem[4 * i + 1], elem[4 * i + 2], elem[4 * i + 3]);
        }
    }

    int nDomains = theGeometry.nDomains;
    fprintf(file, "Number of domains %d\n", nDomains);
    for (int iDomain = 0; iDomain < nDomains; iDomain++)
    {
        femDomain *theDomain = theGeometry.theDomains[iDomain];
        fprintf(file, "  Domain : %6d \n", iDomain);
        fprintf(file, "  Name : %s\n", theDomain->name);
        fprintf(file, "  Number of elements : %6d\n", theDomain->nElem);
        for (int i = 0; i < theDomain->nElem; i++)
        {
            fprintf(file, "%6d", theDomain->elem[i]);
            if ((i + 1) != theDomain->nElem && (i + 1) % 10 == 0)
                fprintf(file, "\n");
        }
        fprintf(file, "\n");
    }

    fclose(file);
}

void geoMeshRead(const char *filename)
{
    FILE *file = fopen(filename, "r");

    int trash, *elem;

    femNodes *theNodes = malloc(sizeof(femNodes));
    theGeometry.theNodes = theNodes;
    ErrorScan(fscanf(file, "Number of nodes %d \n", &theNodes->nNodes));
    theNodes->X = malloc(sizeof(double) * (theNodes->nNodes));
    theNodes->Y = malloc(sizeof(double) * (theNodes->nNodes));
    for (int i = 0; i < theNodes->nNodes; i++)
    {
        ErrorScan(fscanf(file, "%d : %le %le \n", &trash, &theNodes->X[i], &theNodes->Y[i]));
    }

    femMesh *theEdges = malloc(sizeof(femMesh));
    theGeometry.theEdges = theEdges;
    theEdges->nLocalNode = 2;
    theEdges->nodes = theNodes;
    ErrorScan(fscanf(file, "Number of edges %d \n", &theEdges->nElem));
    theEdges->elem = malloc(sizeof(int) * theEdges->nLocalNode * theEdges->nElem);
    for (int i = 0; i < theEdges->nElem; ++i)
    {
        elem = theEdges->elem;
        ErrorScan(fscanf(file, "%6d : %6d %6d \n", &trash, &elem[2 * i], &elem[2 * i + 1]));
    }

    femMesh *theElements = malloc(sizeof(femMesh));
    theGeometry.theElements = theElements;
    theElements->nLocalNode = 0;
    theElements->nodes = theNodes;
    char elementType[MAXNAME];
    ErrorScan(fscanf(file, "Number of %s %d \n", elementType, &theElements->nElem));
    if (strncasecmp(elementType, "triangles", MAXNAME) == 0)
    {
        theElements->nLocalNode = 3;
        theElements->elem = malloc(sizeof(int) * theElements->nLocalNode * theElements->nElem);
        for (int i = 0; i < theElements->nElem; ++i)
        {
            elem = theElements->elem;
            ErrorScan(fscanf(file, "%6d : %6d %6d %6d \n",
                             &trash, &elem[3 * i], &elem[3 * i + 1], &elem[3 * i + 2]));
        }
    }
    if (strncasecmp(elementType, "quads", MAXNAME) == 0)
    {
        theElements->nLocalNode = 4;
        theElements->elem = malloc(sizeof(int) * theElements->nLocalNode * theElements->nElem);
        for (int i = 0; i < theElements->nElem; ++i)
        {
            elem = theElements->elem;
            ErrorScan(fscanf(file, "%6d : %6d %6d %6d %6d \n",
                             &trash, &elem[4 * i], &elem[4 * i + 1], &elem[4 * i + 2], &elem[4 * i + 3]));
        }
    }

    ErrorScan(fscanf(file, "Number of domains %d\n", &theGeometry.nDomains));
    int nDomains = theGeometry.nDomains;
    theGeometry.theDomains = malloc(sizeof(femDomain *) * nDomains);
    for (int iDomain = 0; iDomain < nDomains; iDomain++)
    {
        femDomain *theDomain = malloc(sizeof(femDomain));
        theGeometry.theDomains[iDomain] = theDomain;
        theDomain->mesh = theEdges;
        ErrorScan(fscanf(file, "  Domain : %6d \n", &trash));
        ErrorScan(fscanf(file, "  Name : %[^\n]s \n", (char *)&theDomain->name));
        ErrorScan(fscanf(file, "  Number of elements : %6d\n", &theDomain->nElem));
        theDomain->elem = malloc(sizeof(int) * 2 * theDomain->nElem);
        for (int i = 0; i < theDomain->nElem; i++)
        {
            ErrorScan(fscanf(file, "%6d", &theDomain->elem[i]));
            if ((i + 1) != theDomain->nElem && (i + 1) % 10 == 0)
                ErrorScan(fscanf(file, "\n"));
        }
    }

    fclose(file);
}

void geoSetDomainName(int iDomain, char *name)
{
    if (iDomain >= theGeometry.nDomains)
        Error("Illegal domain number");
    if (geoGetDomain(name) != -1)
        Error("Cannot use the same name for two domains");
    sprintf(theGeometry.theDomains[iDomain]->name, "%s", name);
}

int geoGetDomain(char *name)
{
    int theIndex = -1;
    int nDomains = theGeometry.nDomains;
    for (int iDomain = 0; iDomain < nDomains; iDomain++)
    {
        femDomain *theDomain = theGeometry.theDomains[iDomain];
        if (strncasecmp(name, theDomain->name, MAXNAME) == 0)
            theIndex = iDomain;
    }
    return theIndex;
}

double geoSize(double x, double y)
{
    femGeo *theGeometry = geoGetGeometry();

    double h = theGeometry->h;

    // Hole 1 parameters
    double x0 = theGeometry->x_hole1 + 0.5 * theGeometry->w_hole1; // center X
    double y0 = theGeometry->y_hole1 + 0.5 * theGeometry->h_hole1; // center Y
    double halfW0 = 0.5 * theGeometry->w_hole1;
    double halfH0 = 0.5 * theGeometry->h_hole1;
    double d0 = theGeometry->d;
    double h0 = theGeometry->s;

    // Hole 2 parameters
    double x1 = theGeometry->x_hole2 + 0.5 * theGeometry->w_hole2; // center X
    double y1 = theGeometry->y_hole2 + 0.5 * theGeometry->h_hole2; // center Y
    double halfW1 = 0.5 * theGeometry->w_hole2;
    double halfH1 = 0.5 * theGeometry->h_hole2;
    double d1 = theGeometry->d;
    double h1 = theGeometry->s;

    // Hole 3 parameters
    double x2 = theGeometry->x_hole3 + 0.5 * theGeometry->w_hole3; // center X
    double y2 = theGeometry->y_hole3 + 0.5 * theGeometry->h_hole3; // center Y
    double halfW2 = 0.5 * theGeometry->w_hole3;
    double halfH2 = 0.5 * theGeometry->h_hole3;
    double d2 = theGeometry->d;
    double h2 = theGeometry->s;

    double result = h;

    // Distance to hole 1
    double dx0 = fabs(x - x0) - halfW0;
    double dy0 = fabs(y - y0) - halfH0;
    if (dx0 < 0.0) dx0 = 0.0;
    if (dy0 < 0.0) dy0 = 0.0;
    double dist0 = sqrt(dx0 * dx0 + dy0 * dy0);

    // Distance to hole 2
    double dx1 = fabs(x - x1) - halfW1;
    double dy1 = fabs(y - y1) - halfH1;
    if (dx1 < 0.0) dx1 = 0.0;
    if (dy1 < 0.0) dy1 = 0.0;
    double dist1 = sqrt(dx1 * dx1 + dy1 * dy1);

    // Distance to hole 3
    double dx2 = fabs(x - x2) - halfW2;
    double dy2 = fabs(y - y2) - halfH2;
    if (dx2 < 0.0) dx2 = 0.0;
    if (dy2 < 0.0) dy2 = 0.0;
    double dist2 = sqrt(dx2 * dx2 + dy2 * dy2);

    // Adjust mesh size near hole 1
    if (dist0 <= d0)
    {
        double a2 = (3.0 * (h - h0)) / (d0 * d0);
        double a3 = (2.0 * (h0 - h)) / (d0 * d0 * d0);
        double val0 = h0 + a2 * (dist0 * dist0) + a3 * (dist0 * dist0 * dist0);
        result = val0 < result ? val0 : result;
    }

    // Adjust mesh size near hole 2
    if (dist1 <= d1)
    {
        double a2 = (3.0 * (h - h1)) / (d1 * d1);
        double a3 = (2.0 * (h1 - h)) / (d1 * d1 * d1);
        double val1 = h1 + a2 * (dist1 * dist1) + a3 * (dist1 * dist1 * dist1);
        result = val1 < result ? val1 : result;
    }

    // Adjust mesh size near hole 3
    if (dist2 <= d2)
    {
        double a2 = (3.0 * (h - h2)) / (d2 * d2);
        double a3 = (2.0 * (h2 - h)) / (d2 * d2 * d2);
        double val2 = h2 + a2 * (dist2 * dist2) + a3 * (dist2 * dist2 * dist2);
        result = val2 < result ? val2 : result;
    }

    return result;
}
void geoMeshGenerate()
{
    femGeo *theGeometry = geoGetGeometry();
    geoSetSizeCallback(geoSize);

    int ierr;

    int idRect = gmshModelOccAddRectangle(theGeometry->x_plate, theGeometry->y_plate, 0.0, theGeometry->w_plate, theGeometry->h_plate, -1, 0.1, &ierr);

    int idHoleDown = gmshModelOccAddRectangle(theGeometry->x_hole1, theGeometry->y_hole1, 0.0, theGeometry->w_hole1, theGeometry->h_hole1, -1, 0.1, &ierr);
    int idHoleUp = gmshModelOccAddRectangle(theGeometry->x_hole2, theGeometry->y_hole2, 0.0, theGeometry->w_hole2, theGeometry->h_hole2, -1, 0.1, &ierr);

    // Third hole
    int idHoleThird = gmshModelOccAddRectangle(theGeometry->x_hole3, theGeometry->y_hole3, 0.0,
                                               theGeometry->w_hole3, theGeometry->h_hole3,
                                               -1, 0.1, &ierr);

    int rect[] = {2, idRect};
    int holeDown[] = {2, idHoleDown};
    int holeUp[] = {2, idHoleUp};
    int holeThird[] = {2, idHoleThird};

    gmshModelOccCut(rect, 2, holeDown, 2, NULL, NULL, NULL, NULL, NULL, -1, 1, 1, &ierr);
    ErrorGmsh(ierr);
    gmshModelOccCut(rect, 2, holeUp, 2, NULL, NULL, NULL, NULL, NULL, -1, 1, 1, &ierr);
    ErrorGmsh(ierr);
    gmshModelOccCut(rect, 2, holeThird, 2, NULL, NULL, NULL, NULL, NULL, -1, 1, 1, &ierr);
    ErrorGmsh(ierr);

    gmshModelOccSynchronize(&ierr);
    geoSetSizeCallback(geoSize);

    gmshModelOccSynchronize(&ierr);

    if (theGeometry->elementType == FEM_TRIANGLE)
    {
        gmshOptionSetNumber("Mesh.SaveAll", 1, &ierr);
        gmshModelMeshGenerate(2, &ierr);
    }
    else
    {
        gmshOptionSetNumber("Mesh.SaveAll", 1, &ierr);
        gmshOptionSetNumber("Mesh.RecombineAll", 1, &ierr);
        gmshOptionSetNumber("Mesh.Algorithm", 8, &ierr);
        gmshOptionSetNumber("Mesh.RecombinationAlgorithm", 1.0, &ierr);
        gmshModelGeoMeshSetRecombine(2, 1, 45, &ierr);
        gmshModelMeshGenerate(2, &ierr);
    }
    //gmshFltkInitialize(&ierr);
    //gmshFltkRun(&ierr);
}
double femMin(double *x, int n)
{
    double myMin = x[0];
    int i;
    for (i = 1; i < n; i++)
        myMin = fmin(myMin, x[i]);
    return myMin;
}

double femMax(double *x, int n)
{
    double myMax = x[0];
    int i;
    for (i = 1; i < n; i++)
        myMax = fmax(myMax, x[i]);
    return myMax;
}

void femError(char *text, int line, char *file)
{
    printf("\n-------------------------------------------------------------------------------- ");
    printf("\n  Error in %s at line %d : \n  %s\n", file, line, text);
    printf("--------------------------------------------------------------------- Yek Yek !! \n\n");
    exit(69);
}

void femErrorGmsh(int ierr, int line, char *file)
{
    if (ierr == 0)
        return;
    printf("\n-------------------------------------------------------------------------------- ");
    printf("\n  Error in %s at line %d : \n  error code returned by gmsh %d\n", file, line, ierr);
    printf("--------------------------------------------------------------------- Yek Yek !! \n\n");
    gmshFinalize(NULL);
    exit(69);
}

void femErrorScan(int test, int line, char *file)
{
    if (test >= 0)
        return;

    printf("\n-------------------------------------------------------------------------------- ");
    printf("\n  Error in fscanf or fgets in %s at line %d : \n", file, line);
    printf("--------------------------------------------------------------------- Yek Yek !! \n\n");
    exit(69);
}

void femWarning(char *text, int line, char *file)
{
    printf("\n-------------------------------------------------------------------------------- ");
    printf("\n  Warning in %s at line %d : \n  %s\n", file, line, text);
    printf("--------------------------------------------------------------------- Yek Yek !! \n\n");
}