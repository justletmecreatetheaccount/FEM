/*
 *  main.c
 *  Library for EPL1110 : Finite Elements for dummies
 *  Elasticite lineaire plane
 *  Calcul des densités de force aux noeuds contraints
 *
 *  Copyright (C) 2024 UCL-IMMC : Vincent Legat
 *  All rights reserved.
 *
 */
 
#include "glfem.h"


int main(void)
{  
    geoInitialize();
    femGeo* theGeometry = geoGetGeometry();
    theGeometry->elementType = FEM_QUAD;

    double w = 1.0;
    double h = 2.0;
    theGeometry->w_plate = w;
    theGeometry->h_plate = h;
    theGeometry->x_plate = -w / 2.0;
    theGeometry->y_plate = -h / 2.0;
    
    theGeometry->x_hole1 = -w / 4.0;
    theGeometry->y_hole1 = -3 * h / 8.0;
    theGeometry->w_hole1 = w / 2.0;
    theGeometry->h_hole1 = h / 4.0;
    
    theGeometry->x_hole2 = -w / 4.0;
    theGeometry->y_hole2 = h / 8.0;
    theGeometry->w_hole2 = w / 2.0;
    theGeometry->h_hole2 = h / 4.0;
    
    theGeometry->h = w * 0.1;

    theGeometry->elementType = FEM_QUAD;
  
    geoMeshGenerate();
    geoMeshImport();


    geoSetDomainName(0, "Holeup bottom");
    geoSetDomainName(1, "Holeup bottom-right");
    geoSetDomainName(2, "Holeup right");
    geoSetDomainName(3, "Holeup top-right");
    geoSetDomainName(4, "Holeup top");
    geoSetDomainName(5, "Holeup top-left");
    geoSetDomainName(6, "Holeup left");
    geoSetDomainName(7, "Holeup bottom-left");
    geoSetDomainName(8, "Plate bottom");
    geoSetDomainName(9, "Plate bottom-left");
    geoSetDomainName(10, "Plate bottom-right");
    geoSetDomainName(11, "Plate left");
    geoSetDomainName(12, "Plate right");
    geoSetDomainName(13, "Plate top-left");
    geoSetDomainName(14, "Plate top-right");
    geoSetDomainName(15, "Plate top");
    geoSetDomainName(16, "Holedown bottom");
    geoSetDomainName(17, "Holedown bottom-left");
    geoSetDomainName(18, "Holedown bottom-right");
    geoSetDomainName(19, "Holedown left");
    geoSetDomainName(20, "Holedown right");
    geoSetDomainName(21, "Holedown top-left");
    geoSetDomainName(22, "Holedown top-right");
    geoSetDomainName(23, "Holedown top");


    geoMeshWrite("../data/mesh.txt");
    
        

//
//  -6- Visualisation du maillage
//  
    
    double *meshSizeField = malloc(theGeometry->theNodes->nNodes*sizeof(double));
    femNodes *theNodes = theGeometry->theNodes;
    for(int i=0; i < theNodes->nNodes; ++i)
        meshSizeField[i] = geoSize(theNodes->X[i], theNodes->Y[i]);
    double hMin = femMin(meshSizeField,theNodes->nNodes);  
    double hMax = femMax(meshSizeField,theNodes->nNodes);  
    printf(" ==== Global requested h : %14.7e \n",theGeometry->h);
    printf(" ==== Minimum h          : %14.7e \n",hMin);
    printf(" ==== Maximum h          : %14.7e \n",hMax);
 
    

    free(meshSizeField);  
    geoFinalize();
    glfwTerminate(); 
 
    
    exit(EXIT_SUCCESS);
    return 0;  
}

 
