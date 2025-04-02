/*
 *  main.c
 *  Library for EPL1110 : Finite Elements for dummies
 *  Utilisation de l'API de GMSH pour créer un maillage
 *
 *  Copyright (C) 2023 UCL-IMMC : Vincent Legat
 *  All rights reserved.
 *
 */
 
#include "glfem.h"


int main(void)
{  
    printf("\n\n    V : Mesh and size mesh field \n");
    printf("    D : Domains \n");
    printf("    N : Next domain highlighted\n");



 
    double w = 1.0;
    double h = 2.0;
      
    int ierr;
    
    geoInitialize();
    femGeo* theGeometry = geoGetGeometry();
    

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

   
    geoMeshGenerate();
    geoMeshImport();
    
    

//
//  -2- Creation du fichier du maillage
//
    
    char filename[] = "../data/mesh.txt";
    geoMeshWrite(filename);

//
//  -3- Champ de la taille de référence du maillage
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

 
