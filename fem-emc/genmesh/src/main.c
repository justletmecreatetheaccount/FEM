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


int main(int argc, char *argv[])
{   char tri;
    char plot;
    if (argc < 2)
    {
        printf("\033[34m[INFO]You can to specifiy more arguments when creating the mesh, the first should 1 or 0 to create triangles or quads and the second 1 or 0 to plot it or not\033[0m\n");
        tri = 0;
    } else {
        tri = *argv[1];
        if (argc == 3){
            plot = *argv[2];
        } else {
            plot = 0;
        }
    }


    
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
    
    theGeometry->d = 0.5;// transition zone for hole
    theGeometry->s = 0.05;// refined mesh size near hole 1

    theGeometry->h = w * 0.1;

    if (tri == '1'){
        theGeometry->elementType = FEM_TRIANGLE;
    }
    if (tri == '0'){
        theGeometry->elementType = FEM_QUAD;
    }

  
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
    
        




    double *meshSizeField = malloc(theGeometry->theNodes->nNodes*sizeof(double));
    femNodes *theNodes = theGeometry->theNodes;
    for(int i=0; i < theNodes->nNodes; ++i)
        meshSizeField[i] = geoSize(theNodes->X[i], theNodes->Y[i]);
    double hMin = femMin(meshSizeField,theNodes->nNodes);  
    double hMax = femMax(meshSizeField,theNodes->nNodes);  
    printf(" ==== Global requested h : %14.7e \n",theGeometry->h);
    printf(" ==== Minimum h          : %14.7e \n",hMin);
    printf(" ==== Maximum h          : %14.7e \n",hMax);
    printf("\033[32m[SUCCESS] Le maillage a bien été générer\033[0m\n");
    
    if (plot == '1'){
        char theMessage[MAXNAME];
    
    
        GLFWwindow* window = glfemInit("Maillage Gr 29");
        glfwMakeContextCurrent(window);

        do {
            int w,h;
            glfwGetFramebufferSize(window,&w,&h);
            glfemReshapeWindows(theGeometry->theNodes,w,h);


            glfemPlotField(theGeometry->theElements,meshSizeField);
            glfemPlotMesh(theGeometry->theElements); 
            sprintf(theMessage, "Number of elements : %d ",theGeometry->theElements->nElem);
            glColor3f(1.0,0.0,0.0); glfemMessage(theMessage);

                
            glfwSwapBuffers(window);
            glfwPollEvents();
        } while( glfwGetKey(window,GLFW_KEY_ESCAPE) != GLFW_PRESS &&
                glfwWindowShouldClose(window) != 1 );

        glfwTerminate();    
    }
    
    free(meshSizeField);  
    geoFinalize();

 
    
    exit(EXIT_SUCCESS);
    return 0;  
}

 
