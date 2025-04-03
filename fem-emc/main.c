#include "main.h"
#include <time.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "project/fem.h"
#include "genmesh/src/glfem.h"
#include <unistd.h>



int main(int argc, char** argv){
    // parsing args
    char* meshfile = "../data/mesh.txt";
    char* outfile  = "../data/UV.txt";
    if (argc > 1) meshfile = argv[1];
    if (argc > 2) outfile  = argv[2];
    // small benchmark
    struct timespec t0, t1;
    timespec_get(&t0, TIME_UTC);
    Plot plot = elasticity_solve(meshfile, outfile);
    timespec_get(&t1, TIME_UTC);
    double exec_time = (t1.tv_sec - t0.tv_sec)*1.0 + (t1.tv_nsec - t0.tv_nsec)*1e-9;
    printf("\033[34m[INFO]Your code runs in %.4fs for mesh file '%s'\033[0m\n",exec_time, meshfile);
    
    // extracting arrays for the plot
    femProblem* theProblem = plot.theProblem;
    femGeo* theGeometry = plot.theGeometry;
    double* theSoluce = plot.theSoluce;
    double* theForces = plot.theForces;

    // plot
    femNodes *theNodes = theGeometry->theNodes;
    double maxDeformation = 1e5;
    double deformationFactor = 0.0;
    double deformationStep = 1e4;
    double prev = 0;
    double *normDisplacement = malloc(theNodes->nNodes * sizeof(double));
    double *forcesX = malloc(theNodes->nNodes * sizeof(double));
    double *forcesY = malloc(theNodes->nNodes * sizeof(double));
    

    int mode = 1; 
    int domain = 0;
    int freezingButton = FALSE;
    double t, told = 0;
    char theMessage[MAXNAME];

 
    GLFWwindow* window = glfemInit("Solution Gr 29");
    glfwMakeContextCurrent(window);

    do {
        int w,h;
        glfwGetFramebufferSize(window,&w,&h);
        glfemReshapeWindows(theGeometry->theNodes,w,h);

        t = glfwGetTime();  
        if (glfwGetKey(window,'D') == GLFW_PRESS) { mode = 0;}
        if (glfwGetKey(window,'V') == GLFW_PRESS) { mode = 1;}
        if (glfwGetKey(window,'X') == GLFW_PRESS) { mode = 2;}
        if (glfwGetKey(window,'Y') == GLFW_PRESS) { mode = 3;}
        if (glfwGetKey(window,'N') == GLFW_PRESS && freezingButton == FALSE) { domain++; freezingButton = TRUE; told = t;}
        if (t-told > 0.5) {freezingButton = FALSE; }
        prev = deformationFactor;
        deformationFactor += deformationStep;
        if (deformationFactor >= maxDeformation){
            deformationFactor = maxDeformation;
        }
        int nNodes = theNodes->nNodes;
        for (int i = 0; i < nNodes; i++){
            theNodes->X[i] -= theSoluce[2*i+0]*prev;
            theNodes->Y[i] -= theSoluce[2*i+1]*prev;
        }
        for (int i = 0; i < nNodes; i++){
            theNodes->X[i] += theSoluce[2*i+0]*deformationFactor;
            theNodes->Y[i] += theSoluce[2*i+1]*deformationFactor;
            normDisplacement[i] = sqrt(theSoluce[2*i+0]*theSoluce[2*i+0] + 
                            theSoluce[2*i+1]*theSoluce[2*i+1]);
            forcesX[i] = theForces[2*i+0];
            forcesY[i] = theForces[2*i+1];
        }
        if (mode == 0) {
            domain = domain % theGeometry->nDomains;
            glfemPlotDomain( theGeometry->theDomains[domain]); 
            snprintf(theMessage, sizeof(theMessage), "%.400s : %d",
                    theGeometry->theDomains[domain]->name, domain);
            glColor3f(1.0,0.0,0.0); glfemMessage(theMessage); }
        if (mode == 1) {
            glfemPlotField(theGeometry->theElements,normDisplacement);
            glfemPlotMesh(theGeometry->theElements); 
            sprintf(theMessage, "Number of elements : %d ",theGeometry->theElements->nElem);
            glColor3f(1.0,0.0,0.0); glfemMessage(theMessage); }
        if (mode == 2) {
            glfemPlotField(theGeometry->theElements,forcesX);
            glfemPlotMesh(theGeometry->theElements); 
            sprintf(theMessage, "Number of elements : %d ",theGeometry->theElements->nElem);
            glColor3f(1.0,0.0,0.0); glfemMessage(theMessage); }
        if (mode == 3) {
            glfemPlotField(theGeometry->theElements,forcesY);
            glfemPlotMesh(theGeometry->theElements); 
            sprintf(theMessage, "Number of elements : %d ",theGeometry->theElements->nElem);
            glColor3f(1.0,0.0,0.0); glfemMessage(theMessage); }
         glfwSwapBuffers(window);
         glfwPollEvents();
        usleep(40000);
    } while( glfwGetKey(window,GLFW_KEY_ESCAPE) != GLFW_PRESS &&
             glfwWindowShouldClose(window) != 1 );
    free(normDisplacement);
    free(forcesX);
    free(forcesY);
    femElasticityFree(theProblem);
    geoFinalize(theGeometry);
    return 0;
}
