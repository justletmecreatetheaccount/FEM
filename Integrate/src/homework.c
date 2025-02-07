#include <stdio.h>
#include <math.h>
#include "glfem.h"


double integrate(double x[3], double y[3], double (*f) (double, double))
{
    double I = 0;
    double xLoc[3];
    double yLoc[3];

    // Utilities for the integration
    double eta[3] = {1.0 / 6.0, 2.0 / 3.0, 1.0 / 6.0};
    double xi[3] = {1.0 / 6.0, 1.0 / 6.0, 2.0 / 3.0};

    // Compute the points of integration
    for (int i = 0; i < 3; i++)
    {
        xLoc[i] = x[0] * (1 - xi[i] - eta[i]) + x[2] * xi[i] + x[3] * eta[i];
        yLoc[i] = y[0] * (1 - xi[i] - eta[i]) + y[2] * xi[i] + y[3] * eta[i];
    }

    // Compute the integral
    for (int i = 0; i < 3; i++)
    {
        I += f(xLoc[i], yLoc[i]) / 6.0;
    }

  glfemSetColor(GLFEM_BLACK); glfemDrawElement(x,y,3);
  glfemSetColor(GLFEM_BLUE);  glfemDrawNodes(x,y,3);
  glfemSetColor(GLFEM_RED);   glfemDrawNodes(xLoc,yLoc,3);
    


    return I;
}

double integrateRecursive(double x[3], double y[3], double (*f)(double,double), int n)
{

//
// ... A modifier :-)
// y-compris la ligne juste en dessous :-)
//
    double I = integrate(x,y,f);
    
//
//
//    
     
    return I;
}
