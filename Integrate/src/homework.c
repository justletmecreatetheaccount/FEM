#include <stdio.h>
#include <math.h>
#include "glfem.h"


double integrate(double x[3], double y[3], double (*f) (double, double))
{
    double I = 0;
    double xLoc[3];
    double yLoc[3];

    // Utilities for the integration
    double xi[3] = {1.0 / 6.0, 1.0 / 6.0, 2.0 / 3.0};
    double eta[3] = {1.0 / 6.0, 2.0 / 3.0, 1.0 / 6.0};

    // Compute the points of integration
    for (int i = 0; i < 3; i++)
    {
        xLoc[i] = x[0] * (1 - xi[i] - eta[i]) + x[1] * xi[i] + x[2] * eta[i];
        yLoc[i] = y[0] * (1 - xi[i] - eta[i]) + y[1] * xi[i] + y[2] * eta[i];
    }

    // Compute the integral
    for (int i = 0; i < 3; i++)
    {
        I += f(xLoc[i], yLoc[i]) / 6.0;
    }
glfemSetColor(GLFEM_RED);   glfemDrawNodes(xLoc,yLoc,3);
  glfemSetColor(GLFEM_BLACK); glfemDrawElement(x,y,3);
  glfemSetColor(GLFEM_BLUE);  glfemDrawNodes(x,y,3);
  
    


    return I;
}

double integrateRecursive(double x[3], double y[3], double (*f)(double,double), int n)
{
    double I;

    if (n == 0)
    {
        I = integrate(x, y, f);
    }
    else
    {
        // Compute the mid points (x axis then y axis then hypotenuse)
        double mid_points[3][2];

        double xi[3] = {0.5, 0, 0.5};
        double eta[3] = {0, 0.5, 0.5};

        // Compute the points of integration
        for (int i = 0; i < 3; i++)
        {
            mid_points[i][0] = x[0] * (1 - xi[i] - eta[i]) + x[1] * xi[i] + x[2] * eta[i];
            mid_points[i][1] = y[0] * (1 - xi[i] - eta[i]) + y[1] * xi[i] + y[2] * eta[i];
        }

        // Compute the integral points
        double x0[3] = {mid_points[1][0], mid_points[2][0], x[2]}; //top triangle
        double y0[3] = {mid_points[1][1], mid_points[2][1], y[2]}; 
        double x1[3] = {x[0], mid_points[0][0], mid_points[1][0]}; //corner triangle
        double y1[3] = {y[0], mid_points[0][1], mid_points[1][1]};
        double x2[3] = {mid_points[0][0], x[1], mid_points[2][0]}; //right triangle
        double y2[3] = {mid_points[0][1], y[1], mid_points[2][1]};
        double x3[3] = {mid_points[2][0], mid_points[1][0], mid_points[0][0]}; //center triangle
        double y3[3] = {mid_points[2][1], mid_points[1][1], mid_points[0][1]};

        I = (integrateRecursive(x0, y0, f, n - 1) + integrateRecursive(x1, y1, f, n - 1) + integrateRecursive(x2, y2, f, n - 1) + integrateRecursive(x3, y3, f, n - 1));
        I = I / 4.0;
    }  
     
    return I;
}

