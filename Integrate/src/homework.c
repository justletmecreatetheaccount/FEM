#include <stdio.h>
#include <math.h>
#include "glfem.h"


double integrate(double x[3], double y[3], double (*f) (double, double))
{
    double I = 0;
    double xLoc[3];
    double yLoc[3];


    double xi[3] = {1.0/6.0, 1.0/6.0, 2.0/3.0};
    double eta[3] = {1.0/6.0, 2.0/3.0, 1.0/6.0};


    double jac =  (x[0] - x[1]) * (y[0] - y[2]) - (x[0] - x[2]) * (y[0] - y[1]);
    for (int i = 0; i < 3; i++) {
        xLoc[i] = (1 - xi[i] - eta[i]) * x[0] + xi[i] * x[1] + eta[i] * x[2];
        yLoc[i] = (1 - xi[i] - eta[i]) * y[0] + xi[i] * y[1] + eta[i] * y[2];
    }
    for (int i = 0; i < 3; i++) {
        I += f(xLoc[i], yLoc[i]) * (1.0/6.0);
    }

    glfemSetColor(GLFEM_BLACK); glfemDrawElement(x,y,3);
    glfemSetColor(GLFEM_BLUE);  glfemDrawNodes(x,y,3);
    glfemSetColor(GLFEM_RED);   glfemDrawNodes(xLoc,yLoc,3);
    


    return fabs(jac) * I;
}

double integrateRecursive(double x[3], double y[3], double (*f)(double,double), int n)
{
    double I = 0.0;
    if (n == 0) {
        I = integrate(x, y, f);
    } else {
        double xMid[3], yMid[3];
        for (int i = 0; i < 3; i++) {
            xMid[i] = (x[i] + x[(i + 1) % 3]) / 2.0; // i = [0,1,2]; (i+1) % 3 = [1,2,0]
            yMid[i] = (y[i] + y[(i + 1) % 3]) / 2.0;
        }

        double xSub[4][3] = {
            {x[0], xMid[0], xMid[2]},
            {xMid[0], x[1], xMid[1]},
            {xMid[2], xMid[1], x[2]},
            {xMid[0], xMid[1], xMid[2]}
        };

        double ySub[4][3] = {
            {y[0], yMid[0], yMid[2]},
            {yMid[0], y[1], yMid[1]},
            {yMid[2], yMid[1], y[2]},
            {yMid[0], yMid[1], yMid[2]}
        };

        for (int i = 0; i < 4; i++) {
            I += integrateRecursive(xSub[i], ySub[i], f, n - 1);
        }
    }
    return I;
}
