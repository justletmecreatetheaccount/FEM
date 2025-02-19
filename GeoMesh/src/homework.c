#include "fem.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <strings.h>

int sum_from_to(int from, int to) {
    if (from > to) {
        return (from + to) * (from - to + 1) / 2;
    } else {
        return (from + to) * (to - from + 1) / 2;
    }
}

void dividedDifferenceArray(double* result, double* x, double* y_and_derivatives, int number_of_points, int number_of_derivatives) {
    double* adjusted_x = malloc(number_of_points * (number_of_derivatives + 1) * sizeof(double));

    for (int i = 0; i < number_of_points; i++) {
        for (int j = 0; j < number_of_derivatives + 1; j++) {
            adjusted_x[i * (number_of_derivatives + 1) + j] = x[i];
        }
    }

    for (int i = 0; i < number_of_points; i++) {
        for (int j = 0; j < number_of_derivatives + 1; j++) {
            result[i * (number_of_derivatives + 1) + j] = y_and_derivatives[i];
        }
    }

    for (int l = 1; l < number_of_points * (number_of_derivatives + 1); l++) {
        int index = 0;
        int result_column_start_index = sum_from_to(number_of_points * (number_of_derivatives + 1), number_of_points * (number_of_derivatives + 1) - l + 1);
        int prev_result_column_start_index = result_column_start_index - number_of_points * (number_of_derivatives + 1) + l - 1;
        
        if (l < number_of_derivatives + 1) {
            for (int i = 0; i < number_of_points; i++) {
                for (int d = 0; d < number_of_derivatives + 1 - l; d++) {
                    result[result_column_start_index + index] = y_and_derivatives[i + l * number_of_points] / l;
                    index++;
                }
                if (index < number_of_points * (number_of_derivatives + 1) - l) {
                    for (int nd = 0; nd < l; nd++) {
                        result[result_column_start_index + index] = (result[prev_result_column_start_index + index + 1] - result[prev_result_column_start_index + index]) / (adjusted_x[index + l] - adjusted_x[index]);
                        index++;
                    }
                }
            }
        } else {
            for (int i = 0; i < number_of_points * (number_of_derivatives + 1) - l; i++) {
                result[result_column_start_index + index] = (result[prev_result_column_start_index + index + 1] - result[prev_result_column_start_index + index]) / (adjusted_x[index + l] - adjusted_x[index]);
                index++;
            }
        }
    }
    free(adjusted_x);
}

double computeHermitePolynomeValue(double* x, double* result, int number_of_points, int number_of_derivatives, double t) {

    double* adjusted_x = malloc(number_of_points * (number_of_derivatives + 1) * sizeof(double));

    for (int i = 0; i < number_of_points; i++) {
        for (int j = 0; j < number_of_derivatives + 1; j++) {
            adjusted_x[i * (number_of_derivatives + 1) + j] = x[i];
        }
    }

    double result_value = result[0];
    double last_product;
    for (int i = 1; i < number_of_points * (number_of_derivatives + 1); i++) {
        double product = 1;
        for (int j = 0; j < i; j++) {
            product *= (t - adjusted_x[j]);
        }
        result_value += result[sum_from_to(number_of_points * (number_of_derivatives + 1) - i + 1, number_of_points * (number_of_derivatives + 1))] * product;
    }
    if (result_value < 0) {
        printf("result_valuel = %f\n", result_value);
        printf("t = %f\n", t);
        printf("result[%d] = %f\n", 0, result[0]);
        printf("last_product = %f\n", last_product);
        for (int i = 1; i < number_of_points * (number_of_derivatives + 1); i++) {
            printf("result[%d] = %f\n", i, result[sum_from_to(number_of_points * (number_of_derivatives + 1) - i + 1, number_of_points * (number_of_derivatives + 1))]);
        }
        for (int i = 0; i < number_of_points * (number_of_derivatives + 1); i++) {
            printf("adjusted_x[%d] = %f\n", i, adjusted_x[i]);
        }
        for (int i = 0; i < number_of_points; i++) {
            printf("x[%d] = %f\n", i, x[i]);
        }
        return 0.01;
    }
    free(adjusted_x);
    return result_value;
}

double geoSize(double x, double y){

    femGeo* theGeometry = geoGetGeometry();
    
    double h = theGeometry->h;
    double x0 = theGeometry->xNotch;
    double y0 = theGeometry->yNotch;
    double r0 = theGeometry->rNotch;
    double h0 = theGeometry->hNotch;
    double d0 = theGeometry->dNotch;
  
    
    double x1 = theGeometry->xHole;
    double y1 = theGeometry->yHole;
    double r1 = theGeometry->rHole;
    double h1 = theGeometry->hHole;
    double d1 = theGeometry->dHole;


//
//     A modifier !
//     
// Your contribution starts here ....
//
    
    double holeDist = sqrt((x - x1)*(x - x1) + (y - y1)*(y - y1)) - r1;
    //printf("holeDist = %f\n", holeDist);
    double notchDist = sqrt((x - x0)*(x - x0) + (y - y0)*(y - y0)) - r0;
    //printf("notchDist = %f\n", notchDist);
    if (notchDist < holeDist) {
        if (notchDist > d0) {
            return h;
        }
        double x[] = {0, d0};
        //printf("x[0] = %f, x[1] = %f\n", x[0], x[1]);
        double y[] = {h0, h, 0, 0};
        //printf("y[0] = %f, y[1] = %f, y[2] = %f, y[3] = %f\n", y[0], y[1], y[2], y[3]);
        double* pol = malloc(sum_from_to(0, 2 * (1 + 1)) * sizeof(double) + 1);
        dividedDifferenceArray(pol, x, y, 2, 1);
        double result = computeHermitePolynomeValue(x, pol, 2, 1, notchDist);
        free(pol);
        return result;
    } else {
        if (holeDist > d1) {
            return h;
        }
        double x[] = {0, d1};
        double y[] = {h1, h, 0, 0};
        double* pol = malloc(sum_from_to(0, 2 * (1 + 1)) * sizeof(double) + 1);
        dividedDifferenceArray(pol, x, y, 2, 1);
        double result = computeHermitePolynomeValue(x, pol, 2, 1, holeDist);
        free(pol);
        return result;
    }
    
//   
// Your contribution ends here :-)
//

}


#define ___ 0

void geoMeshGenerate() {

    femGeo* theGeometry = geoGetGeometry();

    double w = theGeometry->LxPlate;
    double h = theGeometry->LyPlate;
     
    double x0 = theGeometry->xNotch;
    double y0 = theGeometry->yNotch;
    double r0 = theGeometry->rNotch;
    
    
    double x1 = theGeometry->xHole;
    double y1 = theGeometry->yHole;
    double r1 = theGeometry->rHole;
 
//
//  -1- Construction de la géométrie avec OpenCascade
//      On crée le rectangle
//      On crée les deux cercles
//      On soustrait les cercles du rectangle :-)
//
 
    int ierr;
    int idPlate = gmshModelOccAddRectangle(-w/2.0, -h/2.0, 0.0, w, h, -1, 0.0, &ierr);
    ErrorGmsh(ierr);
    int idNotch = gmshModelOccAddDisk(x0, y0, 0.0, r0, r0, -1, NULL, 0, NULL, 0, &ierr); 
    ErrorGmsh(ierr);
    int idHole  = gmshModelOccAddDisk(x1, y1, 0.0, r1, r1, -1, NULL, 0, NULL, 0, &ierr);    
    ErrorGmsh(ierr);
    
    int plate[] = {2, idPlate};
    int notch[] = {2, idNotch};
    int hole[]  = {2, idHole};
    
    gmshModelOccCut(plate, 2, notch, 2, NULL, NULL, NULL, NULL, NULL, -1, 1, 1, &ierr); 
    ErrorGmsh(ierr);
    gmshModelOccCut(plate, 2, hole, 2, NULL, NULL, NULL, NULL, NULL, -1, 1, 1, &ierr); 
    ErrorGmsh(ierr);
 
//
//  -2- Définition de la fonction callback pour la taille de référence
//      Synchronisation de OpenCascade avec gmsh
//      Génération du maillage (avec l'option Mesh.SaveAll :-)
                  
   
    geoSetSizeCallback(geoSize);
                                  
    gmshModelOccSynchronize(&ierr);       
    gmshOptionSetNumber("Mesh.SaveAll", 1, &ierr);
    gmshModelMeshGenerate(2, &ierr);  
       
//
//  Generation de quads :-)
//
//    gmshOptionSetNumber("Mesh.SaveAll", 1, &ierr);
//    gmshOptionSetNumber("Mesh.RecombineAll", 1, &ierr);
//    gmshOptionSetNumber("Mesh.Algorithm", 8, &ierr);  chk(ierr);
//    gmshOptionSetNumber("Mesh.RecombinationAlgorithm", 1.0, &ierr);  chk(ierr);
//    gmshModelGeoMeshSetRecombine(2,1,45,&ierr);  chk(ierr);
//    gmshModelMeshGenerate(2, &ierr);  
   
 
//
//  Plot of Fltk
//
//   gmshFltkInitialize(&ierr);
//   gmshFltkRun(&ierr);  chk(ierr);
//
    
}