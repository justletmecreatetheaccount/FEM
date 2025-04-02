#include "fem.h"


/*
#define chk(ierr)                                               \
  if(ierr != 0){                                                \
    fprintf(stderr, "Error on line %i in function '%s': "       \
            "gmsh function returned non-zero error code: %i\n", \
            __LINE__, __FUNCTION__, ierr);                      \
    gmshFinalize(NULL);                                         \
    exit(0);} // ierr  --- for ctest                             
*/


double geoSize(double x, double y)
{
    femGeo *theGeometry = geoGetGeometry();

    double h = theGeometry->h;

    double x0 = theGeometry->x_hole1 + 0.5 * theGeometry->w_hole1;  // center X
    double y0 = theGeometry->y_hole1 + 0.5 * theGeometry->h_hole1;  // center Y
    double halfW0 = 0.5 * theGeometry->w_hole1;
    double halfH0 = 0.5 * theGeometry->h_hole1;
    double d0 = 0.5;      // transition zone for hole 1
    double h0 = 0.005;    // refined mesh size near hole 1

    double x1 = theGeometry->x_hole2 + 0.5 * theGeometry->w_hole2;  // center X
    double y1 = theGeometry->y_hole2 + 0.5 * theGeometry->h_hole2;  // center Y
    double halfW1 = 0.5 * theGeometry->w_hole2;
    double halfH1 = 0.5 * theGeometry->h_hole2;
    double d1 = 0.5;      // transition zone for hole 2
    double h1 = 0.01;     // refined mesh size near hole 2

    // Default mesh size
    double result = h;

    double dx0 = fabs(x - x0) - halfW0;
    double dy0 = fabs(y - y0) - halfH0;
    if (dx0 < 0.0) dx0 = 0.0;
    if (dy0 < 0.0) dy0 = 0.0;
    double dist0 = sqrt(dx0 * dx0 + dy0 * dy0);

    double dx1 = fabs(x - x1) - halfW1;
    double dy1 = fabs(y - y1) - halfH1;
    if (dx1 < 0.0) dx1 = 0.0;
    if (dy1 < 0.0) dy1 = 0.0;
    double dist1 = sqrt(dx1 * dx1 + dy1 * dy1);

    if (dist0 <= d0) {
        double a2 = (3.0 * (h - h0)) / (d0 * d0);
        double a3 = (2.0 * (h0 - h)) / (d0 * d0 * d0);
        double val0 = h0 + a2 * (dist0 * dist0) + a3 * (dist0 * dist0 * dist0);
        result = val0 < result ? val0 : result;
    }

    if (dist1 <= d1) {
        double a2 = (3.0 * (h - h1)) / (d1 * d1);
        double a3 = (2.0 * (h1 - h)) / (d1 * d1 * d1);
        double val1 = h1 + a2 * (dist1 * dist1) + a3 * (dist1 * dist1 * dist1);
        result = val1 < result ? val1 : result;
    }

    return result;
}



void geoMeshGenerate() {
  femGeo *theGeometry = geoGetGeometry();


  geoSetSizeCallback(geoSize);


  int ierr;

    int idRect = gmshModelOccAddRectangle(theGeometry->x_plate, theGeometry->y_plate, 0.0, theGeometry->w_plate, theGeometry->h_plate, -1, 0.1, &ierr);

    int idHoleDown = gmshModelOccAddRectangle(theGeometry->x_hole1, theGeometry->y_hole1, 0.0, theGeometry->w_hole1, theGeometry->h_hole1, -1, 0.1, &ierr);

    int idHoleUp = gmshModelOccAddRectangle(theGeometry->x_hole2, theGeometry->y_hole2, 0.0, theGeometry->w_hole2, theGeometry->h_hole2, -1, 0.1, &ierr);


  int rect[] = {2, idRect};
  int holeDown[] = {2, idHoleDown};
  int holeUp[] = {2, idHoleUp};
  gmshModelOccCut(rect, 2, holeDown, 2, NULL, NULL, NULL, NULL, NULL, -1, 1, 1, &ierr); 
  ErrorGmsh(ierr);
  gmshModelOccCut(rect, 2, holeUp, 2, NULL, NULL, NULL, NULL, NULL, -1, 1, 1, &ierr); 
  ErrorGmsh(ierr);
  gmshModelOccSynchronize(&ierr);
 
//
//  -2- D�finition de la fonction callback pour la taille de r�f�rence
//      Synchronisation de OpenCascade avec gmsh
//      G�n�ration du maillage (avec l'option Mesh.SaveAll :-)
                  
   
    geoSetSizeCallback(geoSize);
                                  
    gmshModelOccSynchronize(&ierr);       
//    gmshOptionSetNumber("Mesh.SaveAll", 1, &ierr);
//    gmshModelMeshGenerate(2, &ierr);  
       
//
//  Generation de quads :-)
//
    gmshOptionSetNumber("Mesh.SaveAll", 1, &ierr);
    gmshOptionSetNumber("Mesh.RecombineAll", 1, &ierr);
    gmshOptionSetNumber("Mesh.Algorithm", 8, &ierr);  //chk(ierr);
    gmshOptionSetNumber("Mesh.RecombinationAlgorithm", 1.0, &ierr);  //chk(ierr);
    gmshModelGeoMeshSetRecombine(2,1,45,&ierr);  //chk(ierr);
    gmshModelMeshGenerate(2, &ierr);  
   
 
//
//  Plot of Fltk
//
   gmshFltkInitialize(&ierr);
   gmshFltkRun(&ierr);  //chk(ierr);
//
    
}