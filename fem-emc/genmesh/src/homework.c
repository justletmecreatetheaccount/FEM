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



static double distanceToRectangleOutside(double x, double y,
                                         double cx, double cy,
                                         double halfW, double halfH)
{
    double dx = fabs(x - cx) - halfW;
    double dy = fabs(y - cy) - halfH;
    if (dx < 0.0) dx = 0.0;
    if (dy < 0.0) dy = 0.0;
    return sqrt(dx * dx + dy * dy);
}

double geoSize(double x, double y)
{
    femGeo* theGeometry = geoGetGeometry();
    
    // Base mesh size
    double h  = theGeometry->h;
    // Domain
    double w  = theGeometry->LxPlate;
    double l  = theGeometry->LyPlate;
    
    // Centered rectangle #1
    double x0  = 0.0;
    double y0  = -l / 4.0;
    double w0  = w / 2.0;
    double ht0 = l / 4.0;
    double d0  = 0.5;   // Transition zone #1
    double h0  = 0.02;  // Smaller mesh near rectangle #1
    
    // Centered rectangle #2
    double x1  = 0.0;
    double y1  = (2.0 * l) / 8.0;
    double w1  = w / 2.0;
    double ht1 = l / 4.0; 
    double d1  = 0.5;   // Transition zone #2
    double h1  = 0.02;  // Smaller mesh near rectangle #2
    
    // Border transition
    double d2 = 0.2;  // Transition zone near domain edges
    double h2 = 0.02;  // Smaller mesh near domain edges

    double result = h;  // Start with base mesh size
    
    // Calculate distance to rectangle #1
    double dist0 = distanceToRectangleOutside(x, y, x0, y0, w0 / 2.0, ht0 / 2.0);
    if (dist0 <= d0) {
        // Hermite interpolation
        double a2 = 3.0 * (h - h0) / (d0 * d0);
        double a3 = 2.0 * (h0 - h) / (d0 * d0 * d0);
        double val0 = h0 + a2 * (dist0 * dist0) + a3 * (dist0 * dist0 * dist0);
        if (val0 < result) {
            result = val0;
        }
    }
    
    // Calculate distance to rectangle #2
    double dist1 = distanceToRectangleOutside(x, y, x1, y1, w1 / 2.0, ht1 / 2.0);
    if (dist1 <= d1) {
        double a2 = 3.0 * (h - h1) / (d1 * d1);
        double a3 = 2.0 * (h1 - h) / (d1 * d1 * d1);
        double val1 = h1 + a2 * (dist1 * dist1) + a3 * (dist1 * dist1 * dist1);
        if (val1 < result) {
            result = val1;
        }
    }
    
    // Distance to domain edges
    double distLeft   = x + w / 2.0;
    double distRight  = (w / 2.0) - x;
    double distBottom = y + l / 2.0;
    double distTop    = (l / 2.0) - y;
    double distEdge   = fmin(fmin(distLeft, distRight), fmin(distBottom, distTop));
    
    // Hermite interpolation near edges
    if (distEdge <= d2) {
        double a2 = (3.0 * (h - h2)) / (d2 * d2);
        double a3 = (2.0 * (h2 - h)) / (d2 * d2 * d2);
        double val2 = h2 + a2 * (distEdge * distEdge) + a3 * (distEdge * distEdge * distEdge);
        if (val2 < result) {
            result = val2;
        }
    }
    
    return result;
}


#define ___ 0

void geoMeshGenerate() {

  femGeo *theGeometry = geoGetGeometry();
  double Lx = 1.0;
  double Ly = 2.0;
  theGeometry->LxPlate = Lx;
  theGeometry->LyPlate = Ly;
  theGeometry->h = Lx * 0.1;

  geoSetSizeCallback(geoSize);

  double w = theGeometry->LxPlate;
  double h = theGeometry->LyPlate;

  int ierr;
  double r = w / 4;
  int idRect = gmshModelOccAddRectangle(-w/2.0, -h/2.0, 0.0, w, h, -1, 0.0, &ierr);
  int idHoleDown = gmshModelOccAddRectangle(-w/4.0, -3*h/8.0, 0.0, w/2.0, h/4.0, -1, 0.1, &ierr);
  int idHoleUp = gmshModelOccAddRectangle(-w/4.0, h/8.0, 0.0, w/2.0, h/4.0, -1, 0.1, &ierr);
    

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