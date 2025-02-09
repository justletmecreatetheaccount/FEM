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


double geoSize(double x, double y) {
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
    
    double result = h;

    // Getting the profile unidimensionally with d
    double _d0 = sqrt(pow(x - x0, 2) + pow(y - y0, 2)) - r0; // Distance to the notch border
    double _d1 = sqrt(pow(x - x1, 2) + pow(y - y1, 2)) - r1; // Distance to the hole border
    if (_d0 <= d0) { // If the points are in the radius of r0 + d0 which is the transition zone
        // Do Hermite Interpolation
        double a2 = (3 * (h - h0)) / pow(d0, 2);
        double a3 = (2 * (h0 - h)) / pow(d0, 3);
        result = h0 + a2 * pow(_d0, 2) + a3 * pow(_d0, 3);
    }
    if (_d1 <= d1) { // Do not use elif in the case of the two transition zones overlapping
        // Second Hermite Interpolation
        double a2 = (3 * (h - h1)) / pow(d1, 2);
        double a3 = (2 * (h1 - h)) / pow(d1, 3);
        result = fmin(result,h1 + a2 * pow(_d1, 2) + a3 * pow(_d1, 3)); // Min for overlapping transition zones
    }
        
    return result;
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
//  -1- Construction de la g�om�trie avec OpenCascade
//      On cr�e le rectangle
//      On cr�e les deux cercles
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
//  -2- D�finition de la fonction callback pour la taille de r�f�rence
//      Synchronisation de OpenCascade avec gmsh
//      G�n�ration du maillage (avec l'option Mesh.SaveAll :-)
                  
   
    geoSetSizeCallback(geoSize);
                                  
    gmshModelOccSynchronize(&ierr);       
    gmshOptionSetNumber("Mesh.SaveAll", 1, &ierr);
    gmshModelMeshGenerate(2, &ierr);  
       
//
//  Generation de quads :-)
//
//    gmshOptionSetNumber("Mesh.SaveAll", 1, &ierr);
//    gmshOptionSetNumber("Mesh.RecombineAll", 1, &ierr);
//    gmshOptionSetNumber("Mesh.Algorithm", 8, &ierr);  //chk(ierr);
//    gmshOptionSetNumber("Mesh.RecombinationAlgorithm", 1.0, &ierr);  //chk(ierr);
//    gmshModelGeoMeshSetRecombine(2,1,45,&ierr);  //chk(ierr);
//    gmshModelMeshGenerate(2, &ierr);  
   
 
//
//  Plot of Fltk
//
//   gmshFltkInitialize(&ierr);
//   gmshFltkRun(&ierr);  //chk(ierr);
//
    
}