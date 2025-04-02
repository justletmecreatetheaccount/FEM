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
    return 0.03;
 
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