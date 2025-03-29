#include "fem.h"
#include <gmshc.h>

//
// Ici, vous pouvez définir votre géométrie :-)
//  (1) Raffiner intelligemment.... (yes )
//  (2) Construire la geometrie avec OpenCascade
//  (3) Construire la geometrie avec les outils de GMSH
//  (4) Obtenir la geometrie en lisant un fichier .geo de GMSH

double geoSize(double x, double y)
{
    femGeo *theGeometry = geoGetGeometry();
    double h = theGeometry->h;
    double w = theGeometry->LxPlate;
    double plateHeight = theGeometry->LyPlate;
    
    // Smoothing parameters (tune to adjust interpolation steepness)
    double holeSmooth = 1.0;    // Amplifies hole transition
    double edgeSmooth = 1.0;    // Amplifies edge transition
    
    // Hole settings
    double lowerHoleW = w / 2.0, lowerHoleH = plateHeight / 4.0;
    double lowerHoleX = 0.0,    lowerHoleY = -plateHeight / 4.0;
    double upperHoleW = w / 2.0, upperHoleH = plateHeight / 4.0;
    double upperHoleX = 0.0,    upperHoleY =  plateHeight / 4.0;
    double cR = w / 20.0;
    double transW = w / 8.0, hHoles = h / 2.5;
    
    // Border settings
    double borderRefine = h / 2.5;
    double borderTransW = w / 6.0;
    
    // Start from base size
    double result = h;
    
    // Distances to lower hole
    double dx = fabs(x - lowerHoleX) - lowerHoleW / 2.0 + cR;
    double dy = fabs(y - lowerHoleY) - lowerHoleH / 2.0 + cR;
    double distLower = (dx > 0 && dy > 0) ? sqrt(dx*dx + dy*dy) - cR
                     : (dx > 0)           ? dx - cR
                     : (dy > 0)           ? dy - cR
                     : -fmin(-dx, -dy);
    
    // Distances to upper hole
    dx = fabs(x - upperHoleX) - upperHoleW / 2.0 + cR;
    dy = fabs(y - upperHoleY) - upperHoleH / 2.0 + cR;
    double distUpper = (dx > 0 && dy > 0) ? sqrt(dx*dx + dy*dy) - cR
                     : (dx > 0)           ? dx - cR
                     : (dy > 0)           ? dy - cR
                     : -fmin(-dx, -dy);
    
    // Interpolation near lower hole
    if (distLower <= transW && distLower > 0) {
        double d = distLower;
        double a2 = holeSmooth * (3.0 * (h - hHoles) / (transW * transW));
        double a3 = holeSmooth * (2.0 * (hHoles - h) / (transW * transW * transW));
        double sizeNear = hHoles + a2 * (d * d) + a3 * (d * d * d);
        result = fmin(result, sizeNear);
    } else if (distLower <= 0) {
        result = hHoles;
    }
    
    // Interpolation near upper hole
    if (distUpper <= transW && distUpper > 0) {
        double d = distUpper;
        double a2 = holeSmooth * (3.0 * (h - hHoles) / (transW * transW));
        double a3 = holeSmooth * (2.0 * (hHoles - h) / (transW * transW * transW));
        double sizeNear = hHoles + a2 * (d * d) + a3 * (d * d * d);
        result = fmin(result, sizeNear);
    } else if (distUpper <= 0) {
        result = hHoles;
    }
    
    // Interpolation near borders
    double distLeft   = fabs(x + w / 2.0);
    double distRight  = fabs(x - w / 2.0);
    double distBottom = fabs(y + plateHeight / 2.0);
    double distTop    = fabs(y - plateHeight / 2.0);
    double distEdge   = fmin(fmin(distLeft, distRight), fmin(distBottom, distTop));
    
    if (distEdge <= borderTransW) {
        double d = distEdge;
        double a2 = edgeSmooth * (3.0 * (h - borderRefine) / (borderTransW * borderTransW));
        double a3 = edgeSmooth * (2.0 * (borderRefine - h) / (borderTransW * borderTransW * borderTransW));
        double sizeEdge = borderRefine + a2 * (d * d) + a3 * (d * d * d);
        result = fmin(result, sizeEdge);
    }
    
    return result;
}

void geoMeshGenerate(double lc) {
  femGeo *theGeometry = geoGetGeometry();
  double Lx = 1.0;
  double Ly = 2.0;
  theGeometry->LxPlate = Lx;
  theGeometry->LyPlate = Ly;
  theGeometry->h = Lx * lc;
  theGeometry->elementType = FEM_TRIANGLE;

  geoSetSizeCallback(geoSize);

  double w = theGeometry->LxPlate;
  double h = theGeometry->LyPlate;

  int ierr;
  double r = w / 4;
  int idRect = gmshModelOccAddRectangle(-w/2.0, -h/2.0, 0.0, w, h, -1, 0.1, &ierr);
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

  // Use a frontal delaunay algorithm
  gmshOptionSetNumber("Mesh.Algorithm", 6, &ierr);
  gmshOptionSetNumber("Mesh.SaveAll", 1, &ierr);
  gmshModelMeshGenerate(2, &ierr);
  
  //gmshFltkInitialize(&ierr);
  //gmshFltkRun(&ierr);  //chk(ierr);
  return;
}

void geoMeshGenerateGeo(void) {
  femGeo *theGeometry = geoGetGeometry();
  double Lx = 1.0;
  double Ly = 1.0;
  theGeometry->LxPlate = Lx;
  theGeometry->LyPlate = Ly;
  theGeometry->h = Lx * 0.05;
  theGeometry->elementType = FEM_TRIANGLE;

  geoSetSizeCallback(geoSize);

  /*
  4 ------------------ 3
  |                    |
  |                    |
  5 ------- 6          |
             \         |
              )        |
             /         |
  8 ------- 7          |
  |                    |
  |                    |
  1 ------------------ 2
  */

  int ierr;
  double w = theGeometry->LxPlate;
  double h = theGeometry->LyPlate;
  double r = w / 4;
  double lc = theGeometry->h;

  int p1 = gmshModelGeoAddPoint(-w / 2, -h / 2, 0., lc, 1, &ierr);
  int p2 = gmshModelGeoAddPoint(w / 2, -h / 2, 0., lc, 2, &ierr);
  int p3 = gmshModelGeoAddPoint(w / 2, h / 2, 0., lc, 3, &ierr);
  int p4 = gmshModelGeoAddPoint(-w / 2, h / 2, 0., lc, 4, &ierr);
  int p5 = gmshModelGeoAddPoint(-w / 2, r, 0., lc, 5, &ierr);
  int p6 = gmshModelGeoAddPoint(0., r, 0., lc, 6, &ierr);
  int p7 = gmshModelGeoAddPoint(0., -r, 0., lc, 7, &ierr);
  int p8 = gmshModelGeoAddPoint(-w / 2, -r, 0., lc, 8, &ierr);
  int p9 = gmshModelGeoAddPoint(0., 0., 0., lc, 9, &ierr); // center of circle

  int l1 = gmshModelGeoAddLine(p1, p2, 1, &ierr);
  int l2 = gmshModelGeoAddLine(p2, p3, 2, &ierr);
  int l3 = gmshModelGeoAddLine(p3, p4, 3, &ierr);
  int l4 = gmshModelGeoAddLine(p4, p5, 4, &ierr);
  int l5 = gmshModelGeoAddLine(p5, p6, 5, &ierr);
  int l6 = gmshModelGeoAddCircleArc(p7, p9, p6, 6, 0., 0., 0., &ierr); // NB : the direction of the curve is reversed
  int l7 = gmshModelGeoAddLine(p7, p8, 7, &ierr);
  int l8 = gmshModelGeoAddLine(p8, p1, 8, &ierr);

  int lTags[] = {l1, l2, l3, l4, l5, -l6, l7, l8}; // NB : "-l6" because the curve is reversed
  int c1[] = {1};
  c1[0] = gmshModelGeoAddCurveLoop(lTags, 8, 1, 0, &ierr);
  int s1 = gmshModelGeoAddPlaneSurface(c1, 1, 1, &ierr);
  gmshModelGeoSynchronize(&ierr);

  gmshOptionSetNumber("Mesh.Algorithm", 3, &ierr);
  gmshOptionSetNumber("Mesh.SaveAll", 1, &ierr);
  gmshModelMeshGenerate(2, &ierr);
}

void geoMeshGenerateGeoFile(const char *filename) {
  femGeo *theGeometry = geoGetGeometry();
  int ierr;
  gmshOpen(filename, &ierr);
  ErrorGmsh(ierr);

  gmshOptionSetNumber("Mesh.Algorithm", 3, &ierr);
  gmshOptionSetNumber("Mesh.SaveAll", 1, &ierr);
  gmshModelMeshGenerate(2, &ierr);
  return;
}

void geoMeshGenerateMshFile(const char *filename) {
  int ierr;
  gmshOpen(filename, &ierr);
  ErrorGmsh(ierr);
  return;
}
