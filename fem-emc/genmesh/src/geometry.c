#include "fem.h"
#include <gmshc.h>

//
// Ici, vous pouvez définir votre géométrie :-)
//  (1) Raffiner intelligemment.... (yes )
//  (2) Construire la geometrie avec OpenCascade
//  (3) Construire la geometrie avec les outils de GMSH
//  (4) Obtenir la geometrie en lisant un fichier .geo de GMSH

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
