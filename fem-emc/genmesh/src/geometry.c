#include "fem.h"
#include <gmshc.h>

//
// Ici, vous pouvez définir votre géométrie :-)
//  (1) Raffiner intelligemment.... (yes )
//  (2) Construire la geometrie avec OpenCascade
//  (3) Construire la geometrie avec les outils de GMSH
//  (4) Obtenir la geometrie en lisant un fichier .geo de GMSH

static double distanceToRectangleInside(double x, double y,
                                        double cx, double cy,
                                        double halfW, double halfH)
{
    // Shift to rectangle center
    double rx = x - cx;
    double ry = y - cy;

    // Distance to each boundary in x and y directions from the inside
    // e.g., halfW - fabs(rx) is how far x is from the vertical boundary
    double dx = halfW - fabs(rx);
    double dy = halfH - fabs(ry);

    // If outside or exactly on boundary, distance is zero (or negative if out)
    if (dx < 0.0) dx = 0.0;
    if (dy < 0.0) dy = 0.0;

    // The distance to the rectangle's boundary for an inside point is the smaller of dx and dy
    return (dx < dy) ? dx : dy;
}

/*
 * Computes the distance from a point (x, y) to the *outside* of a rectangle boundary
 * (e.g., used for holes). If (x, y) is outside the rectangle, this distance
 * can be nonzero. If inside, it is clamped to 0.0.
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
    double h = theGeometry->h;
    double w = theGeometry->LxPlate;  // plate width
    double l = theGeometry->LyPlate;  // plate height

    // 1) First interior rectangle ("Notch")
    double x0  = 0.0;         // center X
    double y0  = -l / 4.0;     // center Y
    double w0  = w / 2.0;      // rectangle width
    double ht0 = l / 4.0;      // rectangle height
    double d0  = 0.5;          // transition zone thickness
    double h0  = 0.05;         // smaller mesh size near the notch

    // 2) Second interior rectangle ("Hole")
    double x1  = 0.0;
    double y1  = (2.0 * l) / 8.0;
    double w1  = w / 2.0;
    double ht1 = l / 4.0; 
    double d1  = 0.5;
    double h1  = 0.05;

    // 3) Outer boundary rectangle (overall plate)
    //    We want to refine the mesh near the edges from inside
    double x2  = 0.0;      // plate center X
    double y2  = 0.0;      // plate center Y
    double w2  = w / 2.0;  // half-width
    double h2_ = l / 2.0;  // half-height
    double d2  = 0.5;      // transition thickness at the boundary
    double r2  = 0.03;     // smaller mesh size near the outer boundary

    // Default mesh size is the larger size
    double result = h;

    // Distance from the first rectangle boundary (treated as an 'outside' type: hole or notch)
    double dist0 = distanceToRectangleOutside(x, y, x0, y0, w0 / 2.0, ht0 / 2.0);
    if (dist0 <= d0) {
        // Hermite interpolation from h (far) to h0 (near)
        double a2 = (3.0 * (h - h0)) / (d0 * d0);
        double a3 = (2.0 * (h0 - h)) / (d0 * d0 * d0);
        double val0 = h0 + a2 * (dist0 * dist0) + a3 * (dist0 * dist0 * dist0);
        result = val0;
    }

    // Distance from the second rectangle boundary (also an 'outside' type)
    double dist1 = distanceToRectangleOutside(x, y, x1, y1, w1 / 2.0, ht1 / 2.0);
    if (dist1 <= d1) {
        double a2 = (3.0 * (h - h1)) / (d1 * d1);
        double a3 = (2.0 * (h1 - h)) / (d1 * d1 * d1);
        double val1 = h1 + a2 * (dist1 * dist1) + a3 * (dist1 * dist1 * dist1);
        if (val1 < result) {
            result = val1;
        }
    }



    return result;
}

void geoMeshGenerate(double lc)
{
    femGeo *theGeometry = geoGetGeometry();
    double Lx = 4.0, Ly = 8.0;
    
    // Setup geometry parameters
    theGeometry->LxPlate = Lx;
    theGeometry->LyPlate = Ly;
    theGeometry->h = Lx * lc;
    theGeometry->elementType = FEM_TRIANGLE;
    geoSetSizeCallback(geoSize);
    
    int ierr;
    
    // Create geometry
    int idRect = gmshModelOccAddRectangle(-Lx/2.0, -Ly/2.0, 0.0, Lx, Ly, -1, 0.0, &ierr);
    int idHoleDown = gmshModelOccAddRectangle(-Lx/4.0, -3*Ly/8.0, 0.0, Lx/2.0, Ly/4.0, -1, 0.1, &ierr);
    int idHoleUp = gmshModelOccAddRectangle(-Lx/4.0, Ly/8.0, 0.0, Lx/2.0, Ly/4.0, -1, 0.1, &ierr);
    
    // Cut the holes from the rectangle
    int rect[] = {2, idRect};
    int holeDown[] = {2, idHoleDown};
    int holeUp[] = {2, idHoleUp};
    gmshModelOccCut(rect, 2, holeDown, 2, NULL, NULL, NULL, NULL, NULL, -1, 1, 1, &ierr);
    ErrorGmsh(ierr);
    gmshModelOccCut(rect, 2, holeUp, 2, NULL, NULL, NULL, NULL, NULL, -1, 1, 1, &ierr);
    ErrorGmsh(ierr);
    gmshModelOccSynchronize(&ierr);
    
    // Generate mesh
    gmshOptionSetNumber("Mesh.Algorithm", 6, &ierr);
    gmshOptionSetNumber("Mesh.SaveAll", 1, &ierr);
    gmshModelMeshGenerate(2, &ierr);
}

void geoMeshGenerateGeo(void)
{
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

void geoMeshGenerateGeoFile(const char *filename)
{
    femGeo *theGeometry = geoGetGeometry();
    int ierr;
    gmshOpen(filename, &ierr);
    ErrorGmsh(ierr);

    gmshOptionSetNumber("Mesh.Algorithm", 3, &ierr);
    gmshOptionSetNumber("Mesh.SaveAll", 1, &ierr);
    gmshModelMeshGenerate(2, &ierr);
    return;
}

void geoMeshGenerateMshFile(const char *filename)
{
    int ierr;
    gmshOpen(filename, &ierr);
    ErrorGmsh(ierr);
    return;
}
