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
    femGeo* theGeometry = geoGetGeometry();
    double h = theGeometry->h;
    double w = theGeometry->LxPlate;
    double plateHeight = theGeometry->LyPlate;
    
    // Hole parameters
    double holeWidth = w/2.0, holeHeight = plateHeight/4.0;
    double lowerY = -plateHeight/4.0, upperY = plateHeight/4.0;
    double cornerRadius = w/20.0, transWidth = w/8.0, hHoles = h/2.5;
    double result = h;
    
    // Calculate distances to both holes
    for (int hole = 0; hole < 2; hole++) {
        double holeY = (hole == 0) ? lowerY : upperY;
        double dx = fabs(x) - holeWidth/2.0 + cornerRadius;
        double dy = fabs(y - holeY) - holeHeight/2.0 + cornerRadius;
        
        double distHole;
        if (dx > 0 && dy > 0)
            distHole = sqrt(dx*dx + dy*dy) - cornerRadius;
        else if (dx > 0)
            distHole = dx - cornerRadius;
        else if (dy > 0)
            distHole = dy - cornerRadius;
        else
            distHole = -fmin(-dx, -dy);
        
        // Apply size transition
        if (distHole <= 0) {
            result = hHoles;
        } else if (distHole <= transWidth) {
            double d = distHole;
            double a2 = 3 * (h - hHoles) / pow(transWidth, 2);
            double a3 = 2 * (hHoles - h) / pow(transWidth, 3);
            double sizeNear = hHoles + a2 * pow(d, 2) + a3 * pow(d, 3);
            result = fmin(result, sizeNear);
        }
    }
    
    return result;
}

void geoMeshGenerate(double lc)
{
    femGeo *theGeometry = geoGetGeometry();
    double Lx = 1.0, Ly = 2.0;
    
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
