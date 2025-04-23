#include "fem.hpp"
#include "file_io.hpp"
#include <string>

/*
static const std::array<double, 3> _gaussTri3Xsi = {
    0.166666666666667, 0.666666666666667, 0.166666666666667};
static const std::array<double, 3> _gaussTri3Eta = {
    0.166666666666667, 0.166666666666667, 0.666666666666667};
static const std::array<double, 3> _gaussTri3Weight = {
    0.166666666666667, 0.166666666666667, 0.166666666666667};
static const std::array<double, 2> _gaussEdge2Xsi = {0.577350269189626,
                                                     -0.577350269189626};
static const std::array<double, 2> _gaussEdge2Weight = {1.000000000000000,
                                                        1.000000000000000};
*/
fem::Node::Node() : x(0), y(0) {};
fem::Node::Node(double x_n, double y_n) : x(x_n), y(y_n) {};

void fem::Mesh<fem::FEM_EDGE>::phi(double _xsi,
                                   std::array<double, FEM_EDGE> &phi) {
  phi[0] = (1 - _xsi) / 2.0;
  phi[1] = (1 + _xsi) / 2.0;
};

void fem::Mesh<fem::FEM_EDGE>::dphidx(double _xsi,
                                      std::array<double, FEM_EDGE> &dphidxsi) {
  dphidxsi[0] = -0.5;
  dphidxsi[1] = 0.5;
}
void fem::Mesh<fem::FEM_TRIANGLE>::phi(double _xsi, double _eta,
                                       std::array<double, FEM_TRIANGLE> &phi) {
  phi[0] = 1 - _xsi - _eta;
  phi[1] = _xsi;
  phi[2] = _eta;
};

void fem::Mesh<fem::FEM_TRIANGLE>::dphidx(
    double _xsi, double _eta, std::array<double, FEM_TRIANGLE> &dphidxsi,
    std::array<double, FEM_TRIANGLE> &dphideta) {
  dphidxsi[0] = -1.0;
  dphidxsi[1] = 1.0;
  dphidxsi[2] = 0.0;
  dphideta[0] = -1.0;
  dphideta[1] = 0.0;
  dphideta[2] = 1.0;
};

fem::Problem::Problem(const char *input_file, double _E, double _nu,
                      double _rho, double _g) {

  fem::read_mesh_file(input_file, (*this));
  E = _E;
  g = _g;
  nu = _nu;
  rho = _rho;

  A = E * (1 - nu) / ((1 + nu) * (1 - 2 * nu));
  B = E * nu / ((1 + nu) * (1 - 2 * nu));
  C = E / (2 * (1 + nu));
}
void fem::assemble_system(Problem problem) {}
