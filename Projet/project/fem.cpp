#include "fem.hpp"
#include "file_io.hpp"
#include <iostream>
#include <string>
#include <tuple>

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

fem::BoundaryCondition::BoundaryCondition(fem::Problem &problem,
                                          fem::Domain *_domain,
                                          fem::BoundaryType _type,
                                          double _value) {
  domain = _domain;
  type = _type;
  value = _value;

  int shift = -1;
  if (type == DIRICHLET_X)
    shift = 0;
  if (type == DIRICHLET_Y)
    shift = 1;
  if (shift == -1)
    return;

  for (unsigned int e = 0; e < domain->number_of_elements; e++) {
    for (unsigned int i = 0; i < fem::FEM_EDGE; i++) {
      int node = domain->elements[e]->element_nodes[i];
      problem.constrained_nodes.emplace_back(2 * node + shift, this);
    }
  }
};

void fem::BoundaryCondition::system_constrain(fem::System &system,
                                              unsigned int myNode,
                                              double myValue) {
  for (unsigned int i = 0; i < system.size; i++) {
    for (unsigned int j = system.row_ptr[i]; j < system.row_ptr[i + 1]; j++) {
      if (system.column[j] == myNode) {
        // Found the column index for the current row
        system.B[i] -= myValue * system.data[j];
        system.data[j] = 0;
        break;
      }
    }
  }

  for (unsigned int i = system.row_ptr[myNode]; i < system.row_ptr[myNode + 1];
       i++) {
    if (system.column[i] == myNode) {
      // Found the column index for the current row
      system.data[i] = 1;
    } else {
      system.data[i] = 0;
    }
  }
  system.B[myNode] = myValue;
}

void fem::Mesh<fem::FEM_EDGE>::phi(double _xsi, double phi[FEM_EDGE]) {
  phi[0] = (1 - _xsi) / 2.0;
  phi[1] = (1 + _xsi) / 2.0;
};

void fem::Mesh<fem::FEM_EDGE>::dphidx(double _xsi, double dphidxsi[FEM_EDGE]) {
  dphidxsi[0] = -0.5;
  dphidxsi[1] = 0.5;
}
void fem::Mesh<fem::FEM_TRIANGLE>::phi(double _xsi, double _eta,
                                       double phi[FEM_TRIANGLE]) {
  phi[0] = 1 - _xsi - _eta;
  phi[1] = _xsi;
  phi[2] = _eta;
};

void fem::Mesh<fem::FEM_TRIANGLE>::dphidx(double _xsi, double _eta,
                                          double dphidxsi[FEM_TRIANGLE],
                                          double dphideta[FEM_TRIANGLE]) {
  dphidxsi[0] = -1.0;
  dphidxsi[1] = 1.0;
  dphidxsi[2] = 0.0;
  dphideta[0] = -1.0;
  dphideta[1] = 0.0;
  dphideta[2] = 1.0;
};

fem::Domain *fem::Geomerty::get_domain_by_name(std::string _name) {
  for (unsigned int i = 0; i < number_of_domains; i++) {
    if (_name == domains_list[i].name) {
      return &domains_list[i];
    }
  }
  std::cout << "cant find domain\n";
  return nullptr;
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

void fem::Problem::femSolutionWrite(int nNodes, int nfields, double *data,
                                    const char *filename) {
  FILE *file = fopen(filename, "w");
  if (!file) {
    printf("Error at %s:%d\nUnable to open file %s\n", __FILE__, __LINE__,
           filename);
    exit(-1);
  }
  fprintf(file, "Size %d,%d\n", nNodes, nfields);
  for (int i = 0; i < nNodes; i++) {
    for (int j = 0; j < nfields - 1; j++) {
      fprintf(file, "%.18le,", data[i * nfields + j]);
    }
    fprintf(file, "%.18le", data[i * nfields + nfields - 1]);
    fprintf(file, "\n");
  }
  fclose(file);
}
