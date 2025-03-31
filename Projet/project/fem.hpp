#pragma once
#include "defines.hpp"
#include <vector>

namespace fem {

typedef enum { FEM_EDGE = 2, FEM_TRIANGLE = 3, FEM_QUAD = 4 } ElementType;
typedef enum { DIRICHLET_X, DIRICHLET_Y, NEUMANN_X, NEUMANN_Y } BoundaryType;
typedef enum { PLANAR_STRESS, PLANAR_STRAIN, AXISYM } ElasticCase;

struct Node {
  double x;
  double y;
};

template <ElementType T> struct Element {
  // stores the nodes id (ie : their position in the mesh array);
  // the nodes are connected to the adjacent ones in the array (the array
  // "loops")
  std::array<int, T> element_nodes;
};

struct Discrete {
  void (*x)(double *xsi, double *eta); // gives the xsi and eta values to
                                       // integrate over the element
  void (*phi)(
      double xsi, double eta,
      double *phi); // the values of the form functions at the given points

  void (*dphidx)(double xsi, double eta, double *dphidxsi, double *dphideta);
};

template <ElementType T> struct Integration {
  const double xsi[T]; // looks too clean to change
  const double eta[T];
  const double weight[T];
};

template <ElementType T> struct Mesh {
  ElementType mesh_type = T;
  std::vector<Element<T>> elements_lists;
  Integration<T> *integration_rule;
  Discrete *functions;
};

struct Domain {
  int number_of_elements;
  std::vector<Element<FEM_EDGE> *> elements;
  char name[MAXNAME];
};

struct BoundaryCondition {
  Domain *domain;
  BoundaryType type;
  double value;
};

template <ElementType T> struct Geomerty {
  int number_of_nodes;
  std::vector<Node> nodes_list;
  Mesh<FEM_EDGE> *the_edge_mesh;
  Mesh<T> *the_inner_mesh;
  int number_of_domains;
  std::vector<Domain> domains_list;
};

struct System {
  std::vector<double> B;
  std::vector<std::vector<double>> A;
  int size;
};

template <ElementType T> struct Problem {
  double E, nu, rho, g;
  double A, B, C;
  int planarStrainStress; // dunno what it does
  int number_of_boundary_conditions;
  BoundaryCondition *conditions;
  int *constrained_nodes; // what the fuck is this
  Geomerty<T> *geometry;
  System *system;
};

} // namespace fem
