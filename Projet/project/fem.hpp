#pragma once
#include <array>
#include <memory>
#include <string>
#include <vector>

namespace fem {

typedef enum { FEM_EDGE = 2, FEM_TRIANGLE = 3 } ElementType;
typedef enum { DIRICHLET_X, DIRICHLET_Y, NEUMANN_X, NEUMANN_Y } BoundaryType;
typedef enum { PLANAR_STRESS, PLANAR_STRAIN, AXISYM } ElasticCase;

struct Node {
  double x;
  double y;

public:
  Node();
  Node(double x_n, double y_n);
};

template <ElementType T> struct Element {
  // stores the nodes id (ie : their position in the mesh array);
  // the nodes are connected to the adjacent ones in the array (the array
  // "loops")
  std::array<int, T> element_nodes;
};

/*Utter Garbaje that should be removed*/
struct Discrete {
  void (*x)(double *xsi, double *eta); // gives the xsi and eta values to
                                       // integrate over the element
  void (*phi)(
      double xsi, double eta,
      double *phi); // the values of the form functions at the given points

  void (*dphidx)(double xsi, double eta, double *dphidxsi, double *dphideta);
};

template <ElementType T> struct Integration {
  // should be const but cpp doesn't allow it :
  // https://stackoverflow.com/questions/14495536/how-do-i-initialize-a-const-data-member
  std::array<double, T> xi;
  std::array<double, T> eta;
  std::array<double, T> weight;
};

template <ElementType T> struct Mesh {
  ElementType mesh_type = T;
  std::vector<Element<T>> elements_lists;
  Integration<T> integration_rule;
  Discrete functions;
};

struct Domain {
  int number_of_elements;
  std::vector<Element<FEM_EDGE> *> elements;
  std::string name;
};

struct BoundaryCondition {
  const Domain &domain; // no smart pointer bc no ownership
  BoundaryType type;
  double value;
};

struct Geomerty {
  int number_of_nodes; // redundant but pretty
  std::vector<Node> nodes_list;
  Mesh<FEM_EDGE> edge_mesh;
  Mesh<FEM_TRIANGLE> full_mesh;
  int number_of_domains;
  std::vector<Domain> domains_list;
};

struct System {
  std::vector<double> B;
  std::vector<double> data;
  std::vector<int> column;
  std::vector<int> row_ptr;
  int size;
};

struct Problem {
  double E, nu, rho, g;
  double A, B, C;
  int planarStrainStress; // dunno what it does
  int number_of_boundary_conditions;
  std::vector<BoundaryCondition> conditions;
  // int *constrained_nodes; // what the fuck is this (array)
  Geomerty geometry;
  System system;
};

} // namespace fem
