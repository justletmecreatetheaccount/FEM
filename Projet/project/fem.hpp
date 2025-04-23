#pragma once
#include "cuda_utils.cuh"
#include <array>
#include <memory>
#include <string>
#include <vector>

namespace fem {

typedef enum { FEM_EDGE = 2, FEM_TRIANGLE = 3 } ElementType;
typedef enum { DIRICHLET_X, DIRICHLET_Y } BoundaryType;

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

template <ElementType T> struct Integration {
  // should be const but cpp doesn't allow it :
  // https://stackoverflow.com/questions/14495536/how-do-i-initialize-a-const-data-member
  std::array<double, T> xi;
  std::array<double, T> eta;
  std::array<double, T> weight;
};

template <ElementType T> struct Mesh {};

template <> struct Mesh<FEM_EDGE> {
  ElementType mesh_type = FEM_EDGE;
  std::vector<Element<FEM_EDGE>> elements_lists;
  Integration<FEM_EDGE> integration_rule;
  // xsi and weights values to integrate over the element
  const std::array<double, FEM_EDGE> xsi = {0.577350269189626,
                                            -0.577350269189626};
  const std::array<double, FEM_EDGE> weights = {1.000000000000000,
                                                1.000000000000000};
  void phi(double _xsi,
           std::array<double, FEM_EDGE>
               &phi); // the values of the form functions at the given points

  void dphidx(double _xsi, std::array<double, FEM_EDGE> &dphidxsi);
};
template <> struct Mesh<FEM_TRIANGLE> {
  ElementType mesh_type = FEM_TRIANGLE;
  std::vector<Element<FEM_TRIANGLE>> elements_lists;
  Integration<FEM_TRIANGLE> integration_rule;
  // xsi, eta and weights values to integrate over the element
  const std::array<double, FEM_TRIANGLE> xsi = {
      0.166666666666667, 0.666666666666667, 0.166666666666667};
  const std::array<double, FEM_TRIANGLE> eta = {
      0.166666666666667, 0.166666666666667, 0.666666666666667};
  const std::array<double, FEM_TRIANGLE> weights = {
      0.166666666666667, 0.166666666666667, 0.166666666666667};
  void
  phi(double _xsi, double _eta,
      std::array<double, FEM_TRIANGLE> &phi); // the values of the form
                                              // functions at the given points

  void dphidx(double _xsi, double _eta,
              std::array<double, FEM_TRIANGLE> &dphidxsi,
              std::array<double, FEM_TRIANGLE> &dphideta);
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
  std::vector<double, cuda::allocator<double>> B;
  std::vector<double, cuda::allocator<double>> data;
  std::vector<int, cuda::allocator<int>> column;
  std::vector<int, cuda::allocator<int>> row_ptr;
  int size;
};

struct Problem {
  double E, nu, rho, g;
  double A, B, C;
  int number_of_boundary_conditions;
  std::vector<BoundaryCondition> conditions;
  // int *constrained_nodes; // what the fuck is this (array)
  Geomerty geometry;
  System system;

  Problem(const char *input_file, double _E, double _nu, double _rho,
          double _g);
};

/* TODO
 */
void assemble_system(Problem problem);
void conjugate_gradient(System system);
void incomplete_cholesky();
} // namespace fem
