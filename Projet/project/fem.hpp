#pragma once
#include "cuda_utils.cuh"
#include <array>
#include <memory>
#include <string>
#include <tuple>
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
  unsigned int element_nodes[T];
};

template <ElementType T> struct Mesh {};

template <> struct Mesh<FEM_EDGE> {
  ElementType mesh_type = FEM_EDGE;
  std::vector<Element<FEM_EDGE>> elements_lists;
  // xsi and weights values to integrate over the element
  double xsi[FEM_EDGE] = {0.577350269189626, -0.577350269189626};
  double weights[FEM_EDGE] = {1.000000000000000, 1.000000000000000};
  void phi(double _xsi,
           double phi[FEM_EDGE]); // the values of the form functions at the
                                  // given points

  void dphidx(double _xsi, double dphidxsi[FEM_EDGE]);
};
template <> struct Mesh<FEM_TRIANGLE> {
  ElementType mesh_type = FEM_TRIANGLE;
  std::vector<Element<FEM_TRIANGLE>> elements_lists;
  // xsi, eta and weights values to integrate over the element
  double xsi[FEM_TRIANGLE] = {0.166666666666667, 0.666666666666667,
                              0.166666666666667};
  double eta[FEM_TRIANGLE] = {0.166666666666667, 0.166666666666667,
                              0.666666666666667};
  double weights[FEM_TRIANGLE] = {0.166666666666667, 0.166666666666667,
                                  0.166666666666667};
  void phi(double _xsi, double _eta,
           double phi[FEM_TRIANGLE]); // the values of the form
                                      // functions at the given points

  void dphidx(double _xsi, double _eta, double dphidxsi[FEM_TRIANGLE],
              double dphideta[FEM_TRIANGLE]);
};

struct System {
  std::vector<double, cuda::allocator<double>> B;
  std::vector<double, cuda::allocator<double>> data;
  std::vector<unsigned int, cuda::allocator<unsigned int>> column;
  std::vector<unsigned int, cuda::allocator<unsigned int>> row_ptr;
  unsigned int size;
};

struct Domain {
  unsigned int number_of_elements;
  std::vector<Element<FEM_EDGE> *> elements;
  std::string name;
};

struct Problem;

struct BoundaryCondition {
  const Domain *domain; // no smart pointer bc no ownership
  BoundaryType type;
  double value;
  BoundaryCondition(Problem &problem, Domain *_domain, BoundaryType _type,
                    double _value);
  static void system_constrain(fem::System &system, unsigned int myNode,
                               double myValue);
};

struct Geomerty {
  unsigned int number_of_nodes; // redundant but pretty
  std::vector<Node> nodes_list;
  Mesh<FEM_EDGE> edge_mesh;
  Mesh<FEM_TRIANGLE> full_mesh;
  unsigned int number_of_domains;
  std::vector<Domain> domains_list;
  Domain *get_domain_by_name(std::string _name);
};

struct Problem {
  double E, nu, rho, g;
  double A, B, C;
  std::vector<BoundaryCondition> conditions;
  std::vector<std::tuple<unsigned int, BoundaryCondition *>>
      constrained_nodes; // array that contains the positions in the system
                         // (node_id x or y) and pointer to associated condition
  Geomerty geometry;
  System system;

  Problem(const char *input_file, double _E, double _nu, double _rho,
          double _g);
};

/* TODO
 */
} // namespace fem
