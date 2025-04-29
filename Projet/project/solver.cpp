#include "../benchmark.h"
#include "cuda_fem.cuh"
#include "cuda_utils.cuh"
#include "fem.hpp"
#include "project/defines.hpp"
#include <bits/stdc++.h>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <vector>

void elasticity_solve(const char *meshfile, const char *outfile, double E,
                      double nu, double rho, double g) {
  cuda::test();
  //  Read the mesh and create the problem
  fem::Problem problem = fem::Problem(meshfile, E, nu, rho, g);

  // Boundary conditions are Dirichlet X and Y
  problem.conditions.emplace_back(
      std::ref(problem), problem.geometry.get_domain_by_name("symmetry"),
      fem::DIRICHLET_X, 0.0);
  problem.conditions.emplace_back(
      std::ref(problem), problem.geometry.get_domain_by_name("symmetry"),
      fem::DIRICHLET_Y, 0.0);
  problem.conditions.emplace_back(std::ref(problem),
                                  problem.geometry.get_domain_by_name("base"),
                                  fem::DIRICHLET_X, 0.0);
  problem.conditions.emplace_back(std::ref(problem),
                                  problem.geometry.get_domain_by_name("base"),
                                  fem::DIRICHLET_Y, 0.0);

  cuda::fem::assemble_system(problem);
  for (unsigned int i = 0; i < problem.constrained_nodes.size(); i++) {
    fem::BoundaryCondition::system_constrain(
        problem.system, std::get<0>(problem.constrained_nodes[i]),
        std::get<1>(problem.constrained_nodes[i])->value);
  }

  // cuda::conjugate_gradient(problem.system.size, problem.system.column.data(),
  //                          problem.system.row_ptr.data(),
  //                          problem.system.data.data(),
  //                          problem.system.B.data());
  // fem::Problem::femSolutionWrite(problem.geometry.number_of_nodes, 2,
  //                                problem.system.B.data(), outfile);

  // Assemble and solve
  // femElasticityPrint(theProblem);
  //

  //   Write out the solution
  // for (unsigned int i = 0; i < problem.system.size; i++) {
  //   std::cout << "B [" << i << "] : " << problem.system.B[i] << "\n";
  // }
}
