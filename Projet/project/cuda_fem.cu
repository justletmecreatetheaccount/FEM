#include "cuda_fem.cuh"
#include "defines.hpp"
#include "fem.hpp"

// Will be called one thread per element
__global__ void
assemble_system_kernel(unsigned int number_of_elements,
                       unsigned int *row_pointer, unsigned int *columns,
                       double *data, double *B, fem::Node *nodes_list,
                       fem::Element<fem::FEM_TRIANGLE> *elements_list,
                       double _A, double _B, double _C, double rho, double g,
                       double *_xsi, double *_eta, double *_weights) {

  if (blockIdx.x * blockDim.x + threadIdx.x < number_of_elements) {
    unsigned int nodes_id[3], nodes_sysX[3], nodes_sysY[3];
    double phi[3], dphidxsi[3], dphideta[3], dphidx[3], dphidy[3], nodes_x[3],
        nodes_y[3];

    for (int i = 0; i < 3; i++) {
      nodes_id[i] =
          elements_list[blockIdx.x * blockDim.x + threadIdx.x].element_nodes[i];
      nodes_sysX[i] = 2 * nodes_id[i];
      nodes_sysY[i] = 2 * nodes_id[i] + 1;
      nodes_x[i] = nodes_list[nodes_id[i]].x;
      nodes_y[i] = nodes_list[nodes_id[i]].y;
    }

    for (int iInteg = 0; iInteg < 3; iInteg++) {
      double xsi = _xsi[iInteg];
      double eta = _eta[iInteg];
      double weight = _weights[iInteg];
      phi[0] = 1 - xsi - eta;
      phi[1] = xsi;
      phi[2] = eta;

      dphidxsi[0] = -1.0;
      dphidxsi[1] = 1.0;
      dphidxsi[2] = 0.0;
      dphideta[0] = -1.0;
      dphideta[1] = 0.0;
      dphideta[2] = 1.0;

      double dxdxsi = 0.0;
      double dxdeta = 0.0;
      double dydxsi = 0.0;
      double dydeta = 0.0;
      for (int i = 0; i < 3; i++) {
        dxdxsi += nodes_x[i] * dphidxsi[i];
        dxdeta += nodes_x[i] * dphideta[i];
        dydxsi += nodes_y[i] * dphidxsi[i];
        dydeta += nodes_y[i] * dphideta[i];
      }
      double jac = fabs(dxdxsi * dydeta - dxdeta * dydxsi);

      for (int i = 0; i < 3; i++) {
        dphidx[i] = (dphidxsi[i] * dydeta - dphideta[i] * dydxsi) / jac;
        dphidy[i] = (dphideta[i] * dxdxsi - dphidxsi[i] * dxdeta) / jac;
      }
      for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
          unsigned int p = row_pointer[nodes_sysX[i]];
          unsigned int q = row_pointer[nodes_sysX[i] + 1];
          for (unsigned int k = p; k < q; k++) {
            if (columns[k] == nodes_sysX[j])
              atomicAdd(&(data[k]), (dphidx[i] * _A * dphidx[j] +
                                     dphidy[i] * _C * dphidy[j]) *
                                        jac * weight);
            if (columns[k] == nodes_sysY[j])
              atomicAdd(&(data[k]), (dphidx[i] * _B * dphidy[j] +
                                     dphidy[i] * _C * dphidx[j]) *
                                        jac * weight);
          }
          p = row_pointer[nodes_sysY[i]];
          q = row_pointer[nodes_sysY[i] + 1];
          for (unsigned int k = p; k < q; k++) {
            if (columns[k] == nodes_sysX[j])
              atomicAdd(&(data[k]), (dphidy[i] * _B * dphidx[j] +
                                     dphidx[i] * _C * dphidy[j]) *
                                        jac * weight);
            if (columns[k] == nodes_sysY[j])
              atomicAdd(&(data[k]), (dphidy[i] * _A * dphidy[j] +
                                     dphidx[i] * _C * dphidx[j]) *
                                        jac * weight);
          }
        }
      }
      for (int i = 0; i < 3; i++) {
        atomicAdd(&(B[nodes_sysY[i]]), -(phi[i] * g * rho * jac * weight));
      }
    }
  }
}

void cuda::fem::assemble_system(::fem::Problem &problem) {
  unsigned int number_of_elements =
      problem.geometry.full_mesh.elements_lists.size();
  unsigned int blocks_number = number_of_elements / THREADS_PER_BLOCK;
  blocks_number++;
  double *_xsi, *_eta, *_weights;
  ::fem::Node *nodes_list;
  ::fem::Element<::fem::FEM_TRIANGLE> *elements_list;
  cudaMalloc(&_xsi, sizeof(double) * 3);
  cudaMalloc(&_eta, sizeof(double) * 3);
  cudaMalloc(&_weights, sizeof(double) * 3);
  cudaMalloc(&nodes_list,
             sizeof(::fem::Node) * problem.geometry.number_of_nodes);
  cudaMalloc(&elements_list,
             sizeof(::fem::Element<::fem::FEM_TRIANGLE>) * number_of_elements);
  cudaMemcpy(_xsi, problem.geometry.full_mesh.xsi, 3 * sizeof(double),
             cudaMemcpyHostToDevice);
  cudaMemcpy(_eta, problem.geometry.full_mesh.eta, 3 * sizeof(double),
             cudaMemcpyHostToDevice);
  cudaMemcpy(_weights, problem.geometry.full_mesh.weights, 3 * sizeof(double),
             cudaMemcpyHostToDevice);
  cudaMemcpy(nodes_list, problem.geometry.nodes_list.data(),
             sizeof(::fem::Node) * problem.geometry.number_of_nodes,
             cudaMemcpyHostToDevice);
  cudaMemcpy(elements_list, problem.geometry.full_mesh.elements_lists.data(),
             sizeof(::fem::Element<::fem::FEM_TRIANGLE>) * number_of_elements,
             cudaMemcpyHostToDevice);
  assemble_system_kernel<<<blocks_number, THREADS_PER_BLOCK>>>(
      number_of_elements, problem.system.row_ptr.data(),
      problem.system.column.data(), problem.system.data.data(),
      problem.system.B.data(), nodes_list, elements_list, problem.A, problem.B,
      problem.C, problem.rho, problem.g, _xsi, _eta, _weights);
  cudaDeviceSynchronize();
  cudaFree(_xsi);
  cudaFree(_eta);
  cudaFree(_weights);
  cudaFree(nodes_list);
  cudaFree(elements_list);
};
