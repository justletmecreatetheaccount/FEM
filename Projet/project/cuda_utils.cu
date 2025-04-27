#include "cuda_utils.cuh"
#include "defines.hpp"
#include <iostream>

// Kernel function to add the elements of two arrays
template <typename T>
__global__ void add_vv(unsigned int n, T *v1, T *v2, T c1, T c2) {
  unsigned int idx = blockIdx.x * blockDim.x + threadIdx.x;
  if (idx < n)
    v2[idx] = c1 * v1[idx] + c2 * v2[idx];
}

template <typename T>
__global__ void multiply_vv(unsigned int n, T *v1, T *v2, T *vout) {
  unsigned int idx = blockIdx.x * blockDim.x + threadIdx.x;
  if (idx < n)
    atomicAdd(vout, v1[idx] * v2[idx]);
}

template <typename T, unsigned int THREADS_PER_ROW>
__global__ void multiply_mv(unsigned int size, unsigned int *row_pointers,
                            unsigned int *columns, T *data, T *vector,
                            T *result, T alpha, T beta) {
  // heavily inspired from
  // https://gpuopen.com/learn/amd-lab-notes/amd-lab-notes-spmv-docs-spmv_part1/
  unsigned int row = threadIdx.y + blockDim.y * blockIdx.x;
  if (row < size) {
    // determine the start and ends of each row
    unsigned int p = row_pointers[row];
    unsigned int q = row_pointers[row + 1];

    // start the sparse row * vector dot product operation
    double sum = 0;
    for (unsigned int i = p + threadIdx.x; i < q; i += THREADS_PER_ROW) {
      sum += data[i] * vector[columns[i]];
    }

    // finish the sparse row * vector dot product operation
#pragma unroll
    for (int i = THREADS_PER_ROW >> 1; i > 0; i >>= 1)
      sum += __shfl_down_sync(sum, i, THREADS_PER_ROW);

    // write to memory
    if (!threadIdx.x) {
      if (beta == 0) {
        result[row] = alpha * sum;
      } else {
        result[row] = alpha * sum + beta * result[row];
      }
    }
  }
}

__host__ __device__ inline unsigned int prevPowerOf2(unsigned int x) {
  if (x == 0)
    return 0;
  x |= (x >> 1);
  x |= (x >> 2);
  x |= (x >> 4);
  x |= (x >> 8);
  x |= (x >> 16);
  return x - (x >> 1);
}

double vectorNorm(int size, double *A) {
  int i;
  double norm = 0;
  for (i = 0; i < size; i++)
    norm += A[i] * A[i];
  norm = sqrt(norm);
  return norm;
}

namespace cuda {
template <class T> T *allocator<T>::allocate(std::size_t n) {
  T *x;
  cudaMallocManaged(&x, n * sizeof(T));
  return x;
}

template <class T> void allocator<T>::deallocate(T *p, std::size_t n) {
  cudaFree(p);
}

// Explicit instantiations:
template class allocator<int>;
template class allocator<unsigned int>;
template class allocator<double>;

template <typename T>
void add_vectors(unsigned int n, T *v1, T *v2, T c1, T c2) {
  int number_of_blocks = n / THREADS_PER_BLOCK + 1;
  add_vv<T><<<number_of_blocks, THREADS_PER_BLOCK>>>(n, v1, v2, c1, c2);
};
// Explicit instantiations:
template void add_vectors<int>(unsigned int n, int *v1, int *v22, int c1,
                               int c2);
template void add_vectors<double>(unsigned int n, double *v1, double *v2,
                                  double c1, double c2);

template <typename T> T multiply_vectors(unsigned int n, T *v1, T *v2) {
  int number_of_blocs = n / THREADS_PER_BLOCK + 1;
  T *result_shared;
  cudaMalloc(&result_shared, sizeof(T));
  cudaMemset(result_shared, 0, sizeof(T));
  multiply_vv<T>
      <<<number_of_blocs, THREADS_PER_BLOCK>>>(n, v1, v2, result_shared);
  T result;
  cudaMemcpy(&result, result_shared, sizeof(T), cudaMemcpyDeviceToHost);
  cudaFree(result_shared);
  return result;
};

// Explicit instantiations:
template int multiply_vectors<int>(unsigned int n, int *v1, int *v2);
template double multiply_vectors<double>(unsigned int n, double *v1,
                                         double *v2);

template <typename T>
void multiply_matrix_vector(unsigned int size, unsigned int *columns,
                            unsigned int *row_pointers, T *data, T *vector,
                            T *result, T alpha, T beta) {

  // Code tombe du camion
  // https://gpuopen.com/learn/amd-lab-notes/amd-lab-notes-spmv-docs-spmv_part1/
  int nnz_per_row = row_pointers[size] / size;
  int threads_per_row = prevPowerOf2(nnz_per_row);
  // limit the number of threads per row to be no larger than the wavefront
  // (warp) size; Couldn't find an exact figure for warp size let's hope 32 is
  // good enough
  threads_per_row = threads_per_row > 32 ? 32 : threads_per_row;
  int rows_per_block = THREADS_PER_BLOCK / threads_per_row;
  int num_blocks = (size + rows_per_block - 1) / rows_per_block;

  dim3 grid(num_blocks, 1, 1);
  dim3 block(threads_per_row, rows_per_block, 1);
  if (threads_per_row <= 2)
    multiply_mv<T, 2><<<grid, block>>>(size, row_pointers, columns, data,
                                       vector, result, alpha, beta);
  else if (threads_per_row <= 4)
    multiply_mv<T, 4><<<grid, block>>>(size, row_pointers, columns, data,
                                       vector, result, alpha, beta);
  else if (threads_per_row <= 8)
    multiply_mv<T, 8><<<grid, block>>>(size, row_pointers, columns, data,
                                       vector, result, alpha, beta);
  else if (threads_per_row <= 16)
    multiply_mv<T, 16><<<grid, block>>>(size, row_pointers, columns, data,
                                        vector, result, alpha, beta);
  else if (threads_per_row <= 32)
    multiply_mv<T, 32><<<grid, block>>>(size, row_pointers, columns, data,
                                        vector, result, alpha, beta);
  else
    multiply_mv<T, 64><<<grid, block>>>(size, row_pointers, columns, data,
                                        vector, result, alpha, beta);
};

// Explicit instantiations:
template void multiply_matrix_vector<int>(unsigned int size,
                                          unsigned int *columns,
                                          unsigned int *row_pointers, int *data,
                                          int *vector, int *result, int alpha,
                                          int beta);
template void multiply_matrix_vector<double>(
    unsigned int size, unsigned int *columns, unsigned int *row_pointers,
    double *data, double *vector, double *result, double alpha, double beta);

void conjugate_gradient(unsigned int size, unsigned int *columns,
                        unsigned int *row_pointers, double *data, double *B) {
  double *Residual, *search_direction, *A_search_direction;
  cudaMallocManaged(&Residual, size * sizeof(double));
  cudaMallocManaged(&search_direction, size * sizeof(double));
  cudaMallocManaged(&A_search_direction, size * sizeof(double));

  // Initialize residual vector
  multiply_matrix_vector<double>(size, columns, row_pointers, data, B, Residual,
                                 1.0, 1.0);
  add_vectors<double>(size, B, Residual, 1.0, -1.0);
  // Skip all if utility is right guess
  if (vectorNorm(size, Residual) < TOL) {
    cudaFree(Residual);
    cudaFree(search_direction);
    cudaFree(A_search_direction);
  }

  // Initialize search direction vector
  memcpy(search_direction, Residual, sizeof(double) * size);
  // Iterate until convergence
  for (int i = 0; i < 100; i++) {
    multiply_matrix_vector<double>(size, columns, row_pointers, data,
                                   search_direction, A_search_direction, 1, 1);

    double old_rz_product = multiply_vectors<double>(size, Residual, Residual);
    double step_size =
        old_rz_product /
        multiply_vectors<double>(size, search_direction, A_search_direction);

    // Update solution
    // dotScalarVector(size, step_size, search_direction, utility);
    // addVector(size, B, utility, B);
    add_vectors<double>(size, search_direction, B, step_size, 1.0);

    // Update residual
    // dotScalarVector(size, step_size, A_search_direction, utility);
    // substracVector(size, Residual, utility, Residual);
    add_vectors<double>(size, A_search_direction, Residual, -step_size, 1.0);

    if (vectorNorm(size, Residual) < TOL) {
      break;
    }

    // Update search direction
    double beta =
        multiply_vectors<double>(size, Residual, Residual) / old_rz_product;
    // dotScalarVector(size, beta, search_direction, utility);
    // addVector(size, Residual, utility, search_direction);
    add_vectors(size, Residual, search_direction, 1.0, beta);
  }
  cudaFree(Residual);
  cudaFree(search_direction);
  cudaFree(A_search_direction);
}

void test() {
  bool good = true;
  int *x, *y, *z;

  // Allocate Unified Memory – accessible from CPU or GPU
  cudaMallocManaged(&x, 1000 * sizeof(int));
  cudaMallocManaged(&y, 1000 * sizeof(int));
  z = (int *)malloc(1000 * sizeof(int));

  for (int i = 0; i < 1000; i++) {
    x[i] = std::rand(); // generate a random number
    y[i] = std::rand();
  }
  memcpy(z, y, 1000 * sizeof(int));

  add_vectors<int>(1000, x, y);

  for (int i = 0; i < 1000; i++) {
    z[i] = x[i] + z[i];
  }
  cudaDeviceSynchronize();
  for (int i = 0; i < 1000; i++) {
    if (z[i] != y[i]) {
      good = false;
      std::cout << "z " << i << " : " << z[i] << " y : " << y[i] << "\n";
    }
  }
  int r1, r2 = 0;

  r1 = multiply_vectors(1000, x, y);
  for (int i = 0; i < 1000; i++) {
    r2 += x[i] * y[i];
  }

  cudaDeviceSynchronize();
  if (r1 != r2) {
    good = false;
    std::cout << "r1 : " << r1 << " " << "r2 : " << r2 << "\n";
  }

  if (good) {
    std::cout << "YAY\n";
  } else {
    std::cout << "NAY\n";
  }
  cudaFree(x);
  cudaFree(y);
  free(z);
}

} // namespace cuda
