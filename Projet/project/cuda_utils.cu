#include "cuda_utils.cuh"
#include "defines.hpp"
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <vector>

__global__ void testing_cuda_vec(double *v) {
  for (int i = 0; i < 6; i++) {
    v[i] = 69;
  }
}

// Kernel function to add the elements of two arrays
template <typename T>
__global__ void add_vv(unsigned int n, T *v1, T *v2, T c1, T c2) {
  unsigned int idx = blockIdx.x * blockDim.x + threadIdx.x;
  if (idx < n)
    v2[idx] = c1 * v1[idx] + c2 * v2[idx];
}

// Kernel to mult two vectors (ettonamment plus rapide que juste le faire sur un
// thread)
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
      sum += __shfl_down_sync(0xffffffff, sum, i, THREADS_PER_ROW);

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
  x |= (x >> 1);
  x |= (x >> 2);
  x |= (x >> 4);
  x |= (x >> 8);
  x |= (x >> 16);
  return x - (x >> 1);
}

double vectorNorm(unsigned int size, double *A) {
  double norm = 0;
  for (unsigned int i = 0; i < size; i++)
    norm += A[i] * A[i];
  norm = sqrt(norm);
  return norm;
}

/* Definition d'un allocateur pour avoir des vecteur avec de la memoire managee
 * cuda
 *
 */
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

/* Additionne deux vecteurs
 * ---
 * v1, v2 > les vecteurs, le resultat de l'operation est mis dans v2
 * c1, c2 > constantes pour multiplier les vecteurs si envie
 */
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
  cudaDeviceSynchronize();
  T result;
  cudaMemcpy(&result, result_shared, sizeof(T), cudaMemcpyDeviceToHost);
  cudaFree(result_shared);
  return result;
};

// Explicit instantiations:
template int multiply_vectors<int>(unsigned int n, int *v1, int *v2);
template double multiply_vectors<double>(unsigned int n, double *v1,
                                         double *v2);

/* Multiplies deux vecteurs (je conseille d'aller voir le lien pour comprendre
 * le code)
 * ---
 * size : la taille de la matrice et du vecteur
 * columns, row_pointers, data : la matrice sous format csr
 * vector : le vecteur avec lequel on multiplies
 * result : le vecteur dans lequel on met le resultat
 * alpha : constante pour multiplier le resultat
 * beta : constante pour preserver la valeur precedente dans result
 */
template <typename T>
void multiply_matrix_vector(unsigned int size, unsigned int *columns,
                            unsigned int *row_pointers, T *data, T *vector,
                            T *result, T alpha, T beta) {

  // Code tombe du camion
  // https://gpuopen.com/learn/amd-lab-notes/amd-lab-notes-spmv-docs-spmv_part1/
  unsigned int nnz_per_row = row_pointers[size] / size;
  unsigned int threads_per_row = prevPowerOf2(nnz_per_row);
  // limit the number of threads per row to be no larger than the wavefront
  // (warp) size; Couldn't find an exact figure for warp size let's hope 32 is
  // good enough
  threads_per_row = threads_per_row > 32 ? 32 : threads_per_row;
  unsigned int rows_per_block = THREADS_PER_BLOCK / threads_per_row;
  unsigned int num_blocks = (size + rows_per_block - 1) / rows_per_block;

  dim3 grid(num_blocks, 1, 1);
  dim3 block(threads_per_row, rows_per_block, 1);
  if (threads_per_row <= 2) {
    multiply_mv<T, 2><<<grid, block>>>(size, row_pointers, columns, data,
                                       vector, result, alpha, beta);
  } else if (threads_per_row <= 4) {
    multiply_mv<T, 4><<<grid, block>>>(size, row_pointers, columns, data,
                                       vector, result, alpha, beta);
  } else if (threads_per_row <= 8) {
    multiply_mv<T, 8><<<grid, block>>>(size, row_pointers, columns, data,
                                       vector, result, alpha, beta);
  } else if (threads_per_row <= 16) {
    multiply_mv<T, 16><<<grid, block>>>(size, row_pointers, columns, data,
                                        vector, result, alpha, beta);
  } else if (threads_per_row <= 32) {
    multiply_mv<T, 32><<<grid, block>>>(size, row_pointers, columns, data,
                                        vector, result, alpha, beta);
  } else {
    multiply_mv<T, 64><<<grid, block>>>(size, row_pointers, columns, data,
                                        vector, result, alpha, beta);
  }
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

// la solution est enregistree dans B
void conjugate_gradient(unsigned int size, unsigned int *columns,
                        unsigned int *row_pointers, double *data, double *B) {
  double *Residual, *search_direction, *A_search_direction;
  cudaMallocManaged(&Residual, size * sizeof(double));
  cudaMallocManaged(&search_direction, size * sizeof(double));
  cudaMallocManaged(&A_search_direction, size * sizeof(double));

  // Initialize residual vector
  multiply_matrix_vector<double>(size, columns, row_pointers, data, B, Residual,
                                 1.0, 0.0);
  cudaDeviceSynchronize();
  add_vectors<double>(size, B, Residual, 1.0, -1.0);
  cudaDeviceSynchronize();
  // Skip all if B is right guess
  double old_res_norm = vectorNorm(size, Residual);
  if (old_res_norm < TOL) {
    cudaFree(Residual);
    cudaFree(search_direction);
    cudaFree(A_search_direction);
  }

  // Initialize search direction vector
  memcpy(search_direction, Residual, sizeof(double) * size);
  // Iterate until convergence
  for (int i = 0; i < 5; i++) {
    if (i == 4) {
      std::cout << "ouf\n";
    }
    multiply_matrix_vector<double>(size, columns, row_pointers, data,
                                   search_direction, A_search_direction, 1, 0);
    cudaDeviceSynchronize();

    double step_size =
        old_res_norm * old_res_norm /
        multiply_vectors<double>(size, search_direction, A_search_direction);

    // Update solution
    add_vectors<double>(size, search_direction, B, step_size, 1.0);
    cudaDeviceSynchronize();

    // Update residual
    add_vectors<double>(size, A_search_direction, Residual, -step_size, 1.0);
    cudaDeviceSynchronize();

    double new_res_norm = vectorNorm(size, Residual);

    if (new_res_norm < TOL) {
      break;
    }

    // Update search direction
    double beta = (new_res_norm * new_res_norm) / (old_res_norm * old_res_norm);
    add_vectors<double>(size, Residual, search_direction, 1, beta);
    cudaDeviceSynchronize();
    old_res_norm = new_res_norm;
  }
  cudaFree(Residual);
  cudaFree(search_direction);
  cudaFree(A_search_direction);
}

// fonction pour tester ce qu'il y a au dessus
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

  add_vectors<int>(1000, x, y, -2, 3);

  for (int i = 0; i < 1000; i++) {
    z[i] = -2 * x[i] + z[i] * 3;
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
  std::cout << "prev pow of 1 : " << prevPowerOf2(1) << "\n";
  std::cout << "prev pow of 0 : " << prevPowerOf2(0) << "\n";
  std::cout << "prev pow of 23 : " << prevPowerOf2(23) << "\n";
  std::cout << "prev pow of 65 : " << prevPowerOf2(65) << "\n";

  std::cout << "matrix mult \n";

  /*
  01 00 00 00 00 .6    69  =  82.2
  03 00 07 02 09 00    42  =  891
  50 00 65 02 .6 65    54  =  8410.4
  00 00 00 01 00 00    00  =  0
  02 52 00 65 00 00    34  =  2322
  00 64 00 74 00 00    22  =  2688
  */

  double _B[6] = {69, 42, 54, 0, 34, 22};
  double _data[17] = {1,   0.6, 3, 7, 2,  9,  50, 65, 2,
                      0.6, 65,  1, 2, 52, 65, 64, 74};
  unsigned int _row_ptr[7] = {0, 2, 6, 11, 12, 15, 17};
  unsigned int _col[17] = {0, 5, 0, 2, 3, 4, 0, 2, 3, 4, 5, 3, 0, 1, 3, 1, 3};
  double *result;
  double *B;
  double *data;
  unsigned int *row_ptr, *col;
  cudaMallocManaged(&result, sizeof(double) * 6);
  cudaMallocManaged(&B, sizeof(double) * 6);
  cudaMallocManaged(&data, sizeof(double) * 17);
  cudaMallocManaged(&row_ptr, sizeof(unsigned int) * 17);
  cudaMallocManaged(&col, sizeof(unsigned int) * 6);

  memcpy(B, _B, sizeof(double) * 6);
  memcpy(data, _data, sizeof(double) * 17);
  memcpy(row_ptr, _row_ptr, sizeof(unsigned int) * 7);
  memcpy(col, _col, sizeof(unsigned int) * 17);

  multiply_matrix_vector<double>(6, col, row_ptr, data, B, result, 1, 1);
  cudaDeviceSynchronize();

  for (int i = 0; i < 6; i++) {
    std::cout << "result : " << result[i] << "\n";
  }

  conjugate_gradient(6, col, row_ptr, data, result);
  for (int i = 0; i < 6; i++) {
    std::cout << "result conj : " << result[i] << "\n";
  }

  cudaFree(result);
  cudaFree(B);
  cudaFree(data);
  cudaFree(row_ptr);
  cudaFree(col);

  std::vector<double, cuda::allocator<double>> test_vec(6);

  testing_cuda_vec<<<1, 1>>>(test_vec.data());
  cudaDeviceSynchronize();
  for (int i = 0; i < 6; i++) {
    std::cout << "test vec " << i << " : " << test_vec[i] << "\n";
  }
}

} // namespace cuda
