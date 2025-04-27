#pragma once
#include <cstddef>
namespace cuda {
// an allocator using unified memory
template <class T> class allocator {
public:
  typedef T value_type;
  allocator() noexcept {}
  template <class U> allocator(const allocator<U> &) noexcept {}
  T *allocate(std::size_t n);
  void deallocate(T *p, std::size_t n);
};

template <class T, class U>
constexpr bool operator==(const allocator<T> &, const allocator<U> &) noexcept {
  return true;
}

template <class T, class U>
constexpr bool operator!=(const allocator<T> &, const allocator<U> &) noexcept {
  return false;
}
template <typename T>
void add_vectors(unsigned int n, T *v1, T *v2, T c1 = 1, T c2 = 1);
template <typename T> T multiply_vectors(unsigned int n, T *v1, T *v2);
template <typename T>
void multiply_matrix_vector(unsigned int size, unsigned int columns,
                            unsigned int row_pointers, T *data, T *vector,
                            T *result);
void conjugate_gradient(unsigned int size, unsigned int *columns,
                        unsigned int *row_pointers, double *data, double *B);
void test();
} // namespace cuda
