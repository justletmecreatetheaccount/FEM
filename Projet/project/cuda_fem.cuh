#pragma once
#include "fem.hpp"

namespace cuda {
namespace fem {
void assemble_system(::fem::Problem &problem);
}
} // namespace cuda
