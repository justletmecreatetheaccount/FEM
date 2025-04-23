#pragma once

#include "fem.hpp"
namespace fem {
/* This function reads the in file and creates all the necessary structures to
 * fill the given problem
 * ======
 * IN : filename -> the name of the file
 *    : problem -> the problem that need to be filled
 */
void read_mesh_file(const char *filename, Problem &problem);
} // namespace fem
