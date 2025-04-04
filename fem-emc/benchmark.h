#pragma once

#ifdef __cplusplus
extern "C" {
#endif

void elasticity_solve(const char *meshfile, const char *outfile, double E, double nu, double rho, double g);
void elasticity_solve_csr(const char *meshfile, const char *outfile, double E, double nu, double rho, double g);

#ifdef __cplusplus
}
#endif // C++
