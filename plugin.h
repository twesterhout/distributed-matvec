#pragma once

#include <complex.h>
#include <stdint.h>

void plugin_init(unsigned number_spins);
void plugin_deinit();
uint64_t plugin_get_dimension();
uint64_t plugin_get_max_nonzero_per_row();
uint64_t const *plugin_get_basis_states();
uint64_t plugin_apply_operator(uint64_t bits, uint64_t *other_spins,
                               _Complex double *other_coeffs);
uint64_t plugin_get_index(uint64_t bits);
void plugin_matvec(double const *x, double *y);
