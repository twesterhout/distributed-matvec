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

unsigned plugin_get_number_spins();
int plugin_get_hamming_weight();
int plugin_get_spin_inversion();
int plugin_is_representative(uint64_t count, uint64_t const spins[],
                             uint8_t output[]);
