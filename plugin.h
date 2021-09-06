#pragma once

#include <stdint.h>

void plugin_init(unsigned const number_spins);
void plugin_deinit();
uint64_t plugin_get_dimension();
uint64_t plugin_get_max_nonzero_per_row();
uint64_t const *plugin_get_basis_states();
