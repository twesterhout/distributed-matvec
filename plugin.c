#include "plugin.h"
#include <lattice_symmetries/lattice_symmetries.h>

#include <complex.h>
#include <stdio.h>
#include <stdlib.h>

typedef struct plugin_state_type {
  ls_spin_basis *basis;
  ls_operator *matrix;
} plugin_state_type;

plugin_state_type plugin_state;

static void print_error_message_and_exit(ls_error_code const status) {
  fprintf(stderr, "Error (%i): %s\n", status, ls_error_to_string(status));
  exit(1);
}

void plugin_init(unsigned const number_spins) {
  int const hamming_weight = number_spins / 2;

  ls_group *group;
  ls_error_code status = ls_create_group(&group, 0, NULL);
  if (status != LS_SUCCESS) {
    goto fail1;
  }

  ls_spin_basis *basis;
  status = ls_create_spin_basis(&basis, group, number_spins, hamming_weight,
                                /*spin_inversion=*/0);
  if (status != LS_SUCCESS) {
    goto fail2;
  }
  status = ls_build(basis);
  if (status != LS_SUCCESS) {
    goto fail3;
  }

  // Get number of representatives
  uint64_t number_states;
  status = ls_get_number_states(basis, &number_states);
  if (status != LS_SUCCESS) {
    goto fail4;
  }
  printf("Hilbert space dimension is %zu\n", (size_t)number_states);

  // Heisenberg Hamiltonian
  // clang-format off
  _Complex double const matrix[] = {1.0,  0.0,  0.0, 0.0,
                                    0.0, -1.0,  2.0, 0.0,
                                    0.0,  2.0, -1.0, 0.0,
                                    0.0,  0.0,  0.0, 1.0};
  // clang-format on
  uint16_t(*edges)[2] = malloc(number_spins * sizeof(uint16_t[2]));
  if (edges == NULL) {
    status = LS_OUT_OF_MEMORY;
    goto fail5;
  }
  for (unsigned i = 0U; i < number_spins; ++i) {
    edges[i][0] = i;
    edges[i][1] = (i + 1) % number_spins;
  }
  ls_interaction *term;
  status = ls_create_interaction2(&term, matrix, number_spins, edges);
  if (status != LS_SUCCESS) {
    goto fail6;
  }
  ls_operator *hamiltonian;
  status = ls_create_operator(&hamiltonian, basis, 1,
                              (ls_interaction const *const *)&term);
  if (status != LS_SUCCESS) {
    goto fail7;
  }

fail7:
  ls_destroy_interaction(term);
fail6:
  free(edges);
fail5:
fail4:
fail3:
  if (status != LS_SUCCESS) {
    ls_destroy_spin_basis(basis);
  }
fail2:
  ls_destroy_group(group);
fail1:
  if (status != LS_SUCCESS) {
    print_error_message_and_exit(status);
  }

  plugin_state = (plugin_state_type){basis, hamiltonian};
}

void plugin_deinit() {
  ls_destroy_spin_basis(plugin_state.basis);
  ls_destroy_operator(plugin_state.matrix);
}

uint64_t plugin_get_dimension() {
  uint64_t dim;
  ls_error_code status = ls_get_number_states(plugin_state.basis, &dim);
  if (status != LS_SUCCESS) {
    print_error_message_and_exit(status);
  }
  return dim;
}

uint64_t plugin_get_max_nonzero_per_row() {
  return ls_operator_max_buffer_size(plugin_state.matrix);
}

uint64_t const *plugin_get_basis_states() {
  ls_states *states;
  ls_get_states(&states, plugin_state.basis);
  uint64_t const *p = ls_states_get_data(states);
  ls_destroy_states(states);
  return p;
}
