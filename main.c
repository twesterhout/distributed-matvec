#include "basis.h"
#include <lattice_symmetries/lattice_symmetries.h>
#include <stdlib.h>

static ls_spin_basis *make_linear_chain_basis(unsigned const number_spins) {
  LATTICE_SYMMETRIES_LOG_DEBUG("Creating basis for a 1D chain of %u spins...\n",
                               number_spins);
  LATTICE_SYMMETRIES_CHECK(number_spins % 2 == 0, "");
  int const hamming_weight = number_spins / 2;

  unsigned *permutation = malloc(number_spins * sizeof(unsigned));
  LATTICE_SYMMETRIES_CHECK(permutation != NULL, "");

  for (unsigned i = 0; i < number_spins; ++i) {
    permutation[i] = (i + 1) % number_spins;
  }
  unsigned sector = (number_spins % 4 == 2) ? number_spins / 2 : 0;
  ls_symmetry *translation;
  ls_error_code status =
      ls_create_symmetry(&translation, number_spins, permutation, sector);
  LATTICE_SYMMETRIES_CHECK(status == LS_SUCCESS, "");

  for (unsigned i = 0; i < number_spins; ++i) {
    permutation[i] = (number_spins - 1 - i);
  }
  sector = (number_spins % 4 == 2) ? 1 : 0;
  ls_symmetry *parity;
  status = ls_create_symmetry(&parity, number_spins, permutation, sector);
  LATTICE_SYMMETRIES_CHECK(status == LS_SUCCESS, "");

  ls_symmetry const *generators[] = {translation, parity};
  ls_group *group;
  status = ls_create_group(
      &group, sizeof(generators) / sizeof(ls_symmetry const *), generators);
  LATTICE_SYMMETRIES_CHECK(status == LS_SUCCESS, "");

  ls_destroy_symmetry(translation);
  ls_destroy_symmetry(parity);

  ls_spin_basis *basis;
  status =
      ls_create_spin_basis(&basis, group, number_spins, hamming_weight,
                           /*spin_inversion=*/(number_spins % 4 == 2) ? 1 : 0);
  LATTICE_SYMMETRIES_CHECK(status == LS_SUCCESS, "");

  ls_destroy_group(group);
  return basis;
}

// static ls_dis

int main() {
  ls_enable_logging();
  ls_spin_basis *basis = make_linear_chain_basis(10);

  ls_destroy_spin_basis(basis);
  return 0;
}
