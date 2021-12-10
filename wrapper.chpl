use CPtr;
use SysCTypes;

require "lattice_symmetries/lattice_symmetries.h";
extern type ls_spin_basis;
extern type ls_operator;
extern type ls_flat_spin_basis;
extern type ls_error_code;

extern proc ls_enable_logging();
extern proc ls_get_number_spins(basis: c_ptr(ls_spin_basis)): c_uint;
extern proc ls_get_hamming_weight(basis: c_ptr(ls_spin_basis)): c_int;
extern proc ls_get_spin_inversion(basis: c_ptr(ls_spin_basis)): c_int;
extern proc ls_is_representative(basis: c_ptr(ls_spin_basis), count: uint(64),
                                 bits: c_ptr(uint(64)), is_repr: c_ptr(uint(8)));
// extern proc ls_convert_to_flat_spin_basis(ptr: c_ptr(c_ptr(ls_flat_spin_basis)),
//                                           basis: c_ptr(ls_spin_basis)): ls_error_code;
// extern proc ls_destroy_flat_spin_basis(ptr: c_ptr(ls_flat_spin_basis));
// extern proc ls_get_buffer_size_for_flat_spin_basis(basis: c_ptr(ls_flat_spin_basis)): uint(64);
// extern proc ls_serialize_flat_spin_basis(basis: c_ptr(ls_flat_spin_basis),
//                                          buffer: [] c_char, size: uint(64)): ls_error_code;
// extern proc ls_deserialize_flat_spin_basis(ptr: c_ptr(c_ptr(ls_flat_spin_basis)),
//                                            buffer: [] c_char, size: uint(64)): ls_error_code ;
// extern proc ls_copy_spin_basis(basis: c_ptr(ls_spin_basis)): c_ptr(ls_spin_basis);

require "lattice_symmetries_haskell.h";
extern record ls_hs_spin_basis_v1 {
  var payload: c_ptr(ls_spin_basis);
  var context: c_void_ptr;
}

extern record ls_hs_operator_v1 {
  var payload: c_ptr(ls_operator);
  var context: c_void_ptr;
}

extern proc ls_hs_init();
extern proc ls_hs_exit();
extern proc ls_hs_basis_and_hamiltonian_from_yaml(path: c_string,
    basis: c_ptr(ls_hs_spin_basis_v1), hamiltonian: c_ptr(ls_hs_operator_v1));
extern proc ls_hs_destroy_spin_basis(basis: c_ptr(ls_hs_spin_basis_v1));
extern proc ls_hs_destroy_operator(basis: c_ptr(ls_hs_operator_v1));
extern proc ls_hs_hdf5_get_dataset_rank(path: c_string, dataset: c_string):c_uint;
extern proc ls_hs_hdf5_get_dataset_shape(path: c_string, dataset: c_string,
                                         shape: c_ptr(uint(64)));
extern proc ls_hs_hdf5_create_dataset_u64(path: c_string,
    dataset: c_string, dim: c_uint, shape: c_ptr(uint(64)));
extern proc ls_hs_hdf5_create_dataset_f64(path: c_string,
    dataset: c_string, dim: c_uint, shape: c_ptr(uint(64)));
extern proc ls_hs_hdf5_write_chunk_u64(path: c_string, dataset: c_string,
    dim: c_uint, offset: c_ptr(uint(64)), shape: c_ptr(uint(64)), data: c_ptr(uint(64)));
extern proc ls_hs_hdf5_write_chunk_f64(path: c_string, dataset: c_string,
    dim: c_uint, offset: c_ptr(uint(64)), shape: c_ptr(uint(64)), data: c_ptr(real(64)));
extern proc ls_hs_hdf5_read_chunk_u64(path: c_string, dataset: c_string,
    dim: c_uint, offset: c_ptr(uint(64)), shape: c_ptr(uint(64)), data: c_ptr(uint(64)));
extern proc ls_hs_hdf5_read_chunk_f64(path: c_string, dataset: c_string,
    dim: c_uint, offset: c_ptr(uint(64)), shape: c_ptr(uint(64)), data: c_ptr(real(64)));
