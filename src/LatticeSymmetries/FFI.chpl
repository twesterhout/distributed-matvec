module FFI {
  // use CTypes;
  use CPtr;
  use SysCTypes;

  pragma "no doc"
  pragma "fn synchronization free"
  extern proc c_pointer_return(const ref x : ?t) : c_ptr(t);

  inline proc c_const_ptrTo(const ref arr: []) {
    if (!arr.isRectangular() || !arr.domain.dist._value.dsiIsLayout()) then
      compilerError("Only single-locale rectangular arrays support c_ptrTo() at present");

    if (arr._value.locale != here) then
      halt("c_ptrTo() can only be applied to an array from the locale on which it lives (array is on locale "
           + arr._value.locale.id:string + ", call was made on locale " + here.id:string + ")");
    return c_pointer_return(arr[arr.domain.low]);
  }
  inline proc c_const_ptrTo(const ref x) {
    return c_pointer_return(x);
  }

  require "lattice_symmetries_haskell.h";

  extern type ls_hs_particle_type = c_int;
  extern const LS_HS_SPIN : ls_hs_particle_type;
  extern const LS_HS_SPINFUL_FERMION : ls_hs_particle_type;
  extern const LS_HS_SPINLESS_FERMION : ls_hs_particle_type;

  extern record ls_hs_basis_kernels {
    var state_index_kernel : c_fn_ptr;
    var state_index_data : c_void_ptr;
  }
  extern record ls_hs_basis {
    var number_sites : c_int;
    var number_particles : c_int;
    var number_up : c_int;
    var particle_type : ls_hs_particle_type;
    var state_index_is_identity : bool;
    var requires_projection : bool;
    var kernels : c_ptr(ls_hs_basis_kernels);
    var representatives : chpl_external_array;
    // ... other stuff ...
  }
  extern record ls_hs_nonbranching_terms {
    var number_terms : c_int;
    var number_bits : c_int;
    // ... other stuff ...
  }
  extern record ls_hs_operator {
    var basis : c_ptr(ls_hs_basis);
    var off_diag_terms : c_ptr(ls_hs_nonbranching_terms);
    var diag_terms : c_ptr(ls_hs_nonbranching_terms);
    // ... other stuff ...
  }

  extern proc ls_hs_init();
  extern proc ls_hs_exit();

  extern proc ls_hs_create_basis(particleType : ls_hs_particle_type, numberSites : c_int,
                                 numberParticles : c_int, numberUp : c_int) : c_ptr(ls_hs_basis);
  extern proc ls_hs_clone_basis(basis : c_ptr(ls_hs_basis)) : c_ptr(ls_hs_basis);
  extern proc ls_hs_destroy_basis_v2(basis : c_ptr(ls_hs_basis));
  extern proc ls_hs_min_state_estimate(basis : c_ptr(ls_hs_basis)) : uint(64);
  extern proc ls_hs_max_state_estimate(basis : c_ptr(ls_hs_basis)) : uint(64);

  extern proc ls_hs_create_spin_basis_from_json(json_string : c_string) : c_ptr(ls_hs_basis);
  extern proc ls_hs_create_spin_basis_from_yaml(yaml_filename : c_string) : c_ptr(ls_hs_basis);
  extern proc ls_hs_create_spinful_fermion_basis_from_json(json_string : c_string) : c_ptr(ls_hs_basis);
  extern proc ls_hs_create_spinless_fermion_basis_from_json(json_string : c_string) : c_ptr(ls_hs_basis);

  extern proc ls_hs_basis_build(basis : c_ptr(ls_hs_basis));

  extern proc ls_hs_state_index(basis : c_ptr(ls_hs_basis), batch_size : c_ptrdiff,
    spins : c_ptr(uint(64)), spins_stride : c_ptrdiff, indices : c_ptr(c_ptrdiff),
    indices_stride : c_ptrdiff);

  extern proc ls_hs_is_representative(basis : c_ptr(ls_hs_basis), batch_size : c_ptrdiff,
    alphas : c_ptr(uint(64)), alphas_stride : c_ptrdiff, are_representatives : c_ptr(uint(8)),
    norms : c_ptr(real(64)));

  extern proc ls_hs_state_info(basis : c_ptr(ls_hs_basis), batch_size : c_ptrdiff,
                               alphas : c_ptr(uint(64)), alphas_stride : c_ptrdiff,
                               betas : c_ptr(uint(64)), betas_stride : c_ptrdiff,
                               characters : c_ptr(complex(128)), norms : c_ptr(real(64)));

  // extern proc ls_hs_create_basis_kernels(basis : c_ptr(ls_hs_basis)) : c_ptr(ls_hs_basis_kernels);
  // extern proc ls_hs_destroy_basis_kernels(kernels : c_ptr(ls_hs_basis_kernels));

  extern proc ls_hs_create_operator(basis : c_ptr(ls_hs_basis),
    s : c_string, numberTuples : c_int, tupleSize : c_int, tuples : c_ptr(c_int)) : c_ptr(ls_hs_operator);
  extern proc ls_hs_operator_plus(a : c_ptr(ls_hs_operator), b : c_ptr(ls_hs_operator)) : c_ptr(ls_hs_operator);
  extern proc ls_hs_print_terms(op : c_ptr(ls_hs_operator));
  extern proc ls_hs_destroy_operator_v2(op : c_ptr(ls_hs_operator));

  extern proc ls_hs_load_hamiltonian_from_yaml(filename : c_string) : c_ptr(ls_hs_operator);

  extern proc ls_hs_operator_apply_diag_kernel(op : c_ptr(ls_hs_operator),
    batchSize : c_ptrdiff, alphas : c_ptr(uint(64)), alphas_stride : c_ptrdiff,
    coeffs : c_ptr(complex(128)));

  extern proc ls_hs_operator_apply_off_diag_kernel(op : c_ptr(ls_hs_operator),
    batchSize : c_ptrdiff, alphas : c_ptr(uint(64)), alphas_stride : c_ptrdiff,
    betas : c_ptr(uint(64)), betas_stride : c_ptrdiff, coeffs : c_ptr(complex(128)));

  extern proc ls_hs_evaluate_wavefunction_via_statevector(
    kernels : c_ptr(ls_hs_basis_kernels), batch_size : c_ptrdiff,
    alphas : c_ptr(uint(64)), alphas_stride : c_ptrdiff, state_vector : c_void_ptr,
    element_size : uint(64), coeffs : c_void_ptr);

  extern proc ls_hs_basis_has_fixed_hamming_weight(basis : c_ptr(ls_hs_basis)) : bool;

  extern record ls_chpl_kernels {
    var enumerate_states : c_fn_ptr;
  }
  extern proc ls_hs_internal_set_chpl_kernels(kernels : c_ptr(ls_chpl_kernels));


  require "helper.h";
  extern proc print_external_string(s : c_string);


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
}
