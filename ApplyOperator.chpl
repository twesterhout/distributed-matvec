module ApplyOperator {
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

  require "lattice_symmetries_haskell.h";

  extern type ls_hs_basis;
  extern type ls_hs_basis_kernels;
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

  extern proc ls_hs_create_spin_basis(numberSites : c_int, hammingWeight : c_int) : c_ptr(ls_hs_basis);
  extern proc ls_hs_destroy_basis_v2(basis : c_ptr(ls_hs_basis));

  extern proc ls_hs_create_operator(basis : c_ptr(ls_hs_basis),
    s : c_string, numberTuples : c_int, tupleSize : c_int, tuples : c_ptr(c_int)) : c_ptr(ls_hs_operator);
  extern proc ls_hs_destroy_operator_v2(op : c_ptr(ls_hs_operator));

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

  record Basis {
    var payload : c_ptr(ls_hs_basis);

    proc init(p : c_ptr(ls_hs_basis)) {
      payload = p;
    }

    proc deinit() { ls_hs_destroy_basis_v2(payload); }
  }
  proc SpinBasis(numberSites : int, hammingWeight : int = -1) {
    return new Basis(ls_hs_create_spin_basis(numberSites:c_int, hammingWeight:c_int));
  }

  record Operator {
    var payload : c_ptr(ls_hs_operator);

    proc init(const ref basis : Basis, expression : string, const ref indices : [?D] ?i)
      where D.rank == 2 {

      const c_indices = indices:c_int;
      payload = ls_hs_create_operator(basis.payload, expression.localize().c_str(),
        indices.dim(0).size:c_int, indices.dim(1).size:c_int, c_const_ptrTo(c_indices));
    }
    proc deinit() { ls_hs_destroy_operator_v2(payload); }

    inline proc numberOffDiagTerms() : int {
      if (payload.off_diag_terms == nil) { return 0; }
      return payload.off_diag_terms.number_terms:int;
    }

    proc localApplyDiag(const ref alphas : [?D] uint(64), ref coeffs : [?D2] complex(128))
        where D.rank == 1 && D2.rank == 1 {
      assert(alphas.size == coeffs.size);
      const batchSize = alphas.size;
      ls_hs_operator_apply_diag_kernel(payload, batchSize, c_const_ptrTo(alphas), 1, c_ptrTo(coeffs));
    }

    proc localApplyOffDiag(const ref alphas : [?D] uint(64), ref betas : [?D2] uint(64),
                           ref coeffs : [?D3] complex(128))
        where D.rank == 1 && D2.rank == 1 && D3.rank == 1 {
      const batchSize = alphas.size;
      const numberTerms = numberOffDiagTerms();
      assert(betas.size == batchSize * numberOffDiagTerms);
      assert(coeffs.size == betas.size);
      ls_hs_operator_apply_off_diag_kernel(payload, batchSize,
          c_const_ptrTo(alphas), 1, c_ptrTo(betas), 1, c_ptrTo(coeffs));
    }
  }

  // proc localApplyDiagonal(const ref alphas : [?D] uint(64), ref coeffs : [?D2] complex(128),
  //                         op : c_ptr(ls_hs_operator))
  //     where D.rank == 1 && D2.rank == 1 {
  //   assert(alphas.size == coeffs.size);
  //   const batchSize = alphas.size;
  //   ls_hs_operator_apply_diag_kernel(

  // }

  proc localEvaluateWavefunction(const ref alphas : [] uint, const ref state_vector : [] ?eltType,
                                 ref coeffs : [] ?eltType2, kernels : c_ptr(ls_hs_basis_kernels))
      where eltType == eltType2 {
    assert(alphas.size == coeffs.size);
    const batchSize = alphas.size;
    ls_hs_evaluate_wavefunction_via_statevector(
      kernels,
      batchSize,
      c_const_ptrTo(alphas), 1,
      c_const_ptrTo(state_vector):c_void_ptr, c_sizeof(eltType),
      c_ptrTo(coeffs):c_void_ptr
    );
  }

  proc main() {
    ls_hs_init();

    var basis = SpinBasis(10); // ls_hs_create_spin_basis(10, -1);

    var tuples = reshape(
        [0, 1, 3, 2, 0, 2], {0 ..# 3, 0 ..# 2});
    var op = new Operator(basis, "σᶻ₀ σᶻ₁", tuples); // ls_hs_create_operator(basis, "σᶻ₀ σᶻ₁", 3, 2, c_ptrTo(tuples));

    // ls_hs_destroy_operator_v2(op);
    // ls_hs_destroy_basis_v2(basis);
    writeln("Hello world!");
  }

  proc deinit() {
    ls_hs_exit();
  }

  export proc bar(): int {
    main();
    return 0;
  }
}

