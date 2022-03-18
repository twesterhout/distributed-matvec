module ApplyOperator {
  use CPtr;
  use DynamicIters;
  use RangeChunk;
  use SysCTypes;
  use Time;

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
    var kernels : c_ptr(ls_hs_basis_kernels);
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
  extern proc ls_hs_destroy_basis_v2(basis : c_ptr(ls_hs_basis));
  extern proc ls_hs_min_state_estimate(basis : c_ptr(ls_hs_basis)) : uint(64);
  extern proc ls_hs_max_state_estimate(basis : c_ptr(ls_hs_basis)) : uint(64);

  extern proc ls_hs_perform_gc();

  extern proc ls_hs_state_index(basis : c_ptr(ls_hs_basis), batch_size : c_ptrdiff,
    spins : c_ptr(uint(64)), spins_stride : c_ptrdiff, indices : c_ptr(c_ptrdiff),
    indices_stride : c_ptrdiff);

  // extern proc ls_hs_create_basis_kernels(basis : c_ptr(ls_hs_basis)) : c_ptr(ls_hs_basis_kernels);
  // extern proc ls_hs_destroy_basis_kernels(kernels : c_ptr(ls_hs_basis_kernels));

  extern proc ls_hs_create_operator(basis : c_ptr(ls_hs_basis),
    s : c_string, numberTuples : c_int, tupleSize : c_int, tuples : c_ptr(c_int)) : c_ptr(ls_hs_operator);
  extern proc ls_hs_operator_plus(a : c_ptr(ls_hs_operator), b : c_ptr(ls_hs_operator)) : c_ptr(ls_hs_operator);
  extern proc ls_hs_print_terms(op : c_ptr(ls_hs_operator));
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

  extern proc ls_hs_basis_has_fixed_hamming_weight(basis : c_ptr(ls_hs_basis)) : bool;

  extern record ls_chpl_kernels {
    var enumerate_states : c_fn_ptr;
  }
  extern proc ls_hs_internal_set_chpl_kernels(kernels : c_ptr(ls_chpl_kernels));


  class Basis {
    var payload : c_ptr(ls_hs_basis);
    var kernels : c_ptr(ls_hs_basis_kernels);
    var owning : bool;

    proc init(p : c_ptr(ls_hs_basis), owning : bool = true) {
      this.payload = p;
      this.kernels = this.payload.deref().kernels;
      this.owning = owning;
    }

    proc deinit() {
      if (owning) then
        ls_hs_destroy_basis_v2(payload);
    }

    proc isSpinBasis() { return payload.deref().particle_type == LS_HS_SPIN; }
    proc isSpinfulFermionicBasis() { return payload.deref().particle_type == LS_HS_SPINFUL_FERMION; }
    proc isSpinlessFermionicBasis() { return payload.deref().particle_type == LS_HS_SPINLESS_FERMION; }
    proc isStateIndexIdentity() { return payload.deref().state_index_is_identity; }
    // proc hasFixedHammingWeight() { return ls_hs_basis_has_fixed_hamming_weight(payload); }

    proc numberSites() : int { return payload.deref().number_sites; }
    proc numberParticles() : int { return payload.deref().number_particles; }
    proc numberUp() : int { return payload.deref().number_up; }

    proc minStateEstimate() : uint(64) { return ls_hs_min_state_estimate(payload); }
    proc maxStateEstimate() : uint(64) { return ls_hs_max_state_estimate(payload); }
  }
  proc SpinBasis(numberSites : int, hammingWeight : int = -1) {
    return new Basis(ls_hs_create_basis(LS_HS_SPIN, numberSites:c_int, 
                                        numberSites:c_int, hammingWeight:c_int));
  }
  proc SpinlessFermionicBasis(numberSites : int, numberParticles : int = -1) {
    return new Basis(ls_hs_create_basis(LS_HS_SPINLESS_FERMION, numberSites:c_int,
                                        numberParticles:c_int, -1));
  }
  proc SpinfulFermionicBasis(numberSites : int, numberParticles : int = -1) {
    return new Basis(ls_hs_create_basis(LS_HS_SPINFUL_FERMION, numberSites:c_int,
                                        numberParticles:c_int, -1));
  }
  proc SpinfulFermionicBasis(numberSites : int, numberUp : int, numberDown : int) {
    return new Basis(ls_hs_create_basis(LS_HS_SPINFUL_FERMION, numberSites:c_int,
                                        (numberUp + numberDown):c_int, numberUp:c_int));
  }

  proc localCompress(numberTerms : int, ref betas : [?D] uint(64),
                     ref coeffs : [?D2] ?eltType, ref offsets : [?D3] int)
      where D.rank == 1 && D2.rank == 1 && D3.rank == 1 {
    assert(betas.size == coeffs.size);
    const batchSize = offsets.size - 1;
    assert(betas.size == numberTerms * batchSize);

    var offset : int = 0;
    var i : int = 0;

    offsets[0] = 0;
    for batchIdx in 0 ..# batchSize {
      for termIdx in 0 ..# numberTerms {
        if (coeffs[i] != 0) {
          if (i > offset) {
            coeffs[offset] = coeffs[i];
            betas[offset] = betas[i];
          }
          offset += 1;
        }
        i += 1;
      }
      offsets[batchIdx + 1] = offset;
    }
  }

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

  class Operator {
    var payload : c_ptr(ls_hs_operator);
    var basis : c_ptr(ls_hs_basis);
    var kernels : c_ptr(ls_hs_basis_kernels);

    proc init(const ref basis : Basis, expression : string, const ref indices : [?D] ?i)
      where D.rank == 2 {

      const c_indices = indices:c_int;
      this.payload = ls_hs_create_operator(basis.payload, expression.localize().c_str(),
        indices.dim(0).size:c_int, indices.dim(1).size:c_int, c_const_ptrTo(c_indices));
      this.basis = this.payload.deref().basis;
      this.kernels = this.basis.deref().kernels;
    }
    proc init(raw : c_ptr(ls_hs_operator)) {
      this.payload = raw;
      this.basis = this.payload.deref().basis;
      this.kernels = this.basis.deref().kernels;
    }
    proc deinit() { ls_hs_destroy_operator_v2(payload); }

    inline proc numberOffDiagTerms() : int {
      const p = payload.deref().off_diag_terms;
      if (p == nil) { return 0; }
      return p.deref().number_terms:int;
    }

    proc writeTerms() {
      ls_hs_print_terms(payload);
    }

    proc _localApplyDiag(const ref alphas : [?D] uint(64), ref coeffs : [?D2] complex(128))
        where D.rank == 1 && D2.rank == 1 {
      assert(alphas.size == coeffs.size);
      const batchSize = alphas.size;
      ls_hs_operator_apply_diag_kernel(payload, batchSize, c_const_ptrTo(alphas), 1, c_ptrTo(coeffs));
    }

    proc _localApplyOffDiag(const ref alphas : [?D] uint(64), ref betas : [?D2] uint(64),
                            ref coeffs : [?D3] complex(128))
        where D.rank == 1 && D2.rank == 1 && D3.rank == 1 {
      const batchSize = alphas.size;
      const numberTerms = numberOffDiagTerms();
      assert(betas.size == batchSize * numberTerms);
      assert(coeffs.size == betas.size);
      ls_hs_operator_apply_off_diag_kernel(payload, batchSize,
          c_const_ptrTo(alphas), 1, c_ptrTo(betas), 1, c_ptrTo(coeffs));
    }

  }

  operator +(const ref a : Operator, const ref b : Operator) {
    return new Operator(ls_hs_operator_plus(a.payload, b.payload));
  }

  proc localMatVecPart(const ref matrix : Operator,
                       const ref alphas : [?D] uint(64),
                       const ref elements : [],
                       const ref stateVector : [?D2] ?eltType,
                       ref dest : [?D3] ?eltType2)
      where D.rank == 1 && D2.rank == 1
              && D3.rank == D2.rank && eltType2 == eltType {
    const batchSize = alphas.size;
    const numberTerms = matrix.numberOffDiagTerms();

    var ds : [0 ..# batchSize] complex(128) = noinit;
    matrix._localApplyDiag(alphas, ds);
    dest = (ds * elements):eltType;

    var betas : [0 ..# (batchSize * numberTerms)] uint(64) = noinit;
    var cs : [0 ..# (batchSize * numberTerms)] complex(128) = noinit;
    var xs : [0 ..# (batchSize * numberTerms)] eltType = noinit;
    var offsets : [0 ..# (batchSize + 1)] int = noinit;
    matrix._localApplyOffDiag(alphas, betas, cs);
    localCompress(numberTerms, betas, cs, offsets);
    localEvaluateWavefunction(betas, stateVector, xs, matrix.kernels);
    for i in 0 ..# batchSize {
      var acc : complex(128) = 0;
      for k in offsets[i] .. offsets[i + 1] - 1 {
        acc += conjg(cs[k]) * xs[k];
      }
      dest[dest.domain.low + i] += acc:eltType;
    }
  }

  config const localMatVecChunkSize = 10;

  proc localMatVec(const ref matrix : Operator,
                   const ref representatives : [] uint(64),
                   const ref xs : [] ?eltType,
                   ref ys : [] ?eltType2)
      where eltType == eltType2 {
   
    const numChunks = (representatives.size + localMatVecChunkSize - 1)
                      / localMatVecChunkSize;
    // TODO: are we going to run into trouble that there are too many tasks?
    coforall indices in chunks(0 ..# representatives.size, numChunks) {
      localMatVecPart(matrix, representatives[indices],
                      xs[indices], xs, ys[indices]);
    }
  }
  proc localMatVec(const ref matrix : Operator,
                   const ref representatives : [] uint(64),
                   const ref xs : [] ?eltType) {
    var ys : [xs.domain] eltType = noinit;
    localMatVec(matrix, representatives, xs, ys);
    return ys;
  }

  export proc bar(): int {
    return 0;
  }
}

