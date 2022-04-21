module Types {
  // use CTypes;
  use CPtr;
  use SysCTypes;
  use FFI;

  record Basis {
    var payload : c_ptr(ls_hs_basis);
    var owning : bool;

    proc init(p : c_ptr(ls_hs_basis), owning : bool = true) {
      this.payload = p;
      this.owning = owning;
    }
    proc init=(const ref from : Basis) {
      assert(here == from.locale);
      this.payload = ls_hs_clone_basis(from.payload);
      this.owning = true;
    }

    proc _destroy() {
      if (owning) then
        ls_hs_destroy_basis_v2(payload);
    }

    proc deinit() {
      _destroy();
    }

    proc build() {
      writeln("calling ls_hs_basis_build ...");
      ls_hs_basis_build(payload); }

    proc isSpinBasis() { return payload.deref().particle_type == LS_HS_SPIN; }
    proc isSpinfulFermionicBasis() { return payload.deref().particle_type == LS_HS_SPINFUL_FERMION; }
    proc isSpinlessFermionicBasis() { return payload.deref().particle_type == LS_HS_SPINLESS_FERMION; }
    proc isStateIndexIdentity() { return payload.deref().state_index_is_identity; }
    proc requiresProjection() { return payload.deref().requires_projection; }
    proc isHammingWeightFixed() { return ls_hs_basis_has_fixed_hamming_weight(payload); }

    proc numberSites() : int { return payload.deref().number_sites; }
    proc numberParticles() : int { return payload.deref().number_particles; }
    proc numberUp() : int { return payload.deref().number_up; }

    proc minStateEstimate() : uint(64) { return ls_hs_min_state_estimate(payload); }
    proc maxStateEstimate() : uint(64) { return ls_hs_max_state_estimate(payload); }

    proc representatives() {
      ref rs = payload.deref().representatives;
      if rs.elts == nil then
        halt("basis is not built");
      return makeArrayFromExternArray(rs, uint(64));
    }
  }

  operator Basis.= (ref lhs : Basis, const ref rhs : Basis) {
    assert(lhs.locale == rhs.locale);
    lhs._destroy();
    lhs.payload = ls_hs_clone_basis(rhs.payload);
    lhs.owning = true;
  }

  proc loadBasisFromYaml(filename : string) {
    var ptr = ls_hs_create_spin_basis_from_yaml(filename.localize().c_str());
    if ptr == nil then
      halt("failed to load Basis from " + filename);
    return new Basis(ptr);
  }

  proc SpinBasis(json : string) {
    return new Basis(ls_hs_create_spin_basis_from_json(json.localize().c_str()));
  }
  proc SpinBasis(numberSites : int, hammingWeight : int = -1) {
    var json  = "{ \"number_spins\": " + numberSites:string;
    if hammingWeight != -1 then
      json += ", \"hamming_weight\": " + hammingWeight:string;
    json += " }";
    return SpinBasis(json);
    // return new Basis(ls_hs_create_basis(LS_HS_SPIN, numberSites:c_int, 
    //                                     numberSites:c_int, hammingWeight:c_int));
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

  proc isRepresentative(const ref basis : Basis, const ref alphas : [?D] uint(64),
                        ref are_representatives : [?D2] uint(8),
                        ref norms : [?D3] real(64))
      where D.rank == 1 && D2.rank == 1 && D3.rank == 1 {
    const batchSize = alphas.size;
    assert(are_representatives.size == batchSize);
    assert(norms.size == batchSize);

    ls_hs_is_representative(
      basis.payload,
      batchSize,
      c_const_ptrTo(alphas), 1,
      c_ptrTo(are_representatives),
      c_ptrTo(norms)
    );
  }
  proc isRepresentative(const ref basis : Basis, const ref alphas : [?D] uint(64))
      where D.rank == 1 {
    const batchSize = alphas.size;
    var areRepresentatives : [0 ..# batchSize] uint(8);
    var norms : [0 ..# batchSize] real(64);
    isRepresentative(basis, alphas, areRepresentatives, norms);
    return (areRepresentatives, norms);
  }

  class Operator {
    var payload : c_ptr(ls_hs_operator);
    var basis : Basis;
    var owning : bool;

    proc init(const ref basis : Basis, expression : string, const ref indices : [?D] ?i)
      where D.rank == 2 {

      const c_indices = indices:c_int;
      const c_expr = expression.localize().c_str();
      writeln("Creating operator from '", c_expr:string,  "'");
      print_external_string(c_expr);
      this.payload = ls_hs_create_operator(basis.payload, c_expr,
        indices.dim(0).size:c_int, indices.dim(1).size:c_int, c_const_ptrTo(c_indices));
      writeln("Done creating operator!");
      this.basis = new Basis(this.payload.deref().basis, owning=false);
    }
    proc init(raw : c_ptr(ls_hs_operator), owning : bool = true) {
      this.payload = raw;
      this.basis = new Basis(this.payload.deref().basis, owning=false);
    }
    proc deinit() {
      if owning then
        ls_hs_destroy_operator_v2(payload);
    }

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
      assert(alphas.size <= coeffs.size);
      const batchSize = alphas.size;
      ls_hs_operator_apply_diag_kernel(payload, batchSize, c_const_ptrTo(alphas), 1, c_ptrTo(coeffs));
    }

    proc _localApplyOffDiag(const ref alphas : [?D] uint(64), ref betas : [?D2] uint(64),
                            ref coeffs : [?D3] complex(128))
        where D.rank == 1 && D2.rank == 1 && D3.rank == 1 {
      const batchSize = alphas.size;
      const numberTerms = numberOffDiagTerms();
      assert(betas.size >= batchSize * numberTerms);
      assert(coeffs.size >= batchSize * numberTerms);
      ls_hs_operator_apply_off_diag_kernel(payload, batchSize,
          c_const_ptrTo(alphas), 1, c_ptrTo(betas), 1, c_ptrTo(coeffs));
    }

  }

  proc loadHamiltonianFromYaml(filename : string) {
    var ptr = ls_hs_load_hamiltonian_from_yaml(filename.localize().c_str());
    if ptr == nil then
      halt("failed to load Hamiltonian from " + filename);
    return new Operator(ptr);
  }

  operator +(const ref a : Operator, const ref b : Operator) {
    return new Operator(ls_hs_operator_plus(a.payload, b.payload));
  }



}
