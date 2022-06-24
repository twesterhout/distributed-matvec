module BatchedOperator {

use CTypes;
use FFI;
use ForeignTypes;

// TODO: this is currently implemented inefficiently
private proc localCompressMultiply(batchSize : int, numberTerms : int,
                                   sigmas : c_ptr(uint(64)), cs : c_ptr(complex(128)),
                                   xs : c_ptr(?eltType), offsets : c_ptr(int)) {
  var offset : int = 0;
  var i : int = 0;

  offsets[0] = 0;
  for batchIdx in 0 ..# batchSize {
    for termIdx in 0 ..# numberTerms {
      if (cs[i] != 0) {
        if (i > offset) {
          cs[offset] = cs[i];
          sigmas[offset] = sigmas[i];
        }
        cs[offset] *= xs[batchIdx];
        offset += 1;
      }
      i += 1;
    }
    offsets[batchIdx + 1] = offset;
  }
}

record BatchedOperator {
  var _matrixPtr : c_ptr(Operator);
  var batchSize : int;
  var _numberOffDiagTerms : int;
  var _dom : domain(1);
  var _spins1 : [_dom] uint(64);
  var _spins2 : [_dom] uint(64);
  var _coeffs1 : [_dom] complex(128);
  var _coeffs2 : [_dom] complex(128);
  var _norms : [_dom] real(64);
  var _offsets : [0 ..# batchSize + 1] int;

  proc init(const ref matrix : Operator, batchSize : int) {
    this._matrixPtr = c_const_ptrTo(matrix);
    this.batchSize = batchSize;
    this._numberOffDiagTerms = matrix.numberOffDiagTerms();
    const numberTerms = max(_numberOffDiagTerms, 1);
    this._dom = {0 ..# (batchSize * (numberTerms + 1))};
    // logDebug("BatchedOperator._dom = ", this._dom);
  }
  proc init=(const ref other : BatchedOperator) {
    assert(other.locale == here);
    this._matrixPtr = other._matrixPtr;
    this.batchSize = other.batchSize;
    this._numberOffDiagTerms = other._numberOffDiagTerms;
    this._dom = other._dom;
  }

  proc computeOffDiag(count : int,
                      alphas : c_ptr(uint(64)),
                      xs : c_ptr(?eltType))
      : (c_ptr(uint(64)), c_ptr(complex(128)), c_ptr(int)) {
    assert(count <= batchSize);
    const ref matrix = _matrixPtr.deref();
    // Simple case when no symmetries are used
    if !matrix.basis.requiresProjection() {
      const betas = c_pointer_return(_spins1[0]);
      const cs = c_pointer_return(_coeffs1[0]);
      const offsets = c_pointer_return(_offsets[0]);
      // We compute H|αᵢ⟩ = ∑ⱼ cᵢⱼ|βᵢⱼ⟩ for each |αᵢ⟩ where i ∈ {0 ..# count}
      // at this point j ∈ {0 ..# _numberOffDiagTerms}. Both cᵢⱼ and |βᵢⱼ⟩ are
      // conceptually 2-dimensional arrays, but we use flattened
      // representations of them.
      ls_hs_operator_apply_off_diag_kernel(
        matrix.payload,
        count,
        alphas, 1,
        betas, 1,
        cs);
      // Since many cᵢⱼ could be zero, we compress |βᵢⱼ⟩ and cᵢⱼ arrays by
      // throwing away all zero elements. Furthermore, we multiply each cᵢⱼ by
      // the corresponding xᵢ since we have to do it at some point anyway.
      localCompressMultiply(
        count, _numberOffDiagTerms,
        betas, cs, xs,
        offsets);
      return (betas, cs, offsets);
    }
    // The tricky case when we have to project betas first
    const tempSpins = c_pointer_return(_spins2[0]);
    const tempCoeffs = c_pointer_return(_coeffs2[0]);
    const offsets = c_pointer_return(_offsets[0]);
    ls_hs_operator_apply_off_diag_kernel(
      matrix.payload,
      count,
      alphas, 1,
      tempSpins, 1,
      tempCoeffs);
    localCompressMultiply(
      count, _numberOffDiagTerms,
      tempSpins, tempCoeffs, xs,
      offsets);
    const totalSize = offsets[count];
    // we are also interested in norms of alphas, so we append them to tempSpins
    foreach i in 0 ..# count {
      tempSpins[totalSize + i] = alphas[i];
    }

    const betas = c_pointer_return(_spins1[0]);
    const cs = c_pointer_return(_coeffs1[0]);
    const norms = c_pointer_return(_norms[0]);
    ls_hs_state_info(
      matrix.basis.payload, totalSize + batchSize,
      tempSpins, 1,
      betas, 1,
      cs,
      norms);
    foreach i in 0 ..# count {
      foreach k in offsets[i] ..< offsets[i + 1] {
        cs[k] *= tempCoeffs[k] * norms[k] / norms[totalSize + i];
      }
    }
    return (betas, cs, offsets);
  }

}

} // end module BatchedOperator
