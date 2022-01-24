
use basis;
use wrapper;
use Search;

class DistributedOperator {
  var _localOperators : [OnePerLocale] ls_hs_operator_v1;

  proc init(path: string) {
    this._localOperators = [loc in Locales] new ls_hs_operator_v1(nil, nil);
    complete();

    coforall loc in Locales do on loc {
      const localPath = path;
      var dummyBasis = new ls_hs_spin_basis_v1(nil, nil);
      ls_hs_basis_and_hamiltonian_from_yaml(
        localPath.c_str(), c_ptrTo(dummyBasis), c_ptrTo(this._localOperators[loc.id]));
      ls_hs_destroy_spin_basis(c_ptrTo(dummyBasis));
    }
  }

  proc deinit() {
    coforall loc in Locales do on loc {
      ls_hs_destroy_operator(c_ptrTo(this._localOperators[loc.id]));
    }
  }

  inline proc rawPtr() { return this._localOperators[here.id].payload; }
  inline proc maxEntriesPerRow() { return ls_operator_max_buffer_size(rawPtr()):int; }
  inline proc requiredBufferSize(count : int = 1) {
    return count * maxEntriesPerRow();
  }

  inline proc batchedApply(
    spins : [?D] uint(64),
    outOffsets : [] uint(64),
    outSpins : [] uint(64),
    outCoeffs : [] complex(128)) {
    assert(spins.domain.rank == 1);
    assert(outOffsets.domain.rank == 1);
    assert(outSpins.domain.rank == 1);
    assert(outCoeffs.domain.rank == 1);
    const count = spins.size;
    const expectedBufferSize = requiredBufferSize(count);
    assert(outOffsets.size >= expectedBufferSize + 1);
    assert(outSpins.size >= expectedBufferSize);
    assert(outCoeffs.size >= expectedBufferSize);
    const status = ls_hs_operator_apply(
      rawPtr(), count:uint(64), c_ptrTo(spins),
      c_ptrTo(offsets), c_ptrTo(outSpins), c_ptrTo(outCoeffs));
    assert(status == 0);
  }
}

proc matvecSerial(const ref operator : DistributedOperator, 
                  const ref basisStates : BasisStates,
                  const ref X, // [OnePerLocale] [?D] ?eltType
                  ref Y) // [OnePerLocale] [?D] ?eltType
{
  const batchSize = 1;
  const bufferSize = operator.requiredBufferSize(batchSize);
  forall loc in Locales do on loc {
    const numberChunks = basisStates.counts[loc.id] / batchSize;
    for chunkId in 0 .. numberChunks - 1 {
      const i = chunkId * batchSize;
      offsetsBatch : [0 .. bufferSize - 1] uint(64);
      spinsBatch : [0 .. bufferSize - 1] uint(64);
      coeffsBatch : [0 .. bufferSize - 1] complex(128);
      operator.batchedApply(
        basisStates.representatives[i ..# batchSize],
        offsetsBatch,
        spinsBatch,
        coeffsBatch);
      for k in 0 .. batchSize - 1 {
        var yk : complex(128) = 0;
        for _j in offsetsBatch[k] .. offsetsBatch[k + 1] {
          const sj = spinsBatch[_j];
          var xj : complex(128);
          const localeId = (hash64_01(sj) % numLocales:uint):int;
          on Locales[localeId] {
            const ref representatives =
              basisStates.representatives[localeId][0 ..# basisStates.counts[localeId]];
            const (found, j) = binarySearch(representatives, sj);
            assert(found);
            xj = X[localeId][j];
          }
          yk += coeffsBuffer[_j] * xj;
        }
        Y[loc.id][i + k] = yk;
      }
    }
    // process the rest
    {
      assert(basisStates.counts[loc.id] % batchSize == 0;
    }
  }
}













