
use basis;
use wrapper;
use Search;
use CPtr;
use Time;
use states only hash64_01;

config const kBatchSize = 1;

class DistributedOperator {
  var _localBases : [OnePerLocale] ls_hs_spin_basis_v1;
  var _localOperators : [OnePerLocale] ls_hs_operator_v1;

  proc init(path: string) {
    this._localBases = [loc in Locales] new ls_hs_spin_basis_v1(nil, nil);
    this._localOperators = [loc in Locales] new ls_hs_operator_v1(nil, nil);
    complete();

    // for o in this._localOperators {
    //   writeln(o);
    // }
    forall loc in Locales do on loc {
      ls_hs_basis_and_hamiltonian_from_yaml(
        path.localize().c_str(),
        c_ptrTo(this._localBases[loc.id]),
        c_ptrTo(this._localOperators[loc.id]));
    }
    // for o in this._localOperators {
    //   writeln(o);
    // }
  }

  proc deinit() {
    forall loc in Locales do on loc {
      // writeln("ls_hs_destroy_operator on ", this._localOperators[loc.id]);
      ls_hs_destroy_operator(c_ptrTo(this._localOperators[loc.id]));
      // writeln("ls_hs_destroy_basis on ", this._localBases[loc.id]);
      ls_hs_destroy_spin_basis(c_ptrTo(this._localBases[loc.id]));
      // writeln("done: operator=", this._localOperators[loc.id],
      //         " basis=", this._localBases[loc.id]);
    }
    // writeln("[Chapel] DistributedOperator.deinit is done");
  }

  proc rawPtr() { return this._localOperators[here.id].payload; }
  proc maxEntriesPerRow() { return ls_operator_max_buffer_size(rawPtr()):int; }
  proc requiredBufferSize(count : int = 1) {
    return count * maxEntriesPerRow();
  }

  proc batchedApply(
    spins : [?D] uint(64),
    outOffsets : [] uint(64),
    outSpins : [] uint(64),
    outCoeffs : [] complex(128))
     where spins.domain.rank == 1 &&
           outOffsets.domain.rank == 1 &&
           outSpins.domain.rank == 1 &&
           outCoeffs.domain.rank == 1 {
    const count = spins.size;
    // writeln("spins=", spins);
    // writeln("outOffsets=", outOffsets);
    // writeln("outSpins=", outSpins);
    // writeln("outCoeffs=", outCoeffs);
    // const expectedBufferSize = requiredBufferSize(count);
    // assert(outOffsets.size >= count + 1);
    // assert(outSpins.size >= expectedBufferSize);
    // assert(outCoeffs.size >= expectedBufferSize);
    // assert(spins.locale == here);
    // assert(outOffsets.locale == here);
    // assert(outSpins.locale == here);
    // assert(outCoeffs.locale == here);
    const status = ls_hs_operator_apply(
      c_ptrTo(this._localOperators[here.id]), count:uint(64), c_ptrTo(spins),
      c_ptrTo(outOffsets), c_ptrTo(outSpins), c_ptrTo(outCoeffs));
    assert(status == 0);
    // writeln("After:");
    // writeln("spins=", spins);
    // writeln("outOffsets=", outOffsets);
    // writeln("outSpins=", outSpins);
    // writeln("outCoeffs=", outCoeffs);
  }
}

record TaskState {
  const batchSize : int;
  const bufferSize : int;
  var offsetsBatch : [0 ..# batchSize + 1] uint(64);
  var spinsBatch : [0 ..# bufferSize] uint(64);
  var coeffsBatch : [0 ..# bufferSize] complex(128);
}

proc processBatch(count : int,
                  ref info : TaskState,
                  const ref basisStates : BasisStates,
                  const ref X,
                  ref _timer : Timer) {
  type eltType = X[0].eltType;
  var Y : [0 ..# count] complex(128);

  for k in 0 ..# count {
    var yk : complex(128) = 0;
    for _j in info.offsetsBatch[k] .. info.offsetsBatch[k + 1] - 1 {
      const sj = info.spinsBatch[_j:int];

      var xj : eltType;
      const localeId = (hash64_01(sj) % numLocales:uint):int;
      // _timer.start();
      // on Locales[localeId] {
        // const ref representatives =
        //   basisStates.representatives[localeId][0 ..# basisStates.counts[localeId]];
        const (found, j) = binarySearch(
          basisStates.representatives[localeId][0 ..# basisStates.counts[localeId]], sj);
        // writeln("found=", found, ", sj=", sj, ", representatives=", representatives);
        assert(found);
        xj = X[localeId][0, j];
      // }
      // _timer.stop();

      yk += info.coeffsBatch[_j:int] * xj;
    }

    Y[k] = yk;
  }
  return Y;
}

proc batchedApply(matrix : DistributedOperator, ref spins, ref s : TaskState) {
  matrix.batchedApply(
    spins,
    s.offsetsBatch,
    s.spinsBatch,
    s.coeffsBatch);
  // writeln("batchedApply(", spins, "):");
  // writeln("  ", s.offsetsBatch);
  // writeln("  ", s.spinsBatch);
  // writeln("  ", s.coeffsBatch);
}
proc batchedApply(matrix : DistributedOperator, ref spins, ref s : TaskState,
                  ref timer : Timer) {
  timer.start();
  batchedApply(matrix, spins, s);
  timer.stop();
}

proc matvecMinimal(matrix : DistributedOperator, 
                   const ref basisStates : BasisStates)
                  // const ref X, // [OnePerLocale] [?D] ?eltType
                  // ref Y) // [OnePerLocale] [?D] ?eltType
{
  // const batchSize = 1; // max(1, basisStates.size / (10 * loc.maxTaskPar));
  // const bufferSize = 13; // matrix.requiredBufferSize(batchSize);
  // type eltType = X[0].eltType;
  coforall loc in Locales do on loc {
    // ref representatives =
    //   basisStates.representatives[loc.id][0 ..# basisStates.counts[loc.id]];
    // const numberChunks = 1; // (representatives.size + batchSize - 1) / batchSize;
    // var offsetsBatch : [0 .. batchSize] uint(64);
    // var spinsBatch : [0 .. bufferSize - 1] uint(64);
    // var coeffsBatch : [0 .. bufferSize - 1] complex(128);
    // for chunkId in 0 ..# numberChunks {
       //  with (in offsetsBatch, in spinsBatch, in coeffsBatch) {
      // const i = chunkId * batchSize;
      // const n = min(batchSize, basisStates.counts[loc.id] - i);
      // assert(n > 0);
      // assert(n == 1);
      // _batchedApplyTimers[loc.id].start();
      // var offsetsBatch : [0 ..# batchSize + 1] uint(64);
      // var spinsBatch : [0 ..# bufferSize] uint(64);
      // var coeffsBatch : [0 ..# bufferSize] complex(128);
      writeln(basisStates.representatives[loc.id][0 ..# 1]);
      writeln(c_ptrTo(basisStates.representatives[loc.id][0 ..# 1]));
      // matrix.batchedApply(
      //   basisStates.representatives[loc.id][0 ..# 1],
      //   offsetsBatch,
      //   spinsBatch,
      //   coeffsBatch);
    // }
  }
}

proc matvecSerial(matrix : DistributedOperator, 
                  const ref basisStates : BasisStates,
                  const ref X, // [OnePerLocale] [?D] ?eltType
                  ref Y) // [OnePerLocale] [?D] ?eltType
{
  var _totalTimer : Timer = new Timer();
  var _batchedApplyTimers : [OnePerLocale] Timer = new Timer();
  var _searchTimers : [OnePerLocale] Timer = new Timer();
  _totalTimer.start();
  assert(basisStates.representatives.size == numLocales);
  assert(basisStates.counts.size == numLocales);
  assert(X.size == numLocales);
  assert(Y.size == numLocales);

  const batchSize = kBatchSize; // max(1, basisStates.size / (10 * loc.maxTaskPar));
  const bufferSize = matrix.requiredBufferSize(batchSize);
  type eltType = X[0].eltType;
  coforall loc in Locales do on loc {
    ref representatives =
      basisStates.representatives[loc.id][0 ..# basisStates.counts[loc.id]];
    const numberChunks = (representatives.size + batchSize - 1) / batchSize;
    var info = new TaskState(batchSize, bufferSize);
    forall chunkId in 0 ..# numberChunks
        with (in info) {
      const i = chunkId * batchSize;
      const n = min(batchSize, basisStates.counts[loc.id] - i);
      assert(n > 0);
      assert(n == 1);
      batchedApply(matrix, basisStates.representatives[loc.id][i ..# n], info);

      const yk = processBatch(n, info, basisStates, X, _searchTimers[loc.id]);
      if isSubtype(Y[loc.id].eltType, complex) {
        Y[loc.id][0, i ..# n] = yk;
      }
      else {
        assert(&& reduce (yk.im == 0));
        Y[loc.id][0, i ..# n] = yk.re;
      }

      // for k in 0 ..# n {
      //   // writeln("si=", basisStates.representatives[loc.id][i + k]);
      //   // writeln(offsetsBatch);
      //   var yk : complex(128) = 0;
      //   for _j in offsetsBatch[k] .. offsetsBatch[k + 1] - 1 {
      //     const sj = spinsBatch[_j:int];
      //     var xj : eltType;
      //     const localeId = (hash64_01(sj) % numLocales:uint):int;
      //     _searchTimers[loc.id].start();
      //     on Locales[localeId] {
      //       const ref representatives =
      //         basisStates.representatives[localeId][0 ..# basisStates.counts[localeId]];
      //       const (found, j) = binarySearch(representatives, sj);
      //       // writeln("found=", found, ", sj=", sj, ", representatives=", representatives);
      //       assert(found);
      //       xj = X[localeId][0, j];
      //     }
      //     _searchTimers[loc.id].stop();
      //     yk += coeffsBatch[_j:int] * xj;
      //   }

      //   if isSubtype(Y[loc.id].eltType, complex) {
      //     Y[loc.id][0, i + k] = yk;
      //   }
      //   else {
      //     assert(yk.im == 0);
      //     Y[loc.id][0, i + k] = yk.re;
      //   }
      // }
    }
  }

  _totalTimer.stop();
  writeln("[Chapel] matvecSerial took ", _totalTimer.elapsed());
  writeln("[Chapel]     ", _batchedApplyTimers.elapsed(), " spent in matrix.batchedApply");
  writeln("[Chapel]     ", _searchTimers.elapsed(), " spent in searching for indices");
}













