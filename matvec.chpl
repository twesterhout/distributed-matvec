
use basis;
use wrapper;
use Search;
use CPtr;
use Time;
use states only hash64_01;
use Memory.Diagnostics;

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
  var ysBatch : [0 ..# batchSize] complex(128);
}

proc processBatch(count : int,
                  ref info : TaskState,
                  const ref basisStates : BasisStates,
                  const ref X) {
  var timer = new Timer();
  type eltType = X[0].eltType;

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
        // ref r = basisStates.getRepresentatives(Locales[localeId])[0 ..# basisStates.getCounts(Locales[localeId])];
        // const (found, j) = binarySearch(r, sj);
        // writeln("found=", found, ", sj=", sj, ", representatives=", representatives);
        // assert(found);
        timer.start();
        const j = basisStates.getIndex(sj, Locales[localeId]);
        timer.stop();
        xj = X[localeId][0, j];
      // }
      // _timer.stop();

      yk += info.coeffsBatch[_j:int] * xj;
    }

    info.ysBatch[k] = yk;
  }
  return timer.elapsed();
}

inline proc batchedApply(matrix : DistributedOperator, ref spins, ref s : TaskState) {
  matrix.batchedApply(
    spins,
    s.offsetsBatch,
    s.spinsBatch,
    s.coeffsBatch);
}

proc matvecSerial(matrix : DistributedOperator, 
                  const ref basisStates : BasisStates,
                  const ref X, // [OnePerLocale] [?D] ?eltType
                  ref Y) // [OnePerLocale] [?D] ?eltType
{
  var _totalTimer : Timer = new Timer();
  var _applyTime : [OnePerLocale] real;
  var _searchTime : [OnePerLocale] real;
  var _processTime : [OnePerLocale] real;
  _totalTimer.start();
  // assert(basisStates.representatives.size == numLocales);
  // assert(basisStates.counts.size == numLocales);
  // assert(X.size == numLocales);
  // assert(Y.size == numLocales);

  const batchSize = kBatchSize; // max(1, basisStates.size / (10 * loc.maxTaskPar));
  const bufferSize = matrix.requiredBufferSize(batchSize);
  type eltType = X[0].eltType;

  startVerboseMem();
  coforall loc in Locales do on loc {
    ref representatives = basisStates.getRepresentatives()[0 ..# basisStates.getCounts()];
    // basisStates.representatives[loc.id][0 ..# basisStates.counts[loc.id]];
    const numberChunks = (representatives.size + batchSize - 1) / batchSize;
    var info = new TaskState(batchSize, bufferSize);

    var applyTime : atomic real = 0;
    var processTime : atomic real = 0;
    var searchTime : atomic real = 0;
    forall chunkId in 0 ..# numberChunks with (in info) {
      const i = chunkId * batchSize;
      const n = min(batchSize, representatives.size - i);
      // assert(n > 0);
      // assert(n == 1);
      var applyTimer = new Timer();
      applyTimer.start();
      batchedApply(matrix, representatives[i ..# n], info);
      applyTimer.stop();
      applyTime.add(applyTimer.elapsed());

      var processTimer = new Timer();
      processTimer.start();
      const t = processBatch(n, info, basisStates, X);
      processTimer.stop();
      processTime.add(processTimer.elapsed());
      searchTime.add(t);

      if isSubtype(Y[loc.id].eltType, complex) {
        Y[loc.id][0, i ..# n] = info.ysBatch[0 ..# n];
      }
      else {
        // assert(&& reduce (info.ysBatch[0 ..# n].im == 0));
        // The following is allocating memory... :(
        // Y[loc.id][0, i ..# n] = info.ysBatch[0 ..# n].re;
        // Hence, we do the loop manually
        for k in 0 ..# n {
          assert(info.ysBatch[k].im == 0);
          Y[loc.id][0, i + k] = info.ysBatch[k].re;
        }
      }
    }
    _applyTime[loc.id] = applyTime.read();
    _searchTime[loc.id] = searchTime.read();
    _processTime[loc.id] = processTime.read();
  }
  stopVerboseMem();

  _totalTimer.stop();
  writeln("[Chapel] matvecSerial took ", _totalTimer.elapsed());
  writeln("[Chapel]     ", _applyTime, " spent in matrix.batchedApply");
  writeln("[Chapel]     ", _searchTime, " spent in searching for indices");
  writeln("[Chapel]     ", _processTime, " spent in processing batches");
}













