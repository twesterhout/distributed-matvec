module ConcurrentAccessor {

use FFI;
use CTypes;
use ChapelLocks;

config const concurrentAccessorNumLocks = 5 * here.maxTaskPar;

// A simple wrapper around an array `y` which synchronizes accesses. The only
// operation which is allowed is `y[i] += x` <-> `accessor.localAdd(i, x)`. We use
// multiple locks such that multiple threads have lower chance of trying to
// access the same block at the same time. Doing so is still thread-safe, but
// could negatively impact perforamnce especially on machines with high core
// count.
record ConcurrentAccessor {
  type eltType;
  var _data : c_ptr(eltType);
  var _numElts : int;
  var _blockSize : int;
  var numLocks : int;
  var _locks : [0 ..# numLocks] chpl_LocalSpinlock;

  proc init(ref arr : [] ?t, in numLocks : int = concurrentAccessorNumLocks)
      where !arr.domain.stridable && arr.domain.rank == 1 {
    this.eltType = t;
    this._data = c_ptrTo(arr[arr.domain.low]);
    this._numElts = arr.size;
    // It makes no sense to have more locks that there are array elements:
    if numLocks > _numElts then
      numLocks = _numElts;
    this._blockSize = (arr.size + numLocks - 1) / numLocks;
    this.numLocks = numLocks;
    assert(_blockSize * numLocks >= _numElts);
  }

  inline proc _blockIdx(i : int) : int {
    return i / _blockSize;
  }

  inline proc localAdd(i : int, x : eltType) {
    const blockIdx = _blockIdx(i);
    ref lock = _locks[blockIdx];
    lock.lock();
    _data[i] += x;
    lock.unlock();
  }
}

}
