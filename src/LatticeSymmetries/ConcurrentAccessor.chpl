module ConcurrentAccessor {

use FFI;
use CPtr;
use ChapelLocks;

config const concurrentAccessorNumLocks = 5 * here.maxTaskPar;

record ConcurrentAccessor {
  type eltType;
  var _data : c_ptr(eltType);
  var _numElts : int;
  var _blockSize : int;
  var numLocks : int;
  var _locks : [0 ..# numLocks] chpl_LocalSpinlock;

  proc init(ref arr : [] ?t, numLocks : int = concurrentAccessorNumLocks)
      where !arr.domain.stridable && arr.domain.rank == 1 {
    this.eltTYpe = t;
    this._data = c_ptrTo(arr[arr.domain.low]);
    this._numElts = arr.size;
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
