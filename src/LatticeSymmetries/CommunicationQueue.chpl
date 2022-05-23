module CommunicationQueue {

use CTypes;

use FFI;
use Types;
use ConcurrentAccessor;

config const communicationQueueBufferSize = 100;
config const stagingBuffersBufferSize = 10;

// Handles the non-trivial communication of matrix elements between locales
class CommunicationQueue {
  type eltType;
  var _bufferSize : int;
  var _dom : domain(1);
  var _sizes : [LocaleSpace] int;
  var _basisStates : [LocaleSpace] [_dom] uint(64);
  var _coeffs : [LocaleSpace] [_dom] complex(128);
  var _locks : [LocaleSpace] sync bool;
  var _basisPtr : c_ptr(Basis);
  var _accessorPtr : c_ptr(ConcurrentAccessor(eltType));

  proc init(const ref basis : Basis,
            ref accessor : ConcurrentAccessor(?t),
            bufferSize : int = communicationQueueBufferSize) {
    this.eltType = t;
    this._bufferSize = bufferSize;
    this._dom = {0 ..# _bufferSize};
    this._sizes = 0;
    this._basisPtr = c_ptrTo(basis);
    this._accessorPtr = c_ptrTo(accessor);
  }

  inline proc _lock(localeIdx : int) { _locks[localeIdx].writeEF(true); }
  inline proc _unlock(localeIdx : int) { _locks[localeIdx].readFE(); }

  proc localProcess(count : int,
                    sigmas : c_ptr(uint(64)),
                    coeffs : c_ptr(complex(128))) {
    if _basisPtr.deref().isStateIndexIdentity() {
      // Simple case when sigmas actually correspond to indices
      ref accessor = _accessorPtr.deref();
      foreach k in 0 ..# count {
        const i = sigmas[k]:int;
        const c = coeffs[k]:eltType;
        accessor.localAdd(i, c);
      }
    }
    // We do not know the indices of `sigmas` yet:
    // TODO: is the fact that we're allocating `indices` over and over again
    // okay performance-wise?
    var indices : [0 ..# count] int = noinit;
    ls_hs_state_index(
      _basisPtr.deref().payload, count,
      sigmas, 1,
      c_ptrTo(indices), 1);
    ref accessor = _accessorPtr.deref();
    foreach k in 0 ..# count {
      const i = indices[k];
      const c = coeffs[k]:eltType;
      // Importantly, the user could have made a mistake and given us an
      // operator which does not respect the basis symmetries. Then we could
      // have that a |σ⟩ was generated that doesn't belong to our basis. In
      // this case, we should throw an error.
      if i >= 0 then accessor.localAdd(i, c:eltType);
                else halt("invalid index");
    }
  }
  inline proc localProcess(const ref sigmas : [] uint(64),
                           const ref coeffs : [] complex(128)) {
    assert(sigmas.size == coeffs.size);
    localProcess(sigmas.size, c_const_ptrTo(sigmas), c_const_ptrTo(coeffs));
  }

  proc processOnRemote(localeIdx : int) {
    ref size = _sizes[localeIdx];
    const ref sigmas = _basisStates[localeIdx][0 ..# size];
    const ref cs = _coeffs[localeIdx][0 ..# size];
    var copyComplete$ : single bool;
    begin on Locales[localeIdx] {
      const basisStates = sigmas;
      const coeffs = cs;
      copyComplete$.writeEF(true);
      localProcess(basisStates, coeffs);
    }
    // We need to wait for the remote copy to complete before we can reuse
    // `sigmas` and `cs`.
    copyComplete$.readFF();
    size = 0;
  }

  inline proc _enqueueUnsafe(localeIdx : int, count : int,
                             sigmas : c_ptr(uint(64)), coeffs : c_ptr(complex(128))) {
    ref offset = _sizes[localeIdx];
    assert(offset + count <= _dom.size);
    // TODO: Why a loop here? Can't we have a simple memcpy call?
    foreach i in 0 ..# count {
      _basisStates[localeIdx][offset + i] = sigmas[i];
      _coeffs[localeIdx][offset + i] = coeffs[i];
    }
    offset += count;
    // So far, everything was done locally. Only `processOnRemote` involves
    // communication.
    if offset == _dom.size then
      processOnRemote(localeIdx);
  }

  // Enqueue elements which need to be added to the output vector |y⟩. This
  // function effectively performs |y⟩ += cᵢ|σᵢ⟩ for each cᵢ in `coeffs` and σᵢ
  // in `sigmas` (`i` is in `0 ..# count`). `localeIdx` specifies the locale to
  // which `sigmas` and `coeffs` should be sent.
  proc enqueue(localeIdx : int, in count : int,
               sigmas : c_ptr(uint(64)), coeffs : c_ptr(complex(128))) {
    if localeIdx == here.id {
      localProcess(count, sigmas, coeffs);
      return;
    }
    _lock(localeIdx);
    while count > 0 {
      const remaining = min(_dom.size - _sizes[localeIdx], count);
      _enqueueUnsafe(localeIdx, remaining, sigmas, coeffs);
      count -= remaining;
    }
    _unlock(localeIdx);
  }

  proc drain() {
    forall localeIdx in LocaleSpace {
      if _sizes[localeIdx] > 0 {
        _lock(localeIdx);
        processOnRemote(localeIdx);
        _unlock(localeIdx);
      }
    }
    // Check
    foreach localeIdx in LocaleSpace do
      assert(_sizes[localeIdx] == 0);
  }
}

record StagingBuffers {
  type eltType;
  var _capacity : int;
  var _sizes : [0 ..# numLocales] int;
  var _basisStates : [0 ..# numLocales, 0 ..# _capacity] uint(64);
  var _coeffs : [0 ..# numLocales, 0 ..# _capacity] eltType;
  var _magic : c_ptr(CommunicationQueue(eltType));

  proc init(ref magic : CommunicationQueue(?t),
            bufferSize : int = stagingBuffersBufferSize) {
    this.eltType = t;
    this._capacity = bufferSize;
    this._sizes = 0;
    this._magic = c_ptrTo(magic);
  }
  proc init=(const ref other : StagingBuffers) {
    assert(other.locale == here); // we do not support assignment from remote
    this.eltType = other.eltType;
    this._capacity = other._capacity;
    this._magic = other._magic;
    foreach i in LocaleSpace {
      const n = other._sizes[i];
      if n != 0 {
        this._sizes[i] = n;
        // TODO: Are these again crazy slow? Do we care?
        this._basisStates[i, 0 ..# n] = other._basisStates[i, 0 ..# n];
        this._coeffs[i, 0 ..# n] = other._coeffs[i, 0 ..# n];
      }
    }
  }

  // Push all elements remaining in staging buffers for locale 'localIdx' to
  // '_magic'.
  proc flush(localeIdx : int) {
    _magic.enqueue(localeIdx, _sizes[localeIdx],
                   c_const_ptrTo(_basisStates[localeIdx, 0]),
                   c_const_ptrTo(_coeffs[localeIdx, 0]));
    _sizes[localeIdx] = 0;
  }
  // Push all elements remaining in staging buffers to '_magic'
  proc flush() {
    foreach localeIdx in LocaleSpace do
      flush(localeIdx);
  }

  proc add(localeIdx : int, basisState : uint(64), coeff : eltType) {
    ref n = _sizes[localeIdx];
    assert(n < _capacity);
    _basisStates[localeIdx, n] = basisState;
    _coeffs[localeIdx, n] = coeff;
    n += 1;
    if n == _capacity then
      flush(localeIdx);
  }
  proc add(batchSize : int,
           basisStatesPtr : c_ptr(uint(64)),
           coeffsPtr : c_ptr(?t)) {
    // This function could potentially be vectorized, because
    // computing hashes should make it compute bound.
    for i in 0 ..# batchSize {
      const localeIdx = localeOf(basisStatesPtr[i]);
      // TODO: what about conversions from/to complex(128)?
      add(localeIdx, basisStatesPtr[i], coeffsPtr[i]);
    }
  }
}

}
