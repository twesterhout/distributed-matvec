module CommunicationQueue {

use FFI;
use Types;
use ConcurrentAccessor;

class Magic {
  type eltType;
  var _bufferSize : int;
  var _dom : domain(1);
  var _sizes : [LocaleSpace] int;
  var _basisStates : [LocaleSpace] [_dom] uint(64);
  var _coeffs : [LocaleSpace] [_dom] complex(128);
  var _locks : [LocaleSpace] sync bool;
  var _basis : Basis;
  var _accessor : unmanaged ConcurrentAccessor(eltType);

  proc init(basis : Basis,
            accessor : unmanaged ConcurrentAccessor(?t),
            bufferSize : int) {
    this.eltType = t;
    this._bufferSize = bufferSize;
    this._dom = {0 ..# _bufferSize};
    this._sizes = 0;
    this._basis = basis;
    this._accessor = accessor;
  }

  inline proc _lock(localeIdx : int) { _locks[localeIdx].writeEF(true); }
  inline proc _unlock(localeIdx : int) { _locks[localeIdx].readFE(); }

  proc localProcess(count : int,
                    sigmas : c_ptr(uint(64)),
                    coeffs : c_ptr(complex(128))) {
    if _basis.isStateIndexIdentity() {
      // Simple case when sigmas actually correspond to indices
      foreach k in 0 ..# count {
        const i = sigmas[k]:int;
        const c = coeffs[k]:eltType;
        _accessor.localAdd(i, c);
      }
    }
    var indices : [0 ..# count] int = noinit;
    ls_hs_state_index(
      _basis.payload, count,
      sigmas, 1,
      c_ptrTo(indices), 1);
    foreach k in 0 ..# count {
      const i = indices[k];
      const c = coeffs[k]:eltType;
      if i >= 0 then _accessor.localAdd(i, c:eltType);
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
    copyComplete$.readFF();
    size = 0;
  }

  proc _enqueueUnsafe(localeIdx : int, count : int,
                      sigmas : c_ptr(uint(64)), coeffs : c_ptr(complex(128))) {
    ref offset = _sizes[localeIdx];
    assert(offset + count <= _dom.size);
    foreach i in 0 ..# count {
      _basisStates[localeIdx][offset + i] = sigmas[i];
      _coeffs[localeIdx][offset + i] = coeffs[i];
    }
    offset += count;
    if offset == _dom.size then
      processOnRemote(localeIdx);
  }

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
}

}
