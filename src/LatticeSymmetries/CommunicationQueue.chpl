module CommunicationQueue {



record ToDoBuffer {
  var _dom : domain(1);
  var sigmas : [_dom] uint(64);
  var cs : [_dom] complex(128);
}


class Magic {  
  var _bufferSize : int;
  var _dom : domain(1);
  var _sizes : [LocaleSpace] int;
  var _basisStates : [LocaleSpace] [_dom] uint(64);
  var _coeffs : [LocaleSpace] [_dom] complex(128);
  var _locks : [LocaleSpace] sync bool;

  inline proc _lock(localeIdx : int) { _locks[localeIdx].writeEF(true); }
  inline proc _unlock(localeIdx : int) { _locks[localeIdx].readFE(); }

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
