use BitOps;
use List;
use CPtr;
use SysCTypes;
use BlockDist;
use CyclicDist;

use states;

extern proc ls_destroy_flat_spin_basis(ptr: c_ptr(ls_flat_spin_basis));
extern proc ls_get_buffer_size_for_flat_spin_basis(basis: c_ptr(ls_flat_spin_basis)): uint(64);
extern proc ls_serialize_flat_spin_basis(basis: c_ptr(ls_flat_spin_basis),
                                         buffer: [] c_char, size: uint(64)): ls_error_code;
extern proc ls_deserialize_flat_spin_basis(ptr: c_ptr(c_ptr(ls_flat_spin_basis)),
                                           buffer: [] c_char, size: uint(64)): ls_error_code ;

export record WidePtr {
  var ptr: c_void_ptr;
  var loc: chpl_localeID_t;
}

class DistributedBasis {
  const OnePerLocale = LocaleSpace dmapped Block(LocaleSpace);
  var localBases: [OnePerLocale] c_ptr(ls_flat_spin_basis);
  var localStates: [OnePerLocale] PerLocaleState;

  proc init(basis: c_ptr(ls_flat_spin_basis)) {
    this.localBases = [loc in Locales] c_nil;
    this.localStates = [loc in Locales] new PerLocaleState(0);
    complete();

    var size = ls_get_buffer_size_for_flat_spin_basis(basis):int;
    var buffer: [0..<size] c_char;
    ls_serialize_flat_spin_basis(basis, buffer, size:uint);
    for loc in Locales do on loc {
      var localBuffer = buffer;
      ls_deserialize_flat_spin_basis(c_ptrTo(this.localBases[loc.id]), localBuffer, localBuffer.size:uint);
    }
    this.localStates = makeStates(this.localBases);
  }

  proc deinit() {
    for loc in Locales do on loc {
      ls_destroy_flat_spin_basis(this.localBases[loc.id]);
    }
  }
};

export proc ls_distributed_broadcast_basis(widePtr: c_ptr(WidePtr), basis: c_ptr(ls_flat_spin_basis)) {
  var distributedBasis = new unmanaged DistributedBasis(basis);
  ref r = widePtr.deref();
  r.ptr = __primitive("_wide_get_addr", distributedBasis);
  r.loc = __primitive("_wide_get_locale", distributedBasis);
}

proc toDistributedBasis(widePtr: WidePtr) {
  var loc = widePtr.loc;
  var ptr = widePtr.ptr;
  return __primitive("_wide_make", unmanaged DistributedBasis, loc, ptr);
}

export proc ls_distributed_destroy_basis(widePtr: c_ptr(WidePtr)) {
  var distributedBasis = toDistributedBasis(widePtr.deref());
  delete distributedBasis;
}
