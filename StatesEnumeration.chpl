module StatesEnumeration {

use ApplyOperator;
use CPtr;
use SysCTypes;
use BitOps;

/* Get next integer with the same hamming weight.

   Semantically equivalent to
   ```
   const m = popcount(v);
   v += 1;
   while (popcount(v) != m) { v += 1; }
   return v;
   ```
 */
inline proc nextStateFixedHamming(v: uint(64)): uint(64) {
  const t = v | (v - 1);
  return (t + 1) | (((~t & (t + 1)) - 1) >> (ctz(v) + 1));
}

inline proc nextStateGeneral(v: uint(64)): uint(64) {
  return v + 1;
}

inline proc nextState(v: uint(64), isHammingWeightFixed : bool): uint(64) {
  return if isHammingWeightFixed then nextStateFixedHamming(v) else nextStateGeneral(v);
}

/* Fill buffer using either nextStateFixedHamming or nextStateGeneral. Returns
   the number of elements written.
 */
inline proc manyNextState(in v: uint(64), bound: uint(64), buffer: [] uint(64),
                          isHammingWeightFixed : bool) {
  assert(v <= bound);
  for i in buffer.domain {
    buffer[i] = v;
    // If v == bound, incrementing v further is unsafe because it can overflow.
    if (v == bound) { return i + 1; }
    v = nextState(v, isHammingWeightFixed);
  }
  return buffer.size;
}

inline proc testBit(bits: uint(64), i: int): bool { return ((bits >> i) & 1):bool; }
inline proc clearBit(bits: uint(64), i: int) { return bits & ~(1:uint(64) << i); }
inline proc setBit(bits: uint(64), i: int) { return bits | (1:uint(64) << i); }

/* Get the closest to `x` integer with Hamming weight `hammingWeight`.
 */
inline proc closestWithFixedHamming(in x: uint(64), hammingWeight: uint): uint(64) {
  assert(hammingWeight <= 64);
  var weight = popcount(x);
  if (weight > hammingWeight) {
      // Keep clearing lowest bits until we reach the desired Hamming weight.
      var i = 0;
      while (weight > hammingWeight) {
        assert(i < 64);
        if (testBit(x, i)) {
            x = clearBit(x, i);
            weight -= 1;
        }
        i += 1;
      }
      // NOTE: why??
      // var maxValue = hammingWeight == 0 ? (0:uint(64)) : (~(0:uint(64)) << (64 - hammingWeight));
      // if (x < max_value) { x = next_state<true>(x); }
  }
  else if (weight < hammingWeight) {
    // Keep setting lowest bits until we reach the desired Hamming weight.
    var i = 0;
    while (weight < hammingWeight) {
      assert(i < 64);
      if (!testBit(x, i)) {
          x = setBit(x, i);
          weight += 1;
      }
      i += 1;
    }
  }
  return x;
}

proc enumerateStatesFixedHamming(lower : uint(64), upper : uint(64), basis : c_ptr(ls_hs_basis)) {
  assert(popcount(lower) == popcount(upper));
// typedef void (*ls_hs_internal_state_index_kernel_type)(
//     ptrdiff_t batch_size, uint64_t const *alphas, ptrdiff_t alphas_stride,
//     ptrdiff_t *indices, ptrdiff_t indices_stride, void const *private_data);
  var _alphas : c_array(uint(64), 2);
  _alphas[0] = lower;
  _alphas[1] = upper;
  var _indices : c_array(c_ptrdiff, 2);
  ls_hs_state_index(basis, 2, _alphas, 1, _indices, 1);
  var lowerIdx = _indices[0]; // ls_hs_internal_rank_via_combinadics(lower, binomials);
  var upperIdx = _indices[1]; // ls_hs_internal_rank_via_combinadics(upper, binomials);
  var rs : [0 ..# (upperIdx - lowerIdx + 1)] uint(64) = noinit;
  var v = lower;
  for i in rs.domain {
    rs[i] = v;
    v = nextStateFixedHamming(v);
  }
  return rs;
}

proc localEnumerateRepresentatives(const ref basis : Basis,
                                   lower : uint(64) = basis.minStateEstimate(),
                                   upper : uint(64) = basis.maxStateEstimate()) : [] uint(64) {
  writeln("lower = ", lower, ", upper = ", upper);
  if (basis.isStateIndexIdentity()) {
    var rs : [0 ..# (upper - lower + 1)] uint(64) = lower .. upper;
    return rs;
  }
  if (ls_hs_basis_has_fixed_hamming_weight(basis.payload)) {
    return enumerateStatesFixedHamming(lower, upper, basis.payload);
  }

  assert(basis.isSpinfulFermionicBasis());
  const numberSites = basis.numberSites();
  const mask = (1 << numberSites) - 1; // isolate the lower numberSites bits
  const numberUp = basis.numberUp();
  const numberDown = basis.numberParticles() - numberUp;
  writeln("numberUp = ", numberUp, ", numberDown = ", numberDown);

  const basisA = SpinBasis(numberSites, numberUp);
  const basisB = SpinBasis(numberSites, numberDown);
  // NOTE: Chapel doesn't seem to support recursive functions returning arrays :(
  // const rsA = localEnumerateRepresentatives(basisA, lower & mask, upper & mask);
  // const rsB = localEnumerateRepresentatives(basisB,
  //   (lower >> numberSites) & mask, (upper >> numberSites) & mask);
  const rsA = enumerateStatesFixedHamming(lower & mask, upper & mask, basisA.payload);
  writeln("rsA: ", rsA);
  const rsB = enumerateStatesFixedHamming((lower >> numberSites) & mask, (upper >> numberSites) & mask,
                                          basisB.payload);
  writeln("rsB: ", rsB);

  var rs : [0 ..# (rsA.size * rsB.size)] uint(64) = noinit;
  var offset : int = 0;
  for b in rsB {
    for a in rsA {
      rs[offset] = (b << numberSites) | a;
      offset += 1;
    }
  }
  return rs;
}

export proc ls_chpl_enumerate_representatives(p : c_ptr(ls_hs_basis),
                                              lower : uint(64),
                                              upper : uint(64)) : chpl_external_array {
  const basis = new Basis(p, owning=false);
  var rs = localEnumerateRepresentatives(basis, lower, upper);
  return convertToExternalArray(rs);
}

proc initExportedKernels() {
  var kernels = new ls_chpl_kernels(
    c_ptrTo(ls_chpl_enumerate_representatives)
  );
  writeln("Setting chpl_kernels ...");
  ls_hs_internal_set_chpl_kernels(c_ptrTo(kernels));
}
initExportedKernels();

}
