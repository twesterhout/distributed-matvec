module StatesEnumeration {

use Types;
use FFI;

use BitOps;
use List;
// use CTypes;
use CPtr;
use SysCTypes;
use IO;
use DynamicIters;

/* Get next integer with the same hamming weight.

   Semantically equivalent to
   ```
   const m = popcount(v);
   v += 1;
   while (popcount(v) != m) { v += 1; }
   return v;
   ```
 */
private inline proc nextStateFixedHamming(v: uint(64)): uint(64) {
  const t = v | (v - 1);
  return (t + 1) | (((~t & (t + 1)) - 1) >> (ctz(v) + 1));
}

private inline proc nextStateGeneral(v: uint(64)): uint(64) {
  return v + 1;
}

private inline proc nextState(v: uint(64), param isHammingWeightFixed : bool): uint(64) {
  return if isHammingWeightFixed then nextStateFixedHamming(v) else nextStateGeneral(v);
}

/* Fill buffer using either nextStateFixedHamming or nextStateGeneral. Returns
   the number of elements written.
 */
private inline proc manyNextStateImpl(in v: uint(64), bound: uint(64), buffer: [] uint(64),
                                      param isHammingWeightFixed : bool) : int {
  assert(v <= bound, "v is greater than bound");
  for i in buffer.domain {
    buffer[i] = v;
    // If v == bound, incrementing v further is unsafe because it can overflow.
    if (v == bound) { return i + 1; }
    v = nextState(v, isHammingWeightFixed);
  }
  return buffer.size;
}
private proc manyNextState(v: uint(64), bound: uint(64), buffer: [] uint(64),
                           isHammingWeightFixed : bool) : int {
  if isHammingWeightFixed
    then return manyNextStateImpl(v, bound, buffer, true);
    else return manyNextStateImpl(v, bound, buffer, false);
}

private inline proc testBit(bits: uint(64), i: int): bool { return ((bits >> i) & 1):bool; }
private inline proc clearBit(bits: uint(64), i: int) { return bits & ~(1:uint(64) << i); }
private inline proc setBit(bits: uint(64), i: int) { return bits | (1:uint(64) << i); }

/* Get the closest to `x` integer with Hamming weight `hammingWeight`.
 */
private inline proc closestWithFixedHamming(in x: uint(64), hammingWeight: uint): uint(64) {
  assert(hammingWeight <= 64, "Hamming weight too big");
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

config const isRepresentativeBatchSize : int = 65536;

private inline iter findStatesInRange(in lower: uint(64), upper: uint(64),
                                      isHammingWeightFixed : bool,
                                      const ref basis : Basis) {
  if lower > upper then
    return;
  var buffer: [0 ..# isRepresentativeBatchSize] uint(64) = noinit;
  var flags: [0 ..# isRepresentativeBatchSize] uint(8) = noinit;
  var norms: [0 ..# isRepresentativeBatchSize] real(64) = noinit;
  while (true) {
    const written = manyNextState(lower, upper, buffer, isHammingWeightFixed);
    isRepresentative(basis, buffer[0 ..# written], flags[0 ..# written], norms[0 ..# written]);
    for i in 0 ..# written {
      if (flags[i]:bool && norms[i] > 0) { yield buffer[i]; }
    }
    if (buffer[written - 1] == upper) { break; }

    if isHammingWeightFixed { lower = nextState(buffer[written - 1], true); }
    else { lower = nextState(buffer[written - 1], false); }
  }
}

config const enumerateStatesWithProjectionNumChunks : int = 10 * here.maxTaskPar;

/* We want to split `r` into `numChunks` non-overlapping ranges but such that for each range
   `chunk` we have that `popcount(chunk.low) == popcount(chunk.high) == popcount(r.low)` if
   `isHammingWeightFixed` was set to `true`.
 */
iter _customChunks(r : range(uint(64)), numChunks : int, isHammingWeightFixed : bool) {
  logDebug("Creating " + numChunks:string + " chunks for a range with " + r.size:string + " values ...");
  assert(r.size >= numChunks, "more chunks than elements");
  // Start by splitting r into numChunks chunks without worrying about the Hamming weight
  const chunkSize : real = r.size:real / numChunks;
  var offsets : [0 ..# numChunks + 1] uint(64) = noinit;
  foreach i in 0 ..# numChunks {
    offsets[i] = r.low + (i * chunkSize):uint;
    assert(offsets[i] <= r.high, "offsets[i] <= r.high");
  }
  offsets[numChunks] = r.high;

  if isHammingWeightFixed {
    assert(popcount(r.low) == popcount(r.high), "Hamming weight error");
    const hamming = popcount(r.low);
    foreach x in offsets[1 ..# numChunks - 1] do
      x = closestWithFixedHamming(x, hamming);
    foreach x in offsets do
      assert(popcount(x) == popcount(r.low), "Hamming weight error");
  }

  var b = offsets[0];
  var e = offsets[1];
  assert(b <= e);
  yield b .. e;
  for i in 1 ..# numChunks - 1 {
    b = offsets[i];
    e = offsets[i + 1];
    assert(b <= e);
    b = if isHammingWeightFixed then nextState(b, true)
                                else nextState(b, false);
    if b > e then yield e + 1 .. e;
             else yield b .. e;
  }
}


private proc enumerateStatesWithProjection(lower : uint(64), upper : uint(64),
                                           const ref basis : Basis) {
  logDebug("Calling enumerateStatesWithProjection lower=" + lower:string +
           ", upper=" + upper:string + " ...");
  const isHammingWeightFixed = basis.isHammingWeightFixed();
  const numChunks : int = min(enumerateStatesWithProjectionNumChunks, (upper - lower + 1):int);
  var todo : [0 ..# numChunks] range(uint(64)) = _customChunks(lower .. upper, numChunks, isHammingWeightFixed);
  var buffers : [0 ..# numChunks] list(uint(64));
  var countAccumulator : atomic int = 0;
  forall i in dynamic(0 ..# numChunks) {
    const ref chunk = todo[i];
    ref buffer = buffers[i];
    for x in findStatesInRange(chunk.low:uint, chunk.high:uint, isHammingWeightFixed, basis) {
      buffer.append(x);
    }
    countAccumulator.add(buffer.size);
    // writeln(chunk, ": ", buffer);
  }
  // var buffer : list(uint(64));
  // for x in findStatesInRange(lower, upper, isHammingWeightFixed, basis) {
  //   buffer.append(x);
  // }
  const count : int = countAccumulator.read();
  var rs : [0 ..# count] uint(64) = noinit;
  var offset = 0;
  for buffer in buffers {
    rs[offset ..# buffer.size] = buffer;
    offset += buffer.size;
  }
  assert(offset == count);
  return rs;
}

private proc enumerateStatesFixedHamming(lower : uint(64), upper : uint(64),
                                         const ref basis : Basis) {
  logDebug("Calling enumerateStatesFixedHamming lower=" + lower:string +
           ", upper=" + upper:string + " ...");
  assert(popcount(lower) == popcount(upper));
  // The following computes indices of lower and upper
  var _alphas : c_array(uint(64), 2);
  _alphas[0] = lower;
  _alphas[1] = upper;
  var _indices : c_array(c_ptrdiff, 2);
  ls_hs_state_index(basis.payload, 2, _alphas, 1, _indices, 1);
  var lowerIdx = _indices[0];
  var upperIdx = _indices[1];
  // Now that we know the lowerIdx and upperIdx, it's trivial to allocate an
  // array of the right size
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
                                   upper : uint(64) = basis.maxStateEstimate()) {
  if (basis.requiresProjection()) {
    return enumerateStatesWithProjection(lower, upper, basis);
  }
  if (basis.isStateIndexIdentity()) {
    var rs : [0 ..# (upper - lower + 1)] uint(64) = lower .. upper;
    return rs;
  }
  if (basis.isHammingWeightFixed()) {
    return enumerateStatesFixedHamming(lower, upper, basis);
  }

  assert(basis.isSpinfulFermionicBasis());
  const numberSites = basis.numberSites();
  const mask = (1 << numberSites) - 1; // isolate the lower numberSites bits
  const numberUp = basis.numberUp();
  const numberDown = basis.numberParticles() - numberUp;
  const basisA = SpinBasis(numberSites, numberUp);
  const basisB = SpinBasis(numberSites, numberDown);
  // NOTE: Chapel doesn't seem to support recursive functions returning arrays :(
  // const rsA = localEnumerateRepresentatives(basisA, lower & mask, upper & mask);
  // const rsB = localEnumerateRepresentatives(basisB,
  //   (lower >> numberSites) & mask, (upper >> numberSites) & mask);
  const rsA = enumerateStatesFixedHamming(lower & mask, upper & mask, basisA);
  const rsB = enumerateStatesFixedHamming((lower >> numberSites) & mask, (upper >> numberSites) & mask, basisB);

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

// proc _generateBucketRanges() {
//   coforall loc in Locales do on loc {
//     ref localRanges = _ranges[loc.id];
//     ref localStates = _representatives[loc.id][0 ..# _counts[loc.id]];
//     const numberBuckets = 1 << _numberBits;
//     var offset = 0;
//     for i in 0 ..# numberBuckets {
//       localRanges[i] = offset;
//       // writeln("ranges[", loc.id, "][", i, "] = ", offset);
//       while offset != localStates.size && bucketIndex(localStates[offset]) == i {
//         offset += 1;
//       }
//     }
//     localRanges[numberBuckets] = offset;
//     assert(offset == localStates.size);
//   }
// }
// proc localOffsetsForStates(numberBits : int, shift : int,
//                            const ref representatives : [?D] uint(64)) {
//   const numberBuckets = 1 << numberBits;
//   var ranges : [0 ..# numberBuckets + 1] int(64);
//   var offset = 0;
//   for i in 0 ..# numberBuckets {
//     ranges[i] = offset;
//     while offset != representatives.size && (representatives[offset] >> shift):int == i {
//       offset += 1;
//     }
//   }
//   ranges[numberBuckets] = offset;
//   assert(offset == representatives.size);
//   return ranges;
// }

export proc ls_chpl_enumerate_representatives(p : c_ptr(ls_hs_basis),
                                              lower : uint(64),
                                              upper : uint(64),
                                              dest : c_ptr(chpl_external_array)) {
  const basis = new Basis(p, owning=false);
  var rs = localEnumerateRepresentatives(basis, lower, upper);
  dest.deref() = convertToExternalArray(rs);
}

proc initExportedKernels() {
  var kernels = new ls_chpl_kernels(
    c_ptrTo(ls_chpl_enumerate_representatives)
  );
  logDebug("Initializing chpl_kernels ...");
  ls_hs_internal_set_chpl_kernels(c_ptrTo(kernels));
}
initExportedKernels();

}
