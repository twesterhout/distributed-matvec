module StatesEnumeration {

use FFI;
use Types;
use Vector;

use BitOps;
use BlockDist;
use ChplConfig;
use CTypes;
use CyclicDist;
use DynamicIters;
use IO;
use List;
use RangeChunk;
use Search;

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

private inline proc unprojectedStateToIndex(basisState : uint(64), isHammingWeightFixed : bool) : int {
  return if isHammingWeightFixed then ls_hs_fixed_hamming_state_to_index(basisState)
                                 else basisState:int;
}
private inline proc unprojectedIndexToState(stateIndex : int, hammingWeight : int) : uint(64) {
  return if hammingWeight != -1 then ls_hs_fixed_hamming_index_to_state(stateIndex, hammingWeight:c_int)
                                else stateIndex:uint(64);
}

/* We want to split `r` into `numChunks` non-overlapping ranges but such that for each range
   `chunk` we have that `popcount(chunk.low) == popcount(chunk.high) == popcount(r.low)` if
   `isHammingWeightFixed` was set to `true`.
 */
proc determineEnumerationRanges(r : range(uint(64)), in numChunks : int,
                                isHammingWeightFixed : bool) {
  const hammingWeight = if isHammingWeightFixed then popcount(r.low):int else -1;
  const lowIdx = unprojectedStateToIndex(r.low, isHammingWeightFixed);
  const highIdx = unprojectedStateToIndex(r.high, isHammingWeightFixed);
  const totalSize = highIdx - lowIdx + 1;
  if numChunks > totalSize then
    numChunks = totalSize;
  var ranges = newCyclicArr({0 ..# numChunks}, range(uint(64)));
  for (r, i) in zip(chunks(lowIdx .. highIdx, numChunks), 0..) {
    ranges[i] = unprojectedIndexToState(r.low, hammingWeight)
                  .. unprojectedIndexToState(r.high, hammingWeight);
  }
  return ranges;
}

/* Hash function which we use to map spin configurations to locale indices.

   Typical usage:
   ```chapel
   var localeIndex = (hash64_01(x) % numLocales:uint):int;
   ```
*/
private inline proc hash64_01(in x: uint(64)): uint(64) {
  x = (x ^ (x >> 30)) * (0xbf58476d1ce4e5b9:uint(64));
  x = (x ^ (x >> 27)) * (0x94d049bb133111eb:uint(64));
  x = x ^ (x >> 31);
  return x;
}

inline proc localeIdxOf(basisState : uint(64)) : int 
    where CHPL_COMM == "" {
  return 0;
}
inline proc localeIdxOf(basisState : uint(64)) : int 
    where CHPL_COMM != "" {
  return (hash64_01(basisState) % numLocales:uint):int;
}

private proc _enumerateStatesUnprojected(r : range(uint(64)), const ref basis : Basis, ref outVectors) {
  const isHammingWeightFixed = basis.isHammingWeightFixed();
  if isHammingWeightFixed && popcount(r.low) != popcount(r.high) then
    halt("r.low=" + r.low:string + " and r.high=" + r.high:string
        + " have different Hamming weight: " + popcount(r.low):string
        + " vs. " + popcount(r.high):string);
  const lowIdx = unprojectedStateToIndex(r.low, isHammingWeightFixed);
  const highIdx = unprojectedStateToIndex(r.high, isHammingWeightFixed);
  const totalCount = highIdx - lowIdx + 1;

  const estimatedCountPerLocale = totalCount / numLocales;
  foreach i in LocaleSpace do
    outVectors[i].reserve(estimatedCountPerLocale);

  var v = r.low;
  for i in 0 ..# totalCount {
    const localeIdx = localeIdxOf(v);
    outVectors[localeIdx].pushBack(v);
    v = if isHammingWeightFixed
          then nextStateFixedHamming(v)
          else nextStateGeneral(v);
  }
}
private proc _enumerateStatesSpinfulFermion(r : range(uint(64)), const ref basis : Basis, ref outVectors) {
  assert(basis.isSpinfulFermionicBasis());
  const numberSites = basis.numberSites();
  const mask = (1 << numberSites) - 1; // isolate the lower numberSites bits
  const numberUp = basis.numberUp();
  const numberDown = basis.numberParticles() - numberUp;
  const countA =
    unprojectedStateToIndex(r.high & mask, isHammingWeightFixed=true) -
      unprojectedStateToIndex(r.low & mask, isHammingWeightFixed=true) + 1;
  const countB =
    unprojectedStateToIndex((r.high >> numberSites) & mask, isHammingWeightFixed=true) -
      unprojectedStateToIndex((r.low >> numberSites) & mask, isHammingWeightFixed=true) + 1;
  const estimatedCountPerLocale = countA * countB / numLocales;
  foreach i in LocaleSpace do
    outVectors[i].reserve(estimatedCountPerLocale);

  var vB = (r.low >> numberSites) & mask;
  for _iB in 0 ..# countB {
    var vA = r.low & mask;
    for _iA in 0 ..# countA {
      const v = (vB << numberSites) | vA;
      const localeIdx = localeIdxOf(v);
      outVectors[localeIdx].pushBack(v);

      vA = nextStateFixedHamming(vA);
    }
    vB = nextStateFixedHamming(vB);
  }
}
private proc _enumerateStatesProjected(r : range(uint(64)), const ref basis : Basis, ref outVectors) {
  const isHammingWeightFixed = basis.isHammingWeightFixed();
  for v in findStatesInRange(r.low, r.high, isHammingWeightFixed, basis) {
    const localeIdx = localeIdxOf(v);
    outVectors[localeIdx].pushBack(v);
  }
}
private proc _enumerateStates(r : range(uint(64)), const ref basis : Basis, ref outVectors) {
  // logDebug("_enumerateStates");
  if basis.requiresProjection() {
    _enumerateStatesProjected(r, basis, outVectors);
  }
  else if basis.isStateIndexIdentity() || basis.isHammingWeightFixed() {
    _enumerateStatesUnprojected(r, basis, outVectors);
  }
  else {
    _enumerateStatesSpinfulFermion(r, basis, outVectors);
  }
}

proc _emptyVectors(type t) : [LocaleSpace] shared Vector(t)
  { return [loc in LocaleSpace] new shared Vector(t); }

config const enumerateStatesNumChunks : int = 10 * numLocales * here.maxTaskPar;

proc enumerateStates(const ref basis : Basis, in numChunks : int, globalRange : range(uint(64))) {
  // logDebug("enumerateStates");
  const isHammingWeightFixed = basis.isHammingWeightFixed();
  const ranges = determineEnumerationRanges(globalRange, numChunks, isHammingWeightFixed);
  assert(ranges.size <= numChunks);
  numChunks = ranges.size; // We could have fewer elements in the range than numChunks.
                           // In such case, we limit numChunks
  const D = {0 ..# numChunks} dmapped Cyclic(startIdx=0);
  var buckets: [D] [LocaleSpace] shared Vector(uint(64)) = [d in D] _emptyVectors(uint(64));
  // noinit;
  // new Vector(uint(64));
 
  forall (r, i) in zip(ranges, 0..) {
    assert(here == buckets[i].locale);
    const localBasis = basis;
    _enumerateStates(r, localBasis, buckets[i]);
  }

  const OnePerLocale = LocaleSpace dmapped Block(LocaleSpace);
  var representatives : [OnePerLocale] shared Vector(uint(64)) =
    [loc in OnePerLocale] new shared Vector(uint(64));
  forall (vector, localeIdx) in zip(representatives, LocaleSpace) {
    vector.reserve(+ reduce [bucket in buckets] bucket[localeIdx].size);
    for bucket in buckets {
      const ref v = bucket[localeIdx];
      vector.append(v);
    }
  }
  return new BlockVector(representatives);
}
proc enumerateStates(const ref basis : Basis, numChunks : int = enumerateStatesNumChunks) {
  const lower = basis.minStateEstimate();
  const upper = basis.maxStateEstimate();
  return enumerateStates(basis, numChunks, lower .. upper);
}

export proc ls_chpl_enumerate_representatives(p : c_ptr(ls_hs_basis),
                                              lower : uint(64),
                                              upper : uint(64),
                                              dest : c_ptr(chpl_external_array)) {
  const basis = new Basis(p, owning=false);
  // var rs = localEnumerateRepresentatives(basis, lower, upper);
  var rs = enumerateStates(basis);
  ref v = rs[here];
  // v.shrink();
  dest.deref() = convertToExternalArray(v);
}

proc initExportedKernels() {
  var kernels = new ls_chpl_kernels(
    c_ptrTo(ls_chpl_enumerate_representatives)
  );
  logDebug("Initializing chpl_kernels ...");
  ls_hs_internal_set_chpl_kernels(c_ptrTo(kernels));
}


proc determineMergeBounds(const ref basisStates : [] uint(64), numChunks : int) {
  assert(basisStates.size >= numChunks);
  var bounds : [0 ..# numChunks] uint(64);
  for (chunk, bound) in zip(chunks(0 ..# basisStates.size, numChunks), bounds) {
    bound = basisStates[chunk.high];
  }
  return bounds;
}
// inline proc determineMergeBounds(const ref basisStates : Vector(uint(64)), numChunks : int) {
//   return determineMergeBounds(basisStates.toArray(), numChunks);
// }

proc mergeBoundsToIndexRanges(const ref basisStates : [] uint(64),
                              const ref bounds : [] uint(64)) {
  var ranges : [bounds.domain] range(int);
  for (r, b, i) in zip(ranges, bounds, 0..) {
    var loIdx = if i == 0 then 0 else ranges[i - 1].high + 1;
    var (found, hiIdx) = binarySearch(basisStates, b, lo=loIdx);
    if !found then
      hiIdx -= 1;
    r = loIdx .. hiIdx;
    // if hiIdx - loIdx >= 0 {
    //   if i != 0 {
    //     assert(basisStates[loIdx] > bounds[i - 1]);
    //   }
    //   assert(basisStates[hiIdx] <= b);
    // }
  }
  return ranges;
}
// inline proc mergeBoundsToIndexRanges(const ref basisStates : Vector(uint(64)),
//                                      const ref bounds : [] uint(64)) {
//   return mergeBoundsToIndexRanges(basisStates.toArray(), bounds);
// }

proc determineMergeRanges(const ref basisStates : BlockVector(uint(64), 1), numChunks : int) {
  const bounds = determineMergeBounds(basisStates[here], numChunks);
  var ranges : [basisStates._outerDom] [0 ..# bounds.size] range(int);
  forall (r, i) in zip(ranges, 0 ..) {
    const localBounds = bounds;
    r = mergeBoundsToIndexRanges(basisStates.getBlock(i), localBounds);
  }
  return ranges;
}

proc mergeRangesToOffsets(const ref ranges) {
  const numRanges = ranges[0].size;
  var offsets : [0 ..# numRanges + 1] int;
  offsets[0] = 0;
  for i in 0 ..# numRanges {
    offsets[i + 1] = offsets[i] + (+ reduce [k in ranges.domain] ranges[k][i].size);
  }
  return offsets;
}

// config const statesFromHashedToBlockNumChunks = 7;

// proc statesFromHashedToBlock(const ref basisStates : BlockVector(uint(64), 1)) {
//   const ranges = determineMergeRanges(basisStates);
//   const offsets = mergeRangesToOffsets(ranges);
//   const numRanges = offsets.size - 1;
//   const numStates = offsets[numRanges];
//   const dom = {0 ..# numStates} dmapped Block(LocaleSpace);
// 
//   var blockBasisStates : [dom] uint(64);
//   coforall rangeIdx in 0 ..# numRanges {
//     const loc = dom.dsiIndexToLocale(offsets[rangeIdx]);
//     on loc {
//       const 
//       var localBasisStates : {0 ..# basisStates.size, 
// 
// 
//     }
//   }
// }


}
