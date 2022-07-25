module StatesEnumeration {

use FFI;
use ForeignTypes;
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
use AllLocalesBarriers;
use CommDiagnostics;
use Time;

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
  var timer = new Timer();
  timer.start();
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
  timer.stop();
  logDebug("determineEnumerationRanges took ", timer.elapsed());
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
  logDebug("ls_chpl_enumerate_representatives ...");
  const basis = new Basis(p, owning=false);
  // var rs = localEnumerateRepresentatives(basis, lower, upper);
  var rs = enumerateStates(basis);
  // ref v = rs[here];
  var v = rs[here];
  // writeln(v.type:string);
  // writeln(getExternalArrayType(v):string);
  // writeln(v._value.isDefaultRectangular():string);
  dest.deref() = convertToExternalArray(v);
}

proc determineMergeBounds(const ref basisStates : [] uint(64), numChunks : int) {
  assert(basisStates.size >= numChunks);
  assert(numChunks > 1);
  var bounds : [0 ..# numChunks] uint(64);
  for (chunk, bound) in zip(chunks(0 ..# basisStates.size, numChunks), bounds) {
    bound = basisStates[chunk.high];
  }
  bounds[bounds.domain.low] = 0;
  bounds[bounds.domain.high] = ~(0:uint(64));
  return bounds;
}
// inline proc determineMergeBounds(const ref basisStates : Vector(uint(64)), numChunks : int) {
//   return determineMergeBounds(basisStates.toArray(), numChunks);
// }

proc mergeBoundsToIndexRanges(const ref basisStates : [] uint(64),
                              const ref bounds : [] uint(64)) {
  assert(bounds[bounds.domain.low] <= basisStates[basisStates.domain.low]);
  assert(bounds[bounds.domain.high] >= basisStates[basisStates.domain.high]);
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
  assert(ranges[bounds.size - 1].high == basisStates.size - 1);
  return ranges;
}
// inline proc mergeBoundsToIndexRanges(const ref basisStates : Vector(uint(64)),
//                                      const ref bounds : [] uint(64)) {
//   return mergeBoundsToIndexRanges(basisStates.toArray(), bounds);
// }

proc determineMergeRanges(const ref basisStates : BlockVector(uint(64), 1), numChunks : int) {
  // logDebug("determineMergeRanges: numStates=" + (+ reduce basisStates._counts):string);
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

iter parallelMergeRangesIterator(basisStates : BlockVector(uint(64), 1), numChunks : int) {
  halt("not implemented, but required for compilation :(");
}

iter parallelMergeRangesIterator(param tag: iterKind, basisStates : BlockVector(uint(64), 1),
                                 in numChunks : int)
    where tag == iterKind.standalone {
  const maxNumChunks = min reduce basisStates._counts;
  if numChunks > maxNumChunks then
    numChunks = maxNumChunks;

  const ranges = determineMergeRanges(basisStates, numChunks);
  const offsets = mergeRangesToOffsets(ranges);
  const numRanges = offsets.size - 1;
  const numStates = offsets[numRanges];
  // logDebug("numStates=" + numStates:string);

  const space = {0 ..# numStates};
  const dom = space dmapped Block(boundingBox=space);

  coforall rangeIdx in 0 ..# numRanges {
    const loc = dom.dist.dsiIndexToLocale(offsets[rangeIdx]);
    on loc {
      const localRange = offsets[rangeIdx] .. offsets[rangeIdx + 1] - 1;
      const localChunks = [i in LocaleSpace] ranges[i][rangeIdx];
      // logDebug("yielding: (" + localRange:string + ", " + localChunks:string + ")");
      yield (localRange, localChunks);
    }
  }
}

config const statesFromHashedToBlockNumChunks = 7;

proc statesFromHashedToBlock(const ref basisStates : BlockVector(uint(64), 1)) {
  const numStates = + reduce basisStates._counts;
  const space = {0 ..# numStates};
  const dom = space dmapped Block(boundingBox=space);
  var blockBasisStates : [dom] uint(64);

  forall (currentRange, currentChunks) in
      parallelMergeRangesIterator(basisStates, statesFromHashedToBlockNumChunks) {
    const currentCounts = [chunk in currentChunks] chunk.size;
    var localBasisStates = new BlockVector(uint(64), currentCounts);
    // Get data from remote
    forall (i, chunk) in zip(LocaleSpace, currentChunks) {
      localBasisStates._data[i][0 ..# chunk.size] = basisStates._data[i][chunk];
    }

    // Perform local merge
    var mergedStates : [0 ..# currentRange.size] uint(64) = noinit;
    for ((chunkIndex, indexInChunk), offset) in zip(kMergeIndices(localBasisStates), 0..) {
      mergedStates[offset] = localBasisStates._data[chunkIndex][indexInChunk];
    }
    
    // Copy to remote
    blockBasisStates[currentRange] = mergedStates;
  }

  return blockBasisStates;
}

proc makeBlockDomainForVectors(numVectors, numStates) {
  const boundingBox = {0 ..# numVectors, 0 ..# numStates};
  const targetLocales = reshape(Locales, {0 ..# 1, 0 ..# numLocales});
  const dom : domain(2) dmapped Block(boundingBox=boundingBox, targetLocales=targetLocales) =
    boundingBox;
  return dom;
}

proc vectorsFromHashedToBlock(const ref basisStates : BlockVector(uint(64), 1),
                              const ref vectors : BlockVector(?eltType, 2)) {
  const numVectors = vectors._innerDom.dim(0).size;
  const numStates = + reduce basisStates._counts;
  const dom = makeBlockDomainForVectors(numVectors, numStates);
  var blockVectors : [dom] eltType;

  forall (currentRange, currentChunks) in
      parallelMergeRangesIterator(basisStates, statesFromHashedToBlockNumChunks) {
    const currentCounts = [chunk in currentChunks] chunk.size;
    var localBasisStates = new BlockVector(uint(64), currentCounts);
    var localVectors = new BlockVector(eltType, numVectors, currentCounts);
    // Get data from remote
    forall (i, chunk) in zip(LocaleSpace, currentChunks) {
      localBasisStates._data[i][0 ..# chunk.size] = basisStates._data[i][chunk];
      localVectors._data[i][.., 0 ..# chunk.size] = vectors._data[i][.., chunk];
    }

    // Perform local merge
    var mergedVectors : [0 ..# numVectors, 0 ..# currentRange.size] eltType = noinit;
    for ((chunkIndex, indexInChunk), offset) in zip(kMergeIndices(localBasisStates), 0..) {
      foreach k in 0 ..# numVectors {
        mergedVectors[k, offset] = localVectors._data[chunkIndex][k, indexInChunk];
      }
    }
    
    // Copy to remote
    blockVectors[.., currentRange] = mergedVectors;
  }

  return blockVectors;
}

/*
proc statesAndVectorsFromHashedToBlock(const ref basisStates : BlockVector(uint(64), 1),
                                       const ref vectors : BlockVector(?eltType, 2)) {
  const ranges = determineMergeRanges(basisStates, statesFromHashedToBlockNumChunks);
  const offsets = mergeRangesToOffsets(ranges);
  const numRanges = offsets.size - 1;
  const numStates = offsets[numRanges];
  const numVectors = vectors._innerDom.dim(0).size;

  const statesDomain = {0 ..# numStates} dmapped Block(LocaleSpace);
  const vectorsDomain = makeBlockDomainForVectors(numVectors, numStates);
  var blockBasisStates : [statesDomain] uint(64);
  var blockVectors : [vectorsDomain] eltType;

  coforall rangeIdx in 0 ..# numRanges {
    const loc = statesDomain.dist.dsiIndexToLocale(offsets[rangeIdx]);
    on loc {
      const localRanges = [i in LocaleSpace] ranges[i][rangeIdx];
      const localCounts = [r in localRanges] r.size;
      var localBasisStates = new BlockVector(uint(64), localCounts);
      var localVectors = new BlockVector(eltType, numVectors, localCounts);
      forall (i, r) in zip(LocaleSpace, localRanges) {
        // ref v = localBasisStates.getBlock(i);
        localBasisStates._data[i][0 ..# r.size] = basisStates._data[i][r];
        localVectors._data[i][.., 0 ..# r.size] = vectors._data[i][.., r];
      }

      var mergedSize = offsets[rangeIdx + 1] - offsets[rangeIdx];
      var mergedStates : [0 ..# mergedSize] uint(64) = noinit;
      var mergedVectors : [0 ..# numVectors, 0 ..# mergedSize] eltType = noinit;


      for ((chunkIndex, indexInChunk), offset) in zip(kMergeIndices(localBasisStates), 0..) {
        mergedStates[offset] = localBasisStates._data[chunkIndex][indexInChunk];
        foreach k in 0 ..# numVectors {
          mergedVectors[k, offset] = localVectors._data[chunkIndex][k, indexInChunk];
        }
      }

      blockBasisStates[offsets[rangeIdx] .. offsets[rangeIdx + 1] - 1] = mergedStates;
      blockVectors[.., offsets[rangeIdx] .. offsets[rangeIdx + 1] - 1] = mergedVectors;
    }
  }
  return (blockBasisStates, blockVectors);
}
*/

private proc distributionCounts(const ref basisStates : [] uint(64)) {
  var histogram : [LocaleSpace] int;
  forall basisState in basisStates with (+ reduce histogram) {
    const localeIdx = localeIdxOf(basisState);
    histogram[localeIdx] += 1;
  }

  const dom = LocaleSpace dmapped Block(LocaleSpace);
  var histogramDist : [dom] int = histogram;
  return histogramDist;
}


config const statesFromBlockToHashedNumChunks = 7;

proc statesFromBlockToHashed(const ref basisStates : [] uint(64)) {
  const counts = distributionCounts(basisStates);
  var hashedBasisStates = new BlockVector(uint(64), counts);

  coforall loc in Locales with (ref hashedBasisStates) do on loc {
    ref dest = hashedBasisStates[loc];
    var offset = 0;

    for r in chunks(0 ..# basisStates.size, statesFromBlockToHashedNumChunks) {
      // Download from remote
      var localBasisStates = basisStates[r];
      // Determine which states to copy
      var mask : [0 ..# r.size] bool = noinit;
      forall (x, shouldInclude) in zip(localBasisStates, mask) {
        shouldInclude = localeIdxOf(x) == loc.id;
      }
      // Perform the copy
      for (x, shouldInclude) in zip(localBasisStates, mask) {
        if shouldInclude {
          dest[offset] = x;
          offset += 1;
        }
      }
    }
    assert(offset == dest.size);
  }
  return hashedBasisStates;
}

/*
proc countingSort(src : [] ?t, dest : [] t, extractKey) {
  var counts : c_array(int, 256);
  foreach x in src {
    counts[extractKey(x)] += 1;
  }

  var total = 0;
  for i in 0 ..# 256 {
    const oldCount = counts[i];
    counts[i] = total;
    total += oldCount;
  }

  for x in src {
    const key = extractKey(x);
    dest[counts[key]] = x;
    counts[key] += 1;
  }
}


proc sortByLocaleIdx(const ref basisStates : [] uint(64),
                     ref vectors : [] ?eltType)
    where basisStates.domain.rank == 1 &&
          vectors.domain.rank == 2 {

  const numStates = basisStates.size;
  const numVectors = vectors.dim(0).size;
  if vectors.dim(1).size != numStates then
    halt("incompatible dimensions in sortByLocaleIdx: " + numStates:string +
         " vs. " + vectors.dim(1).size:string);

  if numLocales > 256 then
    halt("We're currently using uint(8) to store locale index. " +
         "Running with more than 256 locales is not yet supported");
  var masks : [0 ..# numStates] uint(8) = noinit;
  foreach (mask, basisState) in zip(masks, basisStates) {
    mask = localeIdxOf(basisState):uint(8);
  }

  var counts : [0 ..# numLocales] int;
  foreach mask in masks {
    counts[mask:int] += 1;
  }

  var offsets : [0 ..# numLocales] int;
  var sum : int = 0;
  for k in 0 ..# numLocales {
    offsets[k] = sum;
    sum += counts[k];
  }

  for i in 0 ..# numStates {
    
  }
}
*/






proc vectorsFromBlockToHashed(const ref basisStates : [] uint(64),
                              const ref vectors : [] ?eltType) {
  const counts = distributionCounts(basisStates);
  const numVectors = vectors.dim(0).size;
  var hashedVectors = new BlockVector(eltType, numVectors, counts);

  coforall loc in Locales with (ref hashedVectors) do on loc {
    ref dest = hashedVectors[loc];
    var offset = 0;

    for r in chunks(0 ..# basisStates.size, statesFromBlockToHashedNumChunks) {
      // Download from remote
      var localBasisStates = basisStates[r];
      var localVectors = vectors[.., r];
      // Determine which states to copy
      var mask : [0 ..# r.size] bool = noinit;
      forall (x, shouldInclude) in zip(localBasisStates, mask) {
        shouldInclude = localeIdxOf(x) == loc.id;
      }
      // Perform the copy
      for i in 0 ..# r.size {
        const shouldInclude = mask[i];
        if shouldInclude {
          foreach k in 0 ..# numVectors {
            dest[k, offset] = localVectors[k, r.low + i];
          }
          offset += 1;
        }
      }
    }
    assert(offset == dest.dim(1).size);
  }
  return hashedVectors;
}


config const numChunksPerLocale = 3;
config const kUseLowLevelComm : bool = true;
config const kVerboseComm : bool = false; 

inline proc GET(addr, node, rAddr, size) {
  __primitive("chpl_comm_get", addr, node, rAddr, size);
}

inline proc PUT(addr, node, rAddr, size) {
  __primitive("chpl_comm_put", addr, node, rAddr, size);
}

proc getPerLocaleCountAndOffset(masks : [] int) {
  if kVerboseComm then startVerboseCommHere();
  var perLocaleCount : [0 ..# numLocales, 0 ..# numLocales] int;
  const perLocaleCountPtr = c_ptrTo(perLocaleCount[0, 0]);
  coforall loc in Locales do on loc {
    const ref myMasks = masks[masks.localSubdomain()];
    var myCount : [0 ..# numLocales] int;
    forall key in myMasks with (+ reduce myCount) {
      myCount[key] += 1;
    }

    if kUseLowLevelComm {
      const putOffset = loc.id * numLocales;
      const putSize = numLocales:c_size_t * c_sizeof(int);
      PUT(c_ptrTo(myCount[0]):c_ptr(int), 0, perLocaleCountPtr + putOffset, putSize);
    }
    else {
      perLocaleCount[loc.id, ..] = myCount;
    }
  }

  var perLocaleOffset : [0 ..# numLocales, 0 ..# numLocales] int;
  foreach destLocaleIdx in 0 ..# numLocales {
    var total = 0;
    for srcLocaleIdx in 0 ..# numLocales {
      const count = perLocaleCount[srcLocaleIdx, destLocaleIdx];
      perLocaleOffset[srcLocaleIdx, destLocaleIdx] = total;
      total += count;
    }
  }

  if kVerboseComm then stopVerboseCommHere();
  return (perLocaleCount, perLocaleOffset);
}

proc getPerTaskCountAndOffset(masksSize : int, masksPtr : c_ptr(int),
                              perLocaleOffset : [] int) {
  if kVerboseComm then startVerboseCommHere();
  var perTaskCount : [0 ..# numChunksPerLocale, 0 ..# numLocales] int;
  var perTaskOffset : [0 ..# numChunksPerLocale, 0 ..# numLocales] int;
  const ranges : [0 ..# numChunksPerLocale] range(int) =
    chunks(0 ..# masksSize, numChunksPerLocale);

  forall (r, chunkIdx) in zip(ranges, 0 ..# numChunksPerLocale) {
    foreach i in r do
      perTaskCount[chunkIdx, masksPtr[i]] += 1;
  }

  for destLocaleIdx in 0 ..# numLocales {
    var total = perLocaleOffset[destLocaleIdx];
    for chunkIdx in 0 ..# numChunksPerLocale {
        const count = perTaskCount[chunkIdx, destLocaleIdx];
        perTaskOffset[chunkIdx, destLocaleIdx] = total;
        total += count;
    }
  }
  if kVerboseComm then stopVerboseCommHere();
  return (perTaskCount, perTaskOffset);
}

proc prefixSum(arr : [] ?eltType) {
  var sums : [0 ..# arr.size] eltType;
  var total : eltType = 0;
  for i in sums.domain {
    sums[i] = total;
    total += arr[i];
  }
  return sums;
}

proc permuteSmall(arrSize : int, masks : c_ptr(int), arr : c_ptr(?eltType),
                  counts : [] int, destOffsets : [] int, destPtrs : [] c_ptr(eltType)) {
  if kVerboseComm then startVerboseCommHere();
  // var prepareTimer = new Timer();
  // var computeTimer = new Timer();
  var copyTimer = new Timer();

  // prepareTimer.start();
  var offsets : [0 ..# numLocales] int = prefixSum(counts);
  var src : [0 ..# arrSize] eltType = noinit;
  // prepareTimer.stop();

  // computeTimer.start();
  for i in 0 ..# arrSize {
    const key = masks[i];
    src[offsets[key]] = arr[i];
    offsets[key] += 1;
  }
  // computeTimer.stop();

  copyTimer.start();
  var i = 0;
  for localeIdx in 0 ..# numLocales {
    if counts[localeIdx] > 0 {
      const srcPtr = c_ptrTo(src[i]);
      const destPtr = destPtrs[localeIdx] + destOffsets[localeIdx];
      const destSize = counts[localeIdx]:c_size_t * c_sizeof(eltType);
      PUT(srcPtr, localeIdx, destPtr, destSize);
      i += counts[localeIdx];
    }
  }
  copyTimer.stop();

  if kVerboseComm then stopVerboseCommHere();
  return copyTimer.elapsed();
}

proc _makeDestArr(arr : [] ?eltType, perLocaleCount)
    where arr.domain.rank == 1 {
  const dom = LocaleSpace dmapped Block(LocaleSpace, Locales);
  const destCounts = [i in 0 ..# numLocales] (+ reduce perLocaleCount[.., i]);
  return new BlockVector(eltType, destCounts, dom);
}
proc _makeDestArr(arr : [] ?eltType, perLocaleCount)
    where arr.domain.rank == 2 {
  const dom = LocaleSpace dmapped Block(LocaleSpace, Locales);
  const destCounts = [i in 0 ..# numLocales] (+ reduce perLocaleCount[.., i]);
  return new BlockVector(eltType, arr.shape[0], destCounts, dom);
}

proc arrFromBlockToHashed(masks : [] int, arr : [] ?eltType) {
  var setupTimer = new Timer();
  setupTimer.start();

  // var perLocaleCountAndOffsetTimer = new Timer();
  // perLocaleCountAndOffsetTimer.start();
  const (perLocaleCount, perLocaleOffset) = getPerLocaleCountAndOffset(masks);
  const perLocaleOffsetPtr = c_const_ptrTo(perLocaleOffset[0, 0]);
  // perLocaleCountAndOffsetTimer.stop();
  // logDebug("arrFromBlockToHashed getPerLocaleCountAndOffset took ",
  //          perLocaleCountAndOffsetTimer.elapsed());

  var makeDestArrTimer = new Timer();
  makeDestArrTimer.start();
  var destArr = _makeDestArr(arr, perLocaleCount);
  makeDestArrTimer.stop();
  // logDebug("arrFromBlockToHashed _makeDestArr took ",
  //          makeDestArrTimer.elapsed());

  var destPtrs : [0 ..# numLocales] c_ptr(eltType);
  const destPtrsPtr = c_const_ptrTo(destPtrs);

  if kVerboseComm then startVerboseComm();
  const batchSize = if arr.domain.rank == 1 then 1 else arr.shape[0];
  const batchIncrement = if arr.domain.rank == 1
                           then 0 else destArr._innerDom.shape[1];

  setupTimer.stop();
  logDebug("arrFromBlockToHashed setup stage took ", setupTimer.elapsed(), "\n",
           "                      from which ", makeDestArrTimer.elapsed(),
           " were spent in _makeDestArr");

  coforall loc in Locales do on loc {
    const _myPtr : c_ptr(eltType) = c_ptrTo(destArr._data[loc.id][destArr._innerDom.low]);
    if kUseLowLevelComm then
      PUT(c_const_ptrTo(_myPtr), 0, destPtrsPtr + loc.id, c_sizeof(c_ptr(eltType)));
    else
      destPtrs[loc.id] = _myPtr;

    allLocalesBarrier.barrier();

    var myPerLocaleOffset : [0 ..# numLocales] int = noinit;
    if kUseLowLevelComm then
      GET(c_ptrTo(myPerLocaleOffset[0]), 0, perLocaleOffsetPtr + loc.id * numLocales,
          numLocales:c_size_t * c_sizeof(int));
    else
      myPerLocaleOffset = perLocaleOffset[loc.id, ..];

    var myDestPtrs : [0 ..# numLocales] c_ptr(eltType) = noinit;
    if kUseLowLevelComm then
      GET(c_ptrTo(myDestPtrs[0]), 0, destPtrsPtr,
          numLocales:c_size_t * c_sizeof(c_ptr(eltType)));
    else
      myDestPtrs = destPtrs;

    const mySubdomain = masks.localSubdomain();
    const myMasksPtr = c_const_ptrTo(masks[mySubdomain.low]);
    const myMasksSize = mySubdomain.size;

    const (perTaskCount, perTaskOffset) =
      getPerTaskCountAndOffset(myMasksSize, myMasksPtr, myPerLocaleOffset);
    const ranges : [0 ..# numChunksPerLocale] range(int) =
      chunks(0 ..# myMasksSize, numChunksPerLocale);

    var computeTimer = new Timer();
    computeTimer.start();
    // var permuteSmallPrepareTime : atomic real;
    // var permuteSmallComputeTime : atomic real;
    var permuteSmallCopyTime : atomic real;

    for batchIdx in 0 ..# batchSize {
      const myArrPtr = if arr.domain.rank == 1
                         then c_const_ptrTo(arr[mySubdomain.low])
                         else c_const_ptrTo(arr[batchIdx, mySubdomain.low]);

      forall (r, chunkIdx) in zip(ranges, 0 ..# numChunksPerLocale) {
        const _copy =
          permuteSmall(r.size,
                       myMasksPtr + r.low,
                       myArrPtr + r.low,
                       perTaskCount[chunkIdx, ..],
                       perTaskOffset[chunkIdx, ..],
                       myDestPtrs);
        // permuteSmallPrepareTime.add(_prep);
        // permuteSmallComputeTime.add(_compute);
        permuteSmallCopyTime.add(_copy, memoryOrder.relaxed);
      }
      myDestPtrs += batchIncrement;
    }

    computeTimer.stop();
    logDebug("arrFromBlockToHashed main loop took ", computeTimer.elapsed(), "\n",
             "                      from which ", permuteSmallCopyTime.read(),
             " were spent in remote PUTs (divide by the number of tasks because ",
             "they ran in parallel)");
  }
  if kVerboseComm then stopVerboseComm();
  return destArr;
}


}
