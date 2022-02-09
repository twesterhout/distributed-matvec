use List;
use CPtr;
use SysCTypes;
use BlockDist;
use CyclicDist;
use Time;
use VisualDebug;
use RangeChunk;
use Search;

use wrapper;
use states;
use Merge;

const OnePerLocale = LocaleSpace dmapped Block(LocaleSpace);

class DistributedBasis {
  var _localBases: [OnePerLocale] ls_hs_spin_basis_v1;

  proc init(path: string) {
    this._localBases = [loc in Locales] new ls_hs_spin_basis_v1(nil, nil);
    complete();

    coforall loc in Locales do on loc {
      var dummyHamiltonian = new ls_hs_operator_v1(nil, nil);
      ls_hs_basis_and_hamiltonian_from_yaml(
        path.localize().c_str(), c_ptrTo(this._localBases[loc.id]), c_ptrTo(dummyHamiltonian));
      ls_hs_destroy_operator(c_ptrTo(dummyHamiltonian));
    }
  }

  proc deinit() {
    coforall loc in Locales do on loc {
      ls_hs_destroy_spin_basis(c_ptrTo(this._localBases[loc.id]));
    }
  }

  inline proc rawPtr() { return this._localBases[here.id].payload; }
  inline proc numberSpins() { return ls_get_number_spins(this.rawPtr()); }
  inline proc hammingWeight() { return ls_get_hamming_weight(this.rawPtr()); }
  inline proc spinInversion() { return ls_get_spin_inversion(this.rawPtr()); }
  inline proc isHammingWeightFixed() { return this.hammingWeight() != -1; }

  proc statesBounds() {
    var lower: uint(64);
    var upper: uint(64);
    if (this.isHammingWeightFixed()) {
      lower = ~(0:uint(64)) >> (64 - this.hammingWeight());
      upper = lower << (this.numberSpins() - this.hammingWeight());
    }
    else {
      lower = 0;
      if (this.numberSpins() == 64) { upper = ~(0:uint(64)); }
      else { upper = (1:uint(64) << this.numberSpins()) - 1; }
    }
    if spinInversion() != 0 {
      upper = upper / 2;
    }
    assert(lower <= upper);
    return (lower, upper);
  }
}


inline proc makeShift(numberSpins : int, numberBits : int) : int {
  assert(0 < numberBits && numberBits < 64);
  return if numberBits >= numberSpins then 0 else (numberSpins - numberBits);
}

class BasisStates {
  var _size : int;
  var _representatives : [OnePerLocale] [0 ..# _size] uint(64);
  var _counts : [OnePerLocale] int;
  var _numberSpins : int;
  var _numberBits : int = 16;
  var _shift : int;
  var _ranges : [OnePerLocale] [0 ..# 1 + (1 << _numberBits)] int;

  proc _generateBucketRanges() {
    coforall loc in Locales do on loc {
      ref localRanges = _ranges[loc.id];
      ref localStates = _representatives[loc.id][0 ..# _counts[loc.id]];
      const numberBuckets = 1 << _numberBits;
      var offset = 0;
      for i in 0 ..# numberBuckets {
        localRanges[i] = offset;
        // writeln("ranges[", loc.id, "][", i, "] = ", offset);
        while offset != localStates.size && bucketIndex(localStates[offset]) == i {
          offset += 1;
        }
      }
      localRanges[numberBuckets] = offset;
      assert(offset == localStates.size);
    }
  }

  proc init(numberSpins : int, representatives, counts, numberBits : int = 16) {
    this._size = representatives[0].size;
    this._representatives = representatives;
    this._counts = counts;
    this._numberSpins = numberSpins;
    this._numberBits = min(numberSpins, numberBits);
    this._shift = makeShift(numberSpins, numberBits);
    this.complete();
    _generateBucketRanges();
  }

  inline proc totalNumberStates() : int {
    return + reduce _counts;
  }

  inline proc bucketIndex(x : uint(64)) : int {
    return (x >> _shift):int;
  }

  inline proc getCounts(loc : locale = here) : int {
    return _counts[loc.id];
    // ref r = _representatives[loc.id];
    // return r[0 ..# _counts[loc.id]];
  }

  inline proc getRepresentatives(loc : locale = here) ref {
    return _representatives[loc.id];
    // ref r = _representatives[loc.id][0 ..# _counts[loc.id]];
    // return r;
  }

  proc getIndex(x : uint(64), loc : locale = here) : int {
    const i = bucketIndex(x);
    const b = _ranges[loc.id][i];
    const e = _ranges[loc.id][i + 1];
    // The following causes a memory allocation... :(
    // ref r = getRepresentatives(loc)[_ranges[loc.id][i] .. _ranges[loc.id][i + 1] - 1];
    const dataPtr = c_ptrTo(_representatives[loc.id][b]);
    const size = (e - b):uint(64);
    var j : uint(64);
    if loc != here {
      on loc { j = ls_binary_search(dataPtr, size, x); }
    }
    else {
      j = ls_binary_search(dataPtr, size, x);
    }
    assert(j < size);
    // const (found, j) = binarySearch(r, x);
    // writeln("found=", found, ", sj=", sj, ", representatives=", representatives);
    // assert(found);
    return b + j:int;
  }
}

// proc numlocs() param where (CHPL_COMM==none) {
//   return 1;
// }
// proc numlocs() {
//   return numLocales;
// }

/*
proc datasetShape(filename : string, dataset : string) {
  const rank = ls_hs_hdf5_get_dataset_rank(filename.c_str(), dataset.c_str()):int;
  var c_shape : [0 .. rank - 1] uint(64);
  ls_hs_hdf5_get_dataset_shape(filename.c_str(), dataset.c_str(), c_ptrTo(c_shape));
  return [i in c_shape.domain] c_shape[i]:int; 
}

proc _makeDomain(shape : 1 * int) : domain(1) { return {0 .. shape[0]:int - 1}; }
proc _makeDomain(shape : 2 * int) : domain(2) {
  return {0 .. shape[0]:int - 1, 0 .. shape[1]:int - 1};
}

proc readHDF5Chunk(filename : string, dataset : string, type eltType, offset, shape) {
  const dom = _makeDomain(shape);
  var array : [dom] eltType;
  readHDF5Chunk(filename, dataset, offset, array);
  return array;
}

proc readHDF5Chunk(filename : string, dataset : string, offset, array : [] ?eltType) {
  assert(offset.size == array.rank);
  const rank = offset.size;
  var c_offset : [0 .. rank - 1] uint(64) = noinit;
  var c_shape : [0 .. rank - 1] uint(64) = noinit;
  for i in 0 .. rank - 1 {
    c_offset[i] = offset[i]:uint;
    c_shape[i] = array.dim(i).size:uint;
  }
  if (eltType == uint(64)) {
    ls_hs_hdf5_read_chunk_u64(filename.c_str(), dataset.c_str(),
      rank:c_uint, c_ptrTo(c_offset), c_ptrTo(c_shape), c_ptrTo(array));
  }
  else if (eltType == real(64)) {
    ls_hs_hdf5_read_chunk_f64(filename.c_str(), dataset.c_str(),
      rank:c_uint, c_ptrTo(c_offset), c_ptrTo(c_shape), c_ptrTo(array));
  }
  else {
    assert(false);
  }
  // return array;
}
*/

/*
config const loadArrayChunkSize : int = 10 * 1024 * 1024 / numLocales;

proc _histogramFromChunk(const ref array : [] uint(64), ref timer : Timer) {
  timer.start();
  var histogram : [LocaleSpace] int = 0;
  forall x in array with (+ reduce histogram) {
    const hash = (hash64_01(x) % numLocales:uint):int;
    histogram[hash] += 1;
  }
  timer.stop();
  return histogram;
}

proc calculateNumberStatesPerLocale(filename : string, dataset : string) {
  var _totalTimer = new Timer();
  var _readingTimer : [OnePerLocale] Timer = new Timer();
  var _histogramTimer : [OnePerLocale] Timer = new Timer();
  _totalTimer.start();

  const shape = datasetShape(filename, dataset);
  assert(shape.rank == 1);
  const totalNumberStates = shape[0];
  var histogram : [0 .. numLocales - 1] int;
  const workItems = splitIntoChunks(0, totalNumberStates, loadArrayChunkSize);
  coforall loc in Locales
      with (/*const*/ in workItems, in filename, in dataset,
            + reduce histogram) do on loc {
    for ((lo, hi), chunkIndex) in zip(workItems, 0..) do if (chunkIndex % numLocales == loc.id) {
      _readingTimer[loc.id].start();
      const size = hi - lo;
      const statesChunk = readHDF5Chunk(filename, dataset, uint(64), (lo,), (size,));
      _readingTimer[loc.id].stop();
      histogram += _histogramFromChunk(statesChunk, _histogramTimer[loc.id]);
    }
  }
  assert(+ reduce histogram == totalNumberStates);
  _totalTimer.stop();
  writeln("[Chapel] calculateNumberStatesPerLocale spent ", _readingTimer.elapsed(),
    " reading HDF5 files");
  writeln("[Chapel] calculateNumberStatesPerLocale spent ", _histogramTimer.elapsed(),
    " computing histograms");
  writeln("[Chapel] calculateNumberStatesPerLocale took ", _totalTimer.elapsed());
  return histogram;
}
*/

// MemoryInitialization; moveInitialize
// std::vector x{{1, 2, 3}};

/*
proc splitIntoChunks(lower : int, upper : int, chunkSize : int) {
  assert(upper >= lower);
  assert(chunkSize >= 1);
  const numChunks = (upper - lower + (chunkSize - 1)) / chunkSize;
  var chunks : [0 .. numChunks - 1] (int, int); // = noinit;
  var offset = lower;
  for i in 0 .. numChunks - 1 {
    const nextOffset = min(offset + chunkSize, upper);
    chunks[i] = (offset, nextOffset);
    offset = nextOffset;
  }
  return chunks;
}

proc _computeHashes(const ref states : [] uint(64), ref timer : Timer) {
  timer.start();
  var hashes : [states.domain] int = noinit;
  forall i in vectorizeOnly(states.domain) {
    hashes[i] = (hash64_01(states[i]) % numLocales:uint):int;
  }
  timer.stop();
  return hashes;
}

proc _distributeBasedOnHash(const ref states : [] uint(64), const ref array : [?D] ?eltType,
                            const ref hashes : [] int, ref timer : Timer) where (D.rank == 2) {
  timer.start();
  const numVectors = D.dim(0).size;
  const size = D.dim(1).size;
  assert(states.size == size);
  var buffer : [0 .. numLocales - 1, 0 .. size - 1, 0 .. numVectors - 1] eltType = noinit;
  var offsets : [LocaleSpace] int;
  for i in 0 .. size - 1 {
    const hash = hashes[i];
    // ref dest = buffer[hash, offsets[hash], ..];
    // const ref src = array[.., i];
    // for k in 0 .. numVectors - 1 {
    //   buffer[hash, offsets[hash], k] = array[k, i];
    // }
    // ==> 0.876606 0.888626 0.850986 0.829717
    for k in 0 .. numVectors - 1 {
      buffer[hash, offsets[hash], k] = array[k, i];
    }
    // ==> 36.7382 29.0505 25.9959 26.0736
    // buffer[hash, offsets[hash], ..] = array[.., i];
    // ==> 18.815 0.0 0.0 0.0 (for a file 2^6 times smaller)
    // dest = src;
    // 8.57173 0.0 0.0 0.0
    offsets[hash] += 1;
  }
  timer.stop();
  return (buffer, offsets);
}

proc _distributeBasedOnHash(const ref states : [] uint(64), ref timer : Timer) {
  timer.start();
  const size = states.size;
  var buffer : [0 .. numLocales - 1, 0 .. size - 1] uint(64) = noinit;
  var offsets : [LocaleSpace] int;
  for i in 0 .. size - 1 {
    const hash = (hash64_01(states[i]) % numLocales:uint):int;
    buffer[hash, offsets[hash]] = states[i];
    offsets[hash] += 1;
  }
  timer.stop();
  return (buffer, offsets);
}

proc loadArray(filename : string, basisDataset : string, arrayDataset : string,
               type eltType, const ref countPerLocale : [LocaleSpace] int) {
  var _initGlobalBufferTimer = new Timer();
  var _readTimer : [OnePerLocale] Timer = new Timer();
  var _processTimer1 : [OnePerLocale] Timer = new Timer();
  var _processTimer2 : [OnePerLocale] Timer = new Timer();
  var _copyTimer : [OnePerLocale] Timer = new Timer();

  const totalNumberStates = + reduce countPerLocale;
  const maxCountPerLocale = max reduce countPerLocale;
  const arrayDatasetShape = datasetShape(filename, arrayDataset);
  assert(arrayDatasetShape.size == 2);
  assert(arrayDatasetShape[1] == totalNumberStates);
  const numberVectors = arrayDatasetShape[0];
  // TODO: initialization of globalBuffer is quite a bit of time
  _initGlobalBufferTimer.start();
  var globalBuffer : [OnePerLocale] [0 .. numberVectors - 1, 0 .. maxCountPerLocale - 1] eltType;
  _initGlobalBufferTimer.stop();
  var globalOffsets : [LocaleSpace] int = 0;

  // writeln("globalBuffer.size = ", globalBuffer[0].size);
  // var otherTimer = new Timer();
  // otherTimer.start();
  // forall loc in Locales do on loc {
  //   globalBuffer[loc.id] = 10.0:eltType;
  // }
  // otherTimer.stop();
  // writeln("[Chapel] loadArray spent ", otherTimer.elapsed(), " filling globalBuffer");

  var writeIndex$ : atomic int = 0;
  const workItems = splitIntoChunks(0, totalNumberStates, loadArrayChunkSize / numberVectors);
  coforall loc in Locales
    with (/*const*/ in workItems,
          in filename,
          in basisDataset,
          in arrayDataset) do on loc {
    for ((lo, hi), chunkIndex) in zip(workItems, 0..) do if (chunkIndex % numLocales == loc.id) {
      _readTimer[loc.id].start();
      const size = hi - lo;
      const statesChunk = readHDF5Chunk(filename, basisDataset, uint(64), (lo,), (size,));
      const arrayChunk = readHDF5Chunk(filename, arrayDataset, eltType,
                                       (0, lo), (numberVectors, size));
      _readTimer[loc.id].stop();
      const localHashes = _computeHashes(statesChunk, _processTimer1[loc.id]);
      const (localBuffer, localOffsets) =
        _distributeBasedOnHash(statesChunk, arrayChunk, localHashes, _processTimer2[loc.id]);

      _copyTimer[loc.id].start();
      // Comm diagnostics
      // auto aggregation (only with forall loops with scalars)
      // --report-auto-aggregation
      writeIndex$.waitFor(chunkIndex);
      forall i in 0 .. numLocales - 1 {
        const _start = globalOffsets[i];
        const _size = localOffsets[i];
        for k in 0 ..# numberVectors {
          globalBuffer[i][k, _start ..# _size] = localBuffer[i, 0 ..# _size, k];
        }
        globalOffsets[i] += _size;
      }
      writeIndex$.add(1);
      _copyTimer[loc.id].stop();
    }
  }

  writeln("[Chapel] loadArray spent ", _initGlobalBufferTimer.elapsed(),
    " initializing globalBuffer");
  writeln("[Chapel] loadArray spent ", _readTimer.elapsed(), " reading from HDF5 files");
  writeln("[Chapel] loadArray spent ", _processTimer1.elapsed(), " computing hashes");
  writeln("[Chapel] loadArray spent ", _processTimer2.elapsed(), " copying stuff around");
  writeln("[Chapel] loadArray spent ", _copyTimer.elapsed(), " copying to other locales");
  return globalBuffer;
}

proc loadArraySimple(filename : string, basisDataset : string, arrayDataset : string,
                     type eltType, const ref countPerLocale : [LocaleSpace] int) {
  var _initGlobalBufferTimer = new Timer();
  var _readTimer : [OnePerLocale] Timer = new Timer();
  var _processTimer1 : [OnePerLocale] Timer = new Timer();
  var _processTimer2 : [OnePerLocale] Timer = new Timer();
  var _copyTimer : [OnePerLocale] Timer = new Timer();

  const totalNumberStates = + reduce countPerLocale;
  const maxCountPerLocale = max reduce countPerLocale;
  const arrayDatasetShape = datasetShape(filename, arrayDataset);
  assert(arrayDatasetShape.size == 2);
  assert(arrayDatasetShape[1] == totalNumberStates);
  const numberVectors = arrayDatasetShape[0];
  // TODO: initialization of globalBuffer is quite a bit of time
  _initGlobalBufferTimer.start();
  var globalBuffer : [OnePerLocale] [0 .. numberVectors - 1, 0 .. maxCountPerLocale - 1] eltType;
  _initGlobalBufferTimer.stop();

  coforall loc in Locales do on loc {
    const workItems = splitIntoChunks(0, totalNumberStates, loadArrayChunkSize / numberVectors);
    var globalOffset : int = 0;
    for (lo, hi) in workItems {
      _readTimer[loc.id].start();
      const size = hi - lo;
      const statesChunk = readHDF5Chunk(filename, basisDataset, uint(64), (lo,), (size,));
      const arrayChunk = readHDF5Chunk(filename, arrayDataset, eltType,
                                       (0, lo), (numberVectors, size));
      _readTimer[loc.id].stop();
      _processTimer1[loc.id].start();
      var mask : [statesChunk.domain] bool = noinit;
      forall i in 0 ..# size {
        mask[i] = (hash64_01(statesChunk[i]) % numLocales:uint):int == loc.id;
      }
      _processTimer1[loc.id].stop();

      _copyTimer[loc.id].start();
      forall k in 0 ..# numberVectors {
        var offset = globalOffset;
        for i in 0 ..# size {
          if (mask[i]) {
              globalBuffer[loc.id][k, offset] = arrayChunk[k, i];
              offset += 1;
          }
        }
      }
      // for i in 0 ..# size {
      //   if (mask[i]) {
      //     for k in 0 ..# numberVectors {
      //       globalBuffer[loc.id][k, globalOffset] = arrayChunk[k, i];
      //       globalOffset += 1;
      //     }
      //   }
      // }
      _copyTimer[loc.id].stop();
    }
  }

  writeln("[Chapel] loadArray spent ", _initGlobalBufferTimer.elapsed(),
    " initializing globalBuffer");
  writeln("[Chapel] loadArray spent ", _readTimer.elapsed(), " reading from HDF5 files");
  writeln("[Chapel] loadArray spent ", _processTimer1.elapsed(), " computing hashes");
  writeln("[Chapel] loadArray spent ", _processTimer2.elapsed(), " copying stuff around");
  writeln("[Chapel] loadArray spent ", _copyTimer.elapsed(), " copying to other locales");
  return globalBuffer;
}

proc loadStates(filename : string, basisDataset : string,
                const ref countPerLocale : [LocaleSpace] int) {
  var _initGlobalBufferTimer = new Timer();
  var _readTimer : [OnePerLocale] Timer = new Timer();
  var _processTimer1 : [OnePerLocale] Timer = new Timer();
  var _processTimer2 : [OnePerLocale] Timer = new Timer();
  var _copyTimer : [OnePerLocale] Timer = new Timer();

  const totalNumberStates = + reduce countPerLocale;
  const maxCountPerLocale = max reduce countPerLocale;
  const basisDatasetShape = datasetShape(filename, basisDataset);
  assert(basisDatasetShape.size == 1);
  assert(basisDatasetShape[0] == totalNumberStates);
  _initGlobalBufferTimer.start();
  var globalBuffer : [OnePerLocale] [0 .. maxCountPerLocale - 1] uint(64);
  _initGlobalBufferTimer.stop();
  var globalOffsets : [LocaleSpace] int = 0;

  var writeIndex$ : atomic int = 0;
  const workItems = splitIntoChunks(0, totalNumberStates, loadArrayChunkSize);
  coforall loc in Locales
    with (/*const*/ in workItems,
          in filename,
          in basisDataset) do on loc {
    for ((lo, hi), chunkIndex) in zip(workItems, 0..) do if (chunkIndex % numLocales == loc.id) {
      _readTimer[loc.id].start();
      const size = hi - lo;
      const statesChunk = readHDF5Chunk(filename, basisDataset, uint(64), (lo,), (size,));
      _readTimer[loc.id].stop();
      const (localBuffer, localOffsets) =
        _distributeBasedOnHash(statesChunk, _processTimer2[loc.id]);

      _copyTimer[loc.id].start();
      writeIndex$.waitFor(chunkIndex);
      forall i in 0 .. numLocales - 1 {
        const _start = globalOffsets[i];
        const _size = localOffsets[i];
        globalBuffer[i][_start ..# _size] = localBuffer[i, 0 ..# _size];
        globalOffsets[i] += _size;
      }
      writeIndex$.add(1);
      _copyTimer[loc.id].stop();
    }
  }

  writeln("[Chapel] loadArray spent ", _initGlobalBufferTimer.elapsed(),
    " initializing globalBuffer");
  writeln("[Chapel] loadArray spent ", _readTimer.elapsed(), " reading from HDF5 files");
  // writeln("[Chapel] loadArray spent ", _processTimer1.elapsed(), " computing hashes");
  writeln("[Chapel] loadArray spent ", _processTimer2.elapsed(), " copying stuff around");
  writeln("[Chapel] loadArray spent ", _copyTimer.elapsed(), " copying to other locales");
  return globalBuffer;
}
*/

proc makeStates(basis: DistributedBasis, bits : int = 16) {
  // startVdebug("makeStates");

  var timer = new Timer();
  timer.start();
  const (lower, upper) = basis.statesBounds();
  const chunkSize = max((upper - lower) / (numLocales * here.maxTaskPar * 5):uint, 1);
  const ranges = splitIntoRanges(
    lower, upper, chunkSize, basis.isHammingWeightFixed());
  timer.stop();
  writeln("[Chapel] Using chunkSize = ", chunkSize, "; constructed ",
          ranges.size, " ranges in ", timer.elapsed());
  timer.clear();
  timer.start();
  const D = {0 .. ranges.size - 1} dmapped Cyclic(startIdx=0);
  var chunks: [D] [LocaleSpace] list(uint(64));
  timer.stop();
  writeln("[Chapel] Allocated chunks in ", timer.elapsed());

  timer.clear();
  timer.start();
  coforall loc in Locales do on loc {
    var localIndices = chunks.localSubdomain(loc);
    coforall tid in 0 ..# here.maxTaskPar {
      var taskIndices = localIndices by here.maxTaskPar align (tid * numLocales + loc.id);
      for i in taskIndices {
        const (lower, upper) = ranges[i];
        ref chunk = chunks[i];
        for x in findStatesInRange(lower, upper, basis.isHammingWeightFixed(),
                                   basis.rawPtr()) {
          const hash = (hash64_01(x) % numLocales:uint):int;
          chunk[hash].append(x);
        }
      }
    }
  }
  timer.stop();
  writeln("[Chapel] Constructed all chunks in ", timer.elapsed());

  timer.clear();
  timer.start();
  var maxPerChunk = max reduce ([i in D] (
                      max reduce ([j in LocaleSpace]
                        chunks[i][j].size)));
  timer.stop();
  writeln("[Chapel] maxPerChunk = ", maxPerChunk, " in ", timer.elapsed());
  timer.clear();
  timer.start();
  var arrChunks: [D] [LocaleSpace] [0 ..# maxPerChunk] uint(64);
  var arrSizes: [D] [LocaleSpace] int;
  timer.stop();
  writeln("[Chapel] allocated arrChunks and arrSizes in ", timer.elapsed());

  timer.clear();
  timer.start();
  forall i in D {
    forall j in LocaleSpace {
      const n = chunks[i][j].size;
      arrChunks[i][j][0 ..# n] = chunks[i][j];
      arrSizes[i][j] = n;
    }
  }
  timer.stop();
  writeln("[Chapel] Converted chunks to arrays in ", timer.elapsed());

  // forall ((lower, upper), chunk) in zip(ranges, chunks) with (in nextStateFn) {
  //   assert(chunk.locale == here.locale);
  //   for x in findStatesInRange(lower, upper, nextStateFn, basis.rawPtr()) {
  //     const hash = (hash64_01(x) % numLocales:uint):int;
  //     chunk[hash].append(x);
  //   }
  // }

  timer.clear();
  timer.start();
  // const OnePerLocale = LocaleSpace dmapped Block(LocaleSpace);
  var counts: [OnePerLocale] int =
    [j in LocaleSpace] (+ reduce [i in 0..<ranges.size] chunks[i][j].size);
  const maxCount = max reduce counts;
  // For analyzing how uniform the distribution is.
  // writeln("[Chapel] Counts: ", counts);

  var states: [OnePerLocale] [0 ..# maxCount] uint(64);
  coforall loc in Locales do on loc {
    var offset: int = 0;
    for i in 0 ..# ranges.size {
      var c = arrSizes[i][loc.id]; // chunks[i][loc.id].size;
      states[loc.id][offset ..# c] = arrChunks[i][loc.id][0 ..# c];
      offset += c;
      // writeln("chunks[", i, "][", loc.id, "] = ", chunks[i][loc.id]);
    }
  }
  timer.stop();
  writeln("[Chapel] Constructed states in ", timer.elapsed());

  // stopVdebug();

  // NOTE: this creates a copy of states, right? :(
  return new BasisStates(basis.numberSpins(), states, counts, bits);
}

config const yamlPath = "/home/tom/src/spin-ed/example/heisenberg_pyrochlore_32.yaml";
  // "/home/tom/src/spin-ed/example/heisenberg_square_4x4.yaml";
  // "/home/tom/src/spin-ed/example/heisenberg_chain_4.yaml";
config const hdf5Path = "output.h5";
config const representativesDataset = "/representatives";
config const outputVectorsDataset = "/y";

config const testPath = "/home/tom/src/spin-ed/example/data/heisenberg_kagome_12.h5";

proc constructBasisRepresentatives() {
  // const basis = new shared DistributedBasis(yamlPath);
  // var timer = new Timer();
  // timer.start();
  // var states = makeStates(basis);
  // timer.stop();
  // writeln("[Chapel] makeStates took ", timer.elapsed());
  // for i in states.representatives.dim(0) {
  //   writeln("[Chapel] locale ", i, " contains ", states.counts[i], " states");
  // }

  var numberPerLocale = calculateNumberStatesPerLocale(testPath, "/basis/representatives");

  var loadArrayTimer = new Timer();
  loadArrayTimer.start();
  var vectors = loadArray(testPath, "/basis/representatives", "/hamiltonian/eigenvectors",
    real, numberPerLocale);
  loadArrayTimer.stop();
  writeln("[Chapel] loadArray took ", loadArrayTimer.elapsed());

  loadArrayTimer.clear();
  loadArrayTimer.start();
  var vectors2 = loadArraySimple(testPath, "/basis/representatives", "/hamiltonian/eigenvectors",
    real, numberPerLocale);
  loadArrayTimer.stop();
  writeln("[Chapel] loadArray took ", loadArrayTimer.elapsed());

  // var loadStatesTimer = new Timer();
  // loadStatesTimer.start();
  // var states = loadStates(testPath, "/basis/representatives", numberPerLocale);
  // loadStatesTimer.stop();
  // writeln("[Chapel] loadStates took ", loadStatesTimer.elapsed());


  // for i in LocaleSpace {
  //   writeln(vectors[i]);
  //   writeln("---------");
  // }
  // var vectors : [OnePerLocale] [0 .. 2 - 1, 0 .. states.size - 1] real(64);
  // forall i in OnePerLocale {
  //   on Locales[i] {
  //     for j in states.representatives[i].domain {
  //       vectors[i][0, j] = states.representatives[i][j]:real;
  //       vectors[i][1, j] = states.representatives[i][j]:real;
  //     }
  //   }
  // }
  // writeln(vectors[0]);
  // writeln(vectors[1]);
  // mergeVectorsBy(states.representatives, states.counts, vectors,
  //     hdf5Path, outputVectorsDataset);

  // timer.clear();
  // timer.start();
  // mergeStates(states.representatives, states.counts, hdf5Path,
  //   representativesDataset);
  // timer.stop();
  // writeln("[Chapel] mergeStates took ", timer.elapsed());
}


proc real_main() {
  constructBasisRepresentatives();
  // merge_test();
}

proc main() {
  for loc in Locales do on loc {
    ls_enable_logging();
    ls_hs_init();
  }
  real_main();
  for loc in Locales do on loc {
    ls_hs_exit();
  }
  return 0;
}
