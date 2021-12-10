use List;
use CPtr;
use SysCTypes;
use BlockDist;
use CyclicDist;
use Time;
use VisualDebug;
use RangeChunk;

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
      const localPath = path;
      var dummyHamiltonian = new ls_hs_operator_v1(nil, nil);
      ls_hs_basis_and_hamiltonian_from_yaml(
        localPath.c_str(), c_ptrTo(this._localBases[loc.id]), c_ptrTo(dummyHamiltonian));
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
    return (lower, upper);
  }
};

class BasisStates {
  var size : int;
  var representatives : [OnePerLocale] [0 ..# size] uint(64);
  var counts : [OnePerLocale] int;
}

// proc numlocs() param where (CHPL_COMM==none) {
//   return 1;
// }
// proc numlocs() {
//   return numLocales;
// }

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
  return readHDF5Chunk(filename, dataset, offset, array);
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
  return array;
}


record _StatesChunk { 
  var lock$ : atomic bool;
  var dom : domain(1);
  var buffer : [dom] uint(64);
  var size : int;

  proc init(chunkSize : int) {
    this.dom = {0 .. chunkSize - 1};
    this.size = chunkSize;
  }

  proc lock() { while this.lock$.testAndSet() do chpl_task_yield(); }
  proc unlock() { this.lock$.clear(); }
  proc isLocked() { return this.lock$.read(); }
}

proc histogramFromChunk(const ref array : [] uint(64)) {
  var histogram : [LocaleSpace] int = 0;
  forall x in array with (+ reduce histogram) {
    const hash = (hash64_01(x) % numLocales:uint):int;
    histogram[hash] += 1;
  }
  return histogram;
}
proc histogramFromChunk(const ref chunk : _StatesChunk) {
  assert(chunk.isLocked());
  const needCopy = chunk.buffer.locale != here.locale;
  if (needCopy) {
    var localBuffer = chunk.buffer[0 .. chunk.size - 1];
    return histogramFromChunk(localBuffer);
  }
  return histogramFromChunk(chunk.buffer[0 .. chunk.size - 1]);
}

proc lockNextChunk(const ref chunks : [?D] _StatesChunk) {
  var iteration = 0;
  while (true) {
    for i in D {
      if (!chunks[i].isLocked()) {
        chunks[i].lock();
        return i;
      }
    }
    iteration += 1;
    if (iteration >= 5) {
      sleep(1, TimeUnits.milliseconds);
    }
    else {
      chpl_task_yield();
    }
  }
  halt("unreachable");
}

config const readingChunkSize = 1024 * 1024;

proc calculateNumberStatesPerLocale(filename : string, dataset : string) {
  var timer = new Timer();
  var readTimer = new Timer();
  var timeReading : real = 0;
  timer.start();
  const shape = datasetShape(filename, dataset);
  assert(shape.rank == 1);
  const totalNumberStates = shape[0];
  const chunkSize = min(readingChunkSize, totalNumberStates);

  var chunks : [LocaleSpace] _StatesChunk = new _StatesChunk(chunkSize);
  for chunk in chunks { assert(!chunks.isLocked()); }
  var histogram : [0 .. numLocales - 1, 0 .. numLocales - 1] int;
  var offset = 0;
  while (offset < totalNumberStates) {
    const count = min(chunkSize, totalNumberStates - offset);
    const localeIndex = lockNextChunk(chunks);
    ref chunk = chunks[localeIndex];

    // readTimer.clear();
    readTimer.start();
    readHDF5Chunk(filename, dataset, (offset,), chunk.buffer[0 .. count - 1]);
    readTimer.stop();
    // timeReading += readTimer.elapsed();

    chunk.size = count;
    begin with (ref chunk) {
      on Locales[localeIndex] {
        histogram[localeIndex, ..] += histogramFromChunk(chunk);
      }
      chunk.unlock();
    }
    offset += count;
  }
  for chunk in chunks { chunk.lock(); chunk.unlock(); } // wait for all tasks to complete

  assert(offset == totalNumberStates);
  var finalHistogram : [LocaleSpace] int;
  for j in LocaleSpace {
    finalHistogram[j] = + reduce [i in LocaleSpace] histogram[i, j];
  }
  writeln(finalHistogram);
  assert(+ reduce finalHistogram == totalNumberStates);
  timer.stop();
  writeln("[Chapel] calculateNumberStatesPerLocale took ", timer.elapsed());
  writeln("[Chapel] timeReading ", readTimer.elapsed());
  return finalHistogram;
}

proc loadArray(filename : string, basisDataset : string, arrayDataset : string,
               type eltType, countPerLocale : [LocaleSpace] int) {
  const totalNumberStates = + reduce countPerLocale;
  const maxCountPerLocale = max reduce countPerLocale;
  const arrayDatasetShape = datasetShape(filename, arrayDataset);
  assert(arrayDatasetShape.size == 2);
  assert(arrayDatasetShape[1] == totalNumberStates);
  const numberVectors = arrayDatasetShape[0];

  var globalBuffer : [OnePerLocale] [0 .. numberVectors - 1, 0 .. maxCountPerLocale - 1] eltType;
  var globalOffsets : [LocaleSpace] int = 0;

  const numChunks = max(1, totalNumberStates / (numLocales * numLocales));
  for r in chunks(0 ..< totalNumberStates, numChunks) {
    const statesChunk = readHDF5Chunk(filename, basisDataset, uint(64),
                                      (r.first,), (r.size,));
    const arrayChunk = readHDF5Chunk(filename, arrayDataset, eltType,
                                     (0, r.first), (numberVectors, r.size));
    var localBuffer : [0 .. numLocales - 1, 0 .. numberVectors - 1, 0 .. r.size - 1] eltType;
    var localOffsets : [LocaleSpace] int;
    for i in 0 .. r.size - 1 {
      const hash = (hash64_01(statesChunk[i]) % numLocales:uint):int;
      localBuffer[hash, .., localOffsets[hash]] = arrayChunk[.., i];
      localOffsets[hash] += 1;
    }
    for i in 0 .. numLocales - 1 {
      const start = globalOffsets[i];
      const size = localOffsets[i];
      globalBuffer[i][.., start ..# size] = localBuffer[i, .., 0 ..# size];
      globalOffsets[i] += size;
    }
  }

  return globalBuffer;
}


proc makeStates(basis: DistributedBasis) {
  // startVdebug("makeStates");

  var timer = new Timer();
  timer.start();
  const (lower, upper) = basis.statesBounds();
  const chunkSize = max((upper - lower) / 1000, 1);
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
  const nextStateFn: func(uint(64), uint(64));
  if (basis.isHammingWeightFixed()) { nextStateFn = nextStateFixedHamming; }
  else { nextStateFn = nextStateGeneral; }
  forall ((lower, upper), chunk) in zip(ranges, chunks) with (in nextStateFn) {
    assert(chunk.locale == here.locale);
    for x in findStatesInRange(lower, upper, nextStateFn, basis.rawPtr()) {
      const hash = (hash64_01(x) % numLocales:uint):int;
      chunk[hash].append(x);
    }
  }
  timer.stop();
  writeln("[Chapel] Constructed all chunks in ", timer.elapsed());

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
    for i in {0..<ranges.size} {
      var c = chunks[i][loc.id].size;
      states[loc.id][offset .. offset + c - 1] = chunks[i][loc.id];
      offset += c;
    }
  }
  timer.stop();
  writeln("[Chapel] Constructed states in ", timer.elapsed());

  // stopVdebug();

  // NOTE: this creates a copy of states, right? :(
  return new BasisStates(maxCount, states, counts);
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
  var vectors = loadArray(testPath, "/basis/representatives", "/hamiltonian/eigenvectors",
    real, numberPerLocale);

  for i in LocaleSpace {
    writeln(vectors[i]);
    writeln("---------");
  }
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
