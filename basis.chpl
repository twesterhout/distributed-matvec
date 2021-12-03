use List;
use CPtr;
use SysCTypes;
use BlockDist;
use CyclicDist;
use Time;
use VisualDebug;

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

proc constructBasisRepresentatives() {
  const basis = new shared DistributedBasis(yamlPath);
  var timer = new Timer();
  timer.start();
  var states = makeStates(basis);
  timer.stop();
  writeln("[Chapel] makeStates took ", timer.elapsed());
  for i in states.representatives.dim(0) {
    writeln("[Chapel] locale ", i, " contains ", states.counts[i], " states");
  }
  timer.clear();
  timer.start();
  mergeStates(states.representatives, states.counts, hdf5Path,
    representativesDataset);
  timer.stop();
  writeln("[Chapel] mergeStates took ", timer.elapsed());
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
