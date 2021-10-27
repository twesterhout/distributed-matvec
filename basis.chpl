use List;
use CPtr;
use SysCTypes;
use BlockDist;
use CyclicDist;
use Time;

use wrapper;
use states;

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

  proc rawPtr() { return this._localBases[here.id].payload; }
  proc numberSpins() { return ls_get_number_spins(this.rawPtr()); }
  proc hammingWeight() { return ls_get_hamming_weight(this.rawPtr()); }
  proc isHammingWeightFixed() { return this.hammingWeight() != -1; }

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

/* Local chunk of representative spin configurations.
 */
class LocalBasisStates {
  var size: int;
  var representatives: [0..<size] uint(64);

  proc init(size: int) { this.size = size; }
}

// proc numlocs() param where (CHPL_COMM==none) {
//   return 1;
// }
// proc numlocs() {
//   return numLocales;
// }

proc makeStates(basis: DistributedBasis) {
  const (lower, upper) = basis.statesBounds();
  const chunkSize = max((upper - lower) / 1000, 1);
  // writeln("[Chapel] Using chunkSize = ", chunkSize, "...");
  const ranges = splitIntoRanges(
    lower, upper, chunkSize, basis.isHammingWeightFixed());
  const D = {0..<ranges.size} dmapped Cyclic(startIdx=0);
  var chunks: [D] [LocaleSpace] list(uint(64));

  const nextStateFn: func(uint(64), uint(64));
  if (basis.isHammingWeightFixed()) { nextStateFn = nextStateFixedHamming; }
  else { nextStateFn = nextStateGeneral; }
  // writeln("[Chapel] Calling main forall...");
  // writeln("[Chapel] bases = ", bases);
  forall ((lower, upper), chunk) in zip(ranges, chunks) with (in nextStateFn) {
    for x in findStatesInRange(lower, upper, nextStateFn, basis.rawPtr()) {
      const hash = (hash64_01(x) % numLocales:uint):int;
      chunk[hash].append(x);
    }
  }
  // writeln("[Chapel] Constructed all chunks");

  // const OnePerLocale = LocaleSpace dmapped Block(LocaleSpace);
  var counts: [OnePerLocale] int =
    [j in LocaleSpace] (+ reduce [i in 0..<ranges.size] chunks[i][j].size);
  // var maxCount = max reduce counts;
  writeln("[Chapel] Counts: ", counts); // For analyzing how uniform the distribution is.

  var states: [OnePerLocale] owned LocalBasisStates =
    [i in OnePerLocale] new LocalBasisStates(counts[i]);
  coforall loc in Locales do on loc {
    var offset: int = 0;
    for i in {0..<ranges.size} {
      var c = chunks[i][loc.id].size;
      states[loc.id].representatives[offset..<offset + c] = chunks[i][loc.id];
      offset += c;
    }
  }
  return states;
}

config const yamlPath = "/home/tom/src/spin-ed/example/heisenberg_chain_4.yaml";
      // "/home/tom/src/spin-ed/example/heisenberg_square_4x4.yaml"

proc real_main() {
  const basis = new shared DistributedBasis(yamlPath);
  var timer = new Timer();
  timer.start();
  var states = makeStates(basis); // new DistributedStates(context);
  timer.stop();
  writeln("[Chapel] makeStates took ", timer.elapsed());

  writeln(states[0].representatives);
  writeln(states[1].representatives);
  writeln(merge(states[0].representatives, states[1].representatives));
  writeln(merge(states[1].representatives, states[0].representatives));
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
