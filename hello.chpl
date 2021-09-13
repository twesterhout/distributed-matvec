import Random;
use BitOps;
use List;
use CPtr;
use SysCTypes;
use BlockDist;
use CyclicDist;
use CommDiagnostics;

config const n = 4;

extern proc plugin_init(n:uint(32));
extern proc plugin_deinit();
extern proc plugin_get_dimension(): uint(64);
extern proc plugin_get_max_nonzero_per_row(): uint(64);
extern proc plugin_get_basis_states(): c_ptr(uint(64));
extern proc plugin_apply_operator(spin:uint(64), other_spins:c_ptr(uint(64)), other_coeffs:c_ptr(complex(128))): uint(64);
extern proc plugin_get_index(spin:uint(64)): uint(64);
extern proc plugin_matvec(x: c_ptr(real(64)), y: c_ptr(real(64)));

// Functions which we need to construct a list of "representatives" (i.e.
// special binary configurations with which we're going to work)
extern proc plugin_get_number_spins(): c_uint;
extern proc plugin_get_hamming_weight(): c_int;
extern proc plugin_get_spin_inversion(): c_int;
extern proc plugin_is_representative(count: uint(64), spins: c_ptr(uint(64)), output: c_ptr(uint(8))): c_int;

proc nextStateFixedHamming(v: uint(64)): uint(64) {
  var t = v | (v - 1);
  return (t + 1) | (((~t & (t + 1)) - 1) >> (ctz(v) + 1));
}

proc nextStateGeneral(v: uint(64)): uint(64) {
  return v + 1;
}

proc findRepresentativesInRange(lower: uint(64), upper: uint(64)) {
  var m = plugin_get_hamming_weight();
  if (m != -1) {
    assert(popcount(lower) == m);
    assert(popcount(upper) == m);
  }
  var representatives: list(uint(64));
  var x = lower;
  var flag: uint(8);
  while (x < upper) {
    plugin_is_representative(1, c_ptrTo(x), c_ptrTo(flag));
    if (flag:bool) { representatives.append(x); }

    if (m != -1) { x = nextStateFixedHamming(x); }
    else { x = nextStateGeneral(x); }
  }
  assert(x == upper);
  plugin_is_representative(1, c_ptrTo(x), c_ptrTo(flag));
  if (flag:bool) { representatives.append(x); }

  return representatives;
}

proc testBit(bits: uint(64), i: int): bool { return ((bits >> i) & 1):bool; }
proc clearBit(bits: uint(64), i: int) { return bits & ~(1:uint(64) << i); }
proc setBit(bits: uint(64), i: int) { return bits | (1:uint(64) << i); }

proc closestWithFixedHamming(in x: uint(64), hammingWeight: uint): uint(64) {
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

proc splitIntoTasks(in current: uint(64), bound: uint(64), chunkSize: uint(64)) {
  var ranges: list((uint(64), uint(64)));
  var isHammingWeightFixed = plugin_get_hamming_weight() != -1;
  var hammingWeight = popcount(current);
  while (true) {
      if (bound - current <= chunkSize) {
          ranges.append((current, bound));
          break;
      }
      var next: uint(64);
      if (isHammingWeightFixed) { next = closestWithFixedHamming(current + chunkSize, hammingWeight); }
      else { next = current + chunkSize; }

      assert(next >= current);
      if (next >= bound) {
          ranges.append((current, bound));
          break;
      }
      ranges.append((current, next));
      if (isHammingWeightFixed) { current = nextStateFixedHamming(next); }
      else { current = nextStateGeneral(next); }
  }

  var distRanges = newCyclicArr({0..<ranges.size}, (uint(64), uint(64)));
  distRanges = ranges;
  return distRanges;
}

proc processTasks(ranges: [] (uint(64), uint(64))) {
  var states: [0..<ranges.localSubdomain().size] list(uint(64));
  coforall (i, _k) in zip(0.., ranges.localSubdomain()) {
    var (lower, upper) = ranges[_k];
    states[i] = findRepresentativesInRange(lower, upper);
  }
  var combined: list(uint(64));
  for part in states {
    combined.extend(part);
  }
  return combined;
}

proc makeStatesChapel() {
  var numberSpins = plugin_get_number_spins();
  var hammingWeight = plugin_get_hamming_weight();
  var isFixedHammingWeight = hammingWeight != -1;
  var lower: uint(64);
  var upper: uint(64);
  if (isFixedHammingWeight) {
    lower = ~(0:uint(64)) >> (64 - hammingWeight);
    upper = lower << (numberSpins - hammingWeight);
  }
  else {
    lower = 0;
    if (numberSpins == 64) { upper = ~(0:uint(64)); }
    else { upper = (1:uint(64) << numberSpins) - 1; }
  }

  var chunkSize = max((upper - lower) / 1000, 1);
  writeln("Using chunkSize=", chunkSize, "...");
  var ranges = splitIntoTasks(lower, upper, chunkSize);

  var states: [LocaleSpace dmapped Block(LocaleSpace)] list(uint(64));
  coforall loc in Locales do on loc {
    states[loc.id] = processTasks(ranges);
  }
  return states;
}


// Creates a distributed array of basis states.
//
// First, a normal array is obtained from libplugin.so, but since our C code
// has no knowledge about locales, the array is located on `here`. We
// distribute the array to emulate the real-world scenario.
proc makeStates() {
  var dim = plugin_get_dimension():int;
  var _p: c_ptr(uint(64)) = plugin_get_basis_states();
  var _states = makeArrayFromPtr(_p, dim:uint);
  const D = {0..<dim};
  var states: [D dmapped Block(boundingBox={0..<dim})] uint(64) = _states;
  return states;
}

// 2^N --> binom(N, N/2)

proc makeIndexMap(states: [?D] uint(64)) {
  var ranges: LocaleSpace (uint(64), uint(64));
  coforall L in Locales do on L {
    var indices = D.localSubdomain().dim(0);
    ranges[L.id] = (states[indices.first], states[indices.last]);
  }
  return ranges;
}

proc getLocale(ranges: [] (uint(64), uint(64)), spin: uint(64)) {
  for (i, (l, u)) in zip(0.., ranges) {
    if (l <= spin && spin <= u) { return i; }
  }
  assert(false);
  return -1;
}

proc hash64_01(in x: uint(64)): uint(64) {
  x = (x ^ (x >> 30)) * (0xbf58476d1ce4e5b9:uint(64));
  x = (x ^ (x >> 27)) * (0x94d049bb133111eb:uint(64));
  x = x ^ (x >> 31);
  return x;
}

proc hash32_01(in x: uint(32)): uint(32) {
  x += 1;
  x ^= x >> 17;
  x *= 0xed5ad4bb:uint(32);
  x ^= x >> 11;
  x *= 0xac4c1b51:uint(32);
  x ^= x >> 15;
  x *= 0x31848bab:uint(32);
  x ^= x >> 14;
  return x;
}

proc analyzeStatesDistribution(states: [] uint(64)) {
  var histogram: [LocaleSpace] uint(64);
  for x in states {
    // var l = hash64_01(x) % numLocales:uint;
    var l = hash32_01(x:uint(32)) % (numLocales:uint(32));
    histogram[l:int] += 1;
  }
  return histogram;
}

record LocaleState {
  var size: uint;
  var basis: [0..<size] uint;
  var x: [0..<size] real;
  var y: [0..<size] real;
}



// The actual matrix-vector product y <- Ax where A comes from libplugin.
// basis is a list of basis vectors constructed using makeStates.

// {spins} 0101 -> 0

// 0: 011      0
// 1: 010 110  0 1
//               x
// 010 -> 1     0.5
// 011 -> 0      -2
// 110 -> 1      10
// Replicated
// 
// 011 -> 0      -2
// 010 -> 1     0.5
// 110 -> 1      10

// 011 ___     -2
// 010 110     .5 10
//
// getLocale :: 3bits -> {0, 1} uint(64)
proc apply(basis: [] uint(64), ranges: [] (uint(64), uint(64)), x: [] real, y: [?D] real) {
  coforall L in Locales do on L {
    var dim = plugin_get_max_nonzero_per_row();
    var statesBuffer: [{0..<dim}] uint(64);
    var coeffsBuffer: [{0..<dim}] complex(128);

    //
    // var x: int;
    // forall i in 1..n with (in x) { ... }
    forall i in D.localSubdomain() with (in statesBuffer, in coeffsBuffer) {
      var currentSpin = basis[i];
      var written = plugin_apply_operator(currentSpin, c_ptrTo(statesBuffer[0]), c_ptrTo(coeffsBuffer[0]));
      var acc: complex(128);
      for _k in {0..<written} {
        var localeIndex = getLocale(ranges, statesBuffer[_k]);
        var xj: complex(128);
        on Locales[localeIndex] {
          var j = plugin_get_index(statesBuffer[_k]):int;
          xj = x[j];
        }
        acc += coeffsBuffer[_k] * xj;
      }
      assert(acc.im == 0);
      y[i] = acc.re;
    }
  }
}

proc testHash() {
  // Initialize the C library.
  coforall L in Locales do on L {
    plugin_init(n:uint(32));
  }

  var states = makeStates();
  writeln(analyzeStatesDistribution(states));

  // Deinitialize the C library.
  coforall L in Locales do on L {
    plugin_deinit();
  }
}

proc testApply() {
  // Initialize the C library.
  coforall L in Locales do on L {
    plugin_init(n:uint(32));
  }

  var states = makeStates();
  writeln("Basis states: ", states);
  writeln("Max non-zero per row: ", plugin_get_max_nonzero_per_row());

  var ranges = makeIndexMap(states);

  // Construct x and y arrays. x1 and y1 are used by Chapel code and x2 and y2
  // are used by libplugin. y1 should afterwards equal y2.
  const D = {0..<plugin_get_dimension():int};
  var x1: [D dmapped Block(boundingBox=D)] real;
  var y1: [D dmapped Block(boundingBox=D)] real;
  var x2: [D] real;
  var y2: [D] real;

  // Initialize x
  Random.fillRandom(x1, 1235);
  x2 = x1;

  // Perform y <- Ax in C
  plugin_matvec(c_ptrTo(x2[0]), c_ptrTo(y2[0]));
  writeln(x2);
  writeln(y2);

  // Perform y <- Ax in Chapel
  apply(states, ranges, x1, y1);
  writeln(x1);
  writeln(y1);

  // Verify results
  assert(y1 == y2);

  // Deinitialize the C library.
  coforall L in Locales do on L {
    plugin_deinit();
  }
}

proc testStates() {
  coforall L in Locales do on L { plugin_init(n:uint(32)); }

  // var k = plugin_get_number_spins();
  // var m = plugin_get_hamming_weight();
  // var lower = 0xFFFFFFFFFFFFFFFF >> (64 - m);
  // var upper = lower << (k - m);
  // var states = findRepresentativesInRange(lower, upper);
  // writeln(states);
  writeln(makeStates());

  writeln(makeStatesChapel());


  coforall L in Locales do on L { plugin_deinit(); }
}

testStates();
