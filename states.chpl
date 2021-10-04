use BitOps;
use List;
use CPtr;
use SysCTypes;
use BlockDist;
use CyclicDist;

// Functions which we need to construct a list of "representatives" (i.e.
// special binary configurations with which we're going to work)
// extern proc plugin_get_number_spins(cxt: c_void_ptr): c_uint;
// extern proc plugin_get_hamming_weight(cxt: c_void_ptr): c_int;
// extern proc plugin_get_spin_inversion(cxt: c_void_ptr): c_int;
// extern proc plugin_is_representative(
//   cxt: c_void_ptr, count: uint(64), spins: c_ptr(uint(64)), output: c_ptr(uint(8))): c_int;

// extern proc plugin_new_options_example_01(number_spins: c_uint): c_void_ptr;
// extern proc plugin_delete_options_example_01(cxt: c_void_ptr);

require "lattice_symmetries/lattice_symmetries.h";
extern type ls_error_code;
extern type ls_flat_spin_basis;

extern proc ls_flat_spin_basis_number_spins(basis: c_ptr(ls_flat_spin_basis)): c_uint;
extern proc ls_flat_spin_basis_hamming_weight(basis: c_ptr(ls_flat_spin_basis)): c_int;
extern proc ls_flat_spin_basis_spin_inversion(basis: c_ptr(ls_flat_spin_basis)): c_int;
// 
// void ls_flat_spin_basis_state_info(ls_flat_spin_basis const* basis, uint64_t count,
//                                    void const* spin, void* repr,
//                                    LATTICE_SYMMETRIES_COMPLEX128* character, double* norm);
extern proc ls_flat_spin_basis_is_representative(basis: c_ptr(ls_flat_spin_basis), count: uint(64),
                                                 spin: c_void_ptr, is_repr: c_ptr(uint(8)), norm: c_ptr(real));


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

/* Get next integer with the same hamming weight.

   Semantically equivalent to
   ```
   const m = popcount(v);
   v += 1;
   while (popcount(v) != m) { v += 1; }
   return v;
   ```
 */
proc nextStateFixedHamming(v: uint(64)): uint(64) {
  const t = v | (v - 1);
  return (t + 1) | (((~t & (t + 1)) - 1) >> (ctz(v) + 1));
}

proc nextStateGeneral(v: uint(64)): uint(64) {
  return v + 1;
}

/* Fill buffer using either nextStateFixedHamming or nextStateGeneral. Returns
   the number of elements written.
 */
proc manyNextState(in v: uint(64), bound: uint(64), buffer: [] uint(64),
                   nextStateFn: func(uint(64), uint(64))) {
  assert(v <= bound);
  for i in buffer.domain {
    buffer[i] = v;
    // If v == bound, incrementing v further is unsafe because it can overflow.
    if (v == bound) { return i + 1; }
    v = nextStateFn(v);
  }
  return buffer.size;
}

/* C code uses SIMD and it's beneficial to call plugin_is_representative with
   count > 1. isRepresentativeBatchSize controls the size of chunks which will
   be fed to plugin_is_representative.
 */
config const isRepresentativeBatchSize = 1;

iter findStatesInRange(in lower: uint(64), upper: uint(64), nextStateFn: func(uint(64), uint(64)),
                       basis: c_ptr(ls_flat_spin_basis)) {
  const chunkSize = isRepresentativeBatchSize;
  var buffer: [0..<chunkSize] uint(64);
  var flags: [0..<chunkSize] uint(8);
  var norms: [0..<chunkSize] real(64);
  while (true) {
    const written = manyNextState(lower, upper, buffer, nextStateFn);
    ls_flat_spin_basis_is_representative(
        basis, written:uint, c_ptrTo(buffer):c_void_ptr, c_ptrTo(flags), c_ptrTo(norms));
    for i in {0..<written} {
      if (flags[i]:bool) { yield buffer[i]; }
    }
    if (buffer[written - 1] == upper) { break; }
    lower = nextStateFn(lower);
  }
}

proc testBit(bits: uint(64), i: int): bool { return ((bits >> i) & 1):bool; }
proc clearBit(bits: uint(64), i: int) { return bits & ~(1:uint(64) << i); }
proc setBit(bits: uint(64), i: int) { return bits | (1:uint(64) << i); }

/* Get the closest to `x` integer with Hamming weight `hammingWeight`.
 */
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

proc splitIntoRanges(in current: uint(64), bound: uint(64), chunkSize: uint(64), isHammingWeightFixed: bool) {
  var ranges: list((uint(64), uint(64)));
  // var isHammingWeightFixed = plugin_get_hamming_weight() != -1;
  var hammingWeight = popcount(current);
  while (true) {
      if (bound - current <= chunkSize) {
          ranges.append((current, bound));
          break;
      }
      var next = current + chunkSize;
      if isHammingWeightFixed { next = closestWithFixedHamming(next, hammingWeight); }

      assert(next >= current);
      if (next >= bound) {
          ranges.append((current, bound));
          break;
      }
      ranges.append((current, next));
      current = if isHammingWeightFixed then nextStateFixedHamming(next) else nextStateGeneral(next);
  }

  var distRanges = newCyclicArr({0..<ranges.size}, (uint(64), uint(64)));
  distRanges = ranges;
  return distRanges;
}

record PerLocaleState {
  var size: int;
  var representatives: [0..<size] uint(64);

  proc init(size: int) {
    this.size = size;
  }
}

proc makeStates(bases: [] c_ptr(ls_flat_spin_basis)) {
  var numberSpins = ls_flat_spin_basis_number_spins(bases[here.id]);
  var hammingWeight = ls_flat_spin_basis_hamming_weight(bases[here.id]);
  var spinInversion = ls_flat_spin_basis_spin_inversion(bases[here.id]);
  var isHammingWeightFixed = hammingWeight != -1;
  var lower: uint(64);
  var upper: uint(64);
  if (isHammingWeightFixed) {
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
  var ranges = splitIntoRanges(lower, upper, chunkSize, isHammingWeightFixed);
  const D = {0..<ranges.size} dmapped Cyclic(startIdx=0);
  var chunks: [D] [LocaleSpace] list(uint(64));

  var nextStateFn: func(uint(64), uint(64));
  if (isHammingWeightFixed) { nextStateFn = nextStateFixedHamming; }
  else { nextStateFn = nextStateGeneral; }
  forall ((lower, upper), chunk) in zip(ranges, chunks) with (in nextStateFn) {
    for x in findStatesInRange(lower, upper, nextStateFn, bases[here.id]) {
      var hash = (hash64_01(x) % numLocales:uint):int;
      chunk[hash].append(x);
    }
  }

  const OnePerLocale = LocaleSpace dmapped Block(LocaleSpace);
  var counts: [OnePerLocale] int =
    [j in LocaleSpace] (+ reduce [i in 0..<ranges.size] chunks[i][j].size);
  var maxCount = max reduce counts;
  writeln("Counts: ", counts); // For analyzing how uniform the distribution is.

  var states: [OnePerLocale] PerLocaleState =
    [i in OnePerLocale] new PerLocaleState(counts[i]);

  // var states: [OnePerLocale] [0..<maxCount] uint(64);
  coforall loc in Locales do on loc {
    var offset: int = 0;
    for i in {0..<ranges.size} {
      var c = chunks[i][loc.id].size;
      states[loc.id].representatives[offset..<offset + c] = chunks[i][loc.id];
      offset += c;
    }
  }

  return states;



  // var perLocaleStates: [OnePerLocale] PerLocaleState;
  //   this.localStates = [i in OnePerLocale] new PerLocaleState(counts[i]);
  //   coforall loc in Locales do on loc {
  //     this.localStates[loc.id].representatives = states[loc.id][0..<counts[loc.id]];
  //   }
  // return (states, counts);
}


// export proc foo(k: c_uint): opaque {
//   const OnePerLocale = LocaleSpace dmapped Block(LocaleSpace);
//   var contexts: [OnePerLocale] c_void_ptr;
//   for loc in Locales do on loc { contexts[loc.id] = plugin_new_options_example_01(k); }
//   var (states, counts) = makeStates(contexts);
//   for loc in Locales do on loc { plugin_delete_options_example_01(contexts[loc.id]); }
// 
//   var localStates: [OnePerLocale] PerLocaleState = [i in OnePerLocale] new PerLocaleState(counts[i]);
//   coforall loc in Locales do on loc {
//     localStates[loc.id].representatives = states[loc.id][0..<counts[loc.id]];
//   }
// 
//   // var ret = convertToExternalArray(localStates);
//   return localStates;
// }


// Creates a distributed array of basis states.
//
// First, a normal array is obtained from libplugin.so, but since our C code
// has no knowledge about locales, the array is located on `here`. We
// distribute the array to emulate the real-world scenario.
// proc makeStates() {
//   var dim = plugin_get_dimension():int;
//   var _p: c_ptr(uint(64)) = plugin_get_basis_states();
//   var _states = makeArrayFromPtr(_p, dim:uint);
//   const D = {0..<dim};
//   var states: [D dmapped Block(boundingBox={0..<dim})] uint(64) = _states;
//   return states;
// }
