import Random;
use BitOps;
use List;
use CPtr;
use SysCTypes;
use BlockDist;
use CyclicDist;
use CommDiagnostics;

use states;

config const n = 4;

extern proc plugin_init(n:uint(32));
extern proc plugin_deinit();
extern proc plugin_get_dimension(): uint(64);
extern proc plugin_get_max_nonzero_per_row(): uint(64);
extern proc plugin_get_basis_states(): c_ptr(uint(64));
extern proc plugin_apply_operator(spin:uint(64), other_spins:c_ptr(uint(64)), other_coeffs:c_ptr(complex(128))): uint(64);
extern proc plugin_get_index(spin:uint(64)): uint(64);
extern proc plugin_matvec(x: c_ptr(real(64)), y: c_ptr(real(64)));

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
  // writeln(makeStates());

  var states = makeStatesChapel();
  for i in states.domain {
    writeln(states[i]);
  }
  // var buckets = shuffleWithHash(states);
  // for i in LocaleSpace {
  //   writeln(buckets[i]);
  // }

  coforall L in Locales do on L { plugin_deinit(); }
}

testStates();
