import Random;
use CPtr;
use BlockDist;
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

// Creates a distributed array of basis states.
//
// First, a normal array is obtained from libplugin.so, but since our C code
// has no knowledge about locales, the array is located on `here`. We
// distribute the array to emulate the real-world scenario.
proc makeStates() {
  var dim = plugin_get_dimension():int;
  var _p: c_ptr(uint(64)) = plugin_get_basis_states();
  var _states = makeArrayFromPtr(_p, dim:uint);
  var states: [{0..<dim} dmapped Block(boundingBox={0..<dim})] uint(64) = _states;
  return states;
}

// The actual matrix-vector product y <- Ax where A comes from libplugin.
// basis is a list of basis vectors constructed using makeStates.
proc apply(basis: [] uint(64), x: [] real, y: [?D] real) {
  coforall L in Locales do on L {
    var dim = plugin_get_max_nonzero_per_row();
    var statesBuffer: [{0..<dim}] uint(64);
    var coeffsBuffer: [{0..<dim}] complex(128);

    for i in D.localSubdomain() {
      var currentSpin = basis[i];
      var written: uint;
      written = plugin_apply_operator(currentSpin, c_ptrTo(statesBuffer[0]), c_ptrTo(coeffsBuffer[0]));
      var acc: complex(128);
      for _k in {0..<written} {
        var j: int;
        // TODO: officially one can only call get_index on the corresponding
        // Locale. In other words, we will need another function "get_locale"
        // which will be called before get_index.
        j = plugin_get_index(statesBuffer[_k]):int;
        acc += coeffsBuffer[_k] * x[j];
      }
      assert(acc.im == 0);
      y[i] = acc.re;
    }
  }
}

// Initialize the C library.
coforall L in Locales do on L {
  plugin_init(n:uint(32));
}

var states = makeStates();
writeln("Basis states: ", states);
writeln("Max non-zero per row: ", plugin_get_max_nonzero_per_row());

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
apply(states, x1, y1);
writeln(x1);
writeln(y1);

// Verify results
assert(y1 == y2);

// Deinitialize the C library.
coforall L in Locales do on L {
  plugin_deinit();
}
