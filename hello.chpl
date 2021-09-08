import Random;
use CPtr;
use LinearAlgebra;
use BlockDist, ReplicatedDist;
use CommDiagnostics;

config const n = 4;
config const p = 0.6;

extern proc plugin_init(n:uint(32));
extern proc plugin_deinit();
extern proc plugin_get_dimension():uint(64);
extern proc plugin_get_max_nonzero_per_row():uint(64);
extern proc plugin_get_basis_states(): c_ptr(uint(64));
extern proc plugin_apply_operator(spin:uint(64), other_spins:c_ptr(uint(64)), other_coeffs:c_ptr(complex(128))): uint(64);


/* This class will be swapped by external C functions in the future.
 */
record CSR {
  var numberNonZero: int;
  var dimension: int;
  var elements: [{0..<numberNonZero}] real; // [{0..<numberNonZero} dmapped Replicated()] real;
  var columnIndices: [{0..<numberNonZero}] int; // [{0..<numberNonZero} dmapped Replicated()] int;
  var rowIndices: [{0..dimension}] int; // [{0..dimension} dmapped Replicated()] int;

  proc init(elements, columnIndices, rowIndices) {
    this.numberNonZero = elements.size;
    this.dimension = rowIndices.size - 1;
    this.complete();
    // for L in Locales do on L {
    on Locales[0] {
      this.elements = elements;
      this.columnIndices = columnIndices;
      this.rowIndices = rowIndices;
    }
  }
}

/* For testing purposes only, converts dense matrices to CSR format.

   NOTE: why is it crazy slow?
 */
proc dense2csr(A: [?D] real) {
  writeln("Building CSR...");
  var dimension = D.dim(0).size;
  var numberNonZero = + reduce (A != 0);
  var elements: [0..<numberNonZero] real;
  var columnIndices: [0..<numberNonZero] int;
  var rowIndices: [0..dimension] int;

  var offset = 0;
  rowIndices[0] = 0;
  for i in D.dim(0) {
    rowIndices[i + 1] = rowIndices[i];
    for j in D.dim(1) {
      if A[i, j] != 0 {
        elements[offset] = A[i, j];
        columnIndices[offset] = j;
        offset += 1;
        rowIndices[i + 1] += 1;
      }
    }
  }

  writeln("Constructing CSR...");
  return new CSR(elements, columnIndices, rowIndices);
}

proc computeRow(matrix: CSR, i: int) {
  var start = matrix.rowIndices[i];
  var end = matrix.rowIndices[i + 1];
  var cs = matrix.elements[start..<end];
  var js = matrix.columnIndices[start..<end];
  return (cs, js);
}

proc computeLocale(j: int, _x) {
  return _x[j].locale;
}

proc computeIndex(j: int) {
  return j;
}

proc matvec(H: CSR, x: [?D1] real, y: [?D2] real) {
  writeln("matvec...");
  coforall L in Locales {
    on L {
      var acc: [0..<numLocales] real;
      for i in D2.localSubdomain() {
        acc = 0;

        var (coeffs, basisElements) = computeRow(H, i);
        coforall (c, b) in zip(coeffs, basisElements) {
          on computeLocale(b, x) {
            var j = computeIndex(b);
            acc[here.id] += c * x[j];
          }
        }
        y[i] = + reduce acc;
      }
    }
  }
  // forall i in D2 {
  //   var (coeffs, basisElements) = computeRow(H, i);
  //   var js = computeIndex(basisElements);
  //   var r: real;
  //   for (c, j) in zip(coeffs, js) {
  //     writeln("Accessing j=", j:string, " on ", x[j].locale);
  //     r += c * x[j];
  //   }
  //   y[i] = r;
  // }
}

plugin_init(n:uint(32));
var dim = plugin_get_dimension():int;
writeln("Dimension: ", dim);
writeln("Max non-zero per row: ", plugin_get_max_nonzero_per_row());

proc makeStates() {
  var _p: c_ptr(uint(64)) = plugin_get_basis_states();
  var _states = makeArrayFromPtr(_p, dim:uint);
  var states: [{0..<dim} dmapped Block(boundingBox={0..<dim})] uint(64) = _states;
  return states;
}

var states = makeStates();
writeln(states);

var otherStatesBuffer: [{0..<plugin_get_max_nonzero_per_row()}] uint(64);
var otherCoeffsBuffer: [{0..<plugin_get_max_nonzero_per_row()}] complex(128);
var written: uint(64);
var status: int(32);

for x in states {
  written = plugin_apply_operator(x, c_ptrTo(otherStatesBuffer[0]), c_ptrTo(otherCoeffsBuffer[0]));
  writeln(written);
  writeln(otherStatesBuffer[0..<written]);
  writeln(otherCoeffsBuffer[0..<written]);
}

plugin_deinit();

/*
const D1 = {0..<n};
const D2 = {0..<n, 0..<n};
var x: [D1 dmapped Block(boundingBox=D1)] real;
var x2: [D1] real;
var y: [D1 dmapped Block(boundingBox=D1)] real;
var H: [D2] real;

Random.fillRandom(H, 1234);
Random.fillRandom(x, 1235); // this is not optimal, right?

x2 = x;

for h in H do
  if abs(h) < p then h = 0;

writeln(H);
writeln(dot(H, x2));

var implicit = dense2csr(H);
matvec(implicit, x, y);
writeln("Printing...");
writeln(y);



// var A: [D1 dmapped Replicated()] real;
// var B: [{0..<numLocales} dmapped Block(boundingBox={0..<numLocales})] real;

for L in Locales {
  on L {
    writeln("On locale ", here.id, " elements=", implicit.elements); // .replicand(L));
  }
}
*/
