use matvec;
use basis;
use wrapper;
use Distribute;
use Merge;
use MatVec.IO;
use CPtr;

config const kInputBasisPath = "data/heisenberg_chain_10.yaml";
config const kInputDataPath = "data/matvec/heisenberg_chain_10.h5";
config const kOutputDataPath = "output.h5";
config const kDebugLatticeSymmetries = false;

proc main() {
  initRuntime(false);

  var basis = new DistributedBasis(kInputBasisPath);
  var states = makeStates(basis);
  var matrix = new DistributedOperator(kInputBasisPath);

  var mask = distributionMask(states, chunkSize = 10);
  var X = distributeVectors(loadVectors(kInputDataPath, "/x")[0 ..# 1, ..], mask);
  var Y : [OnePerLocale] [X[0].domain] real(64);
  matvecSerial(matrix, states, X, Y);
  writeln(Y);

  // coforall loc in Locales do on loc {
  //   writeln(matrix.requiredBufferSize());
  // }

  // if (false) {
  //   const batchSize = 1;
  //   const bufferSize = 13;
  //   coforall loc in Locales do on loc {
  //     var spins : [0 ..# 1] uint(64) = 0;
  //     var offsetsBatch : [0 ..# batchSize + 1] uint(64);
  //     var spinsBatch : [0 ..# bufferSize] uint(64);
  //     var coeffsBatch : [0 ..# bufferSize] complex(128);
  //     matrix.batchedApply(
  //       states.representatives[loc.id][0 ..# 1],
  //       offsetsBatch,
  //       spinsBatch,
  //       coeffsBatch);
  //   }
  // }

  var YExpected = distributeVectors(loadVectors(kInputDataPath, "/y")[0 ..# 1, ..], mask);
  writeln(YExpected);

  var status : int = 0;
  for loc in Locales {
    const equal = && reduce (Y[loc.id] == YExpected[loc.id]);
    if !equal {
      writeln("[Chapel] ERROR: Y and YExpected do not match!"); 
      status = 1;
    }
  }
  return status;
}

// proc deinit() {
//   deinitRuntime();
// }
