use matvec;
use basis;
use wrapper;
use Distribute;
use Merge;
use MatVec.IO;
use CPtr;
use Time;

config const kInputBasisPath = "data/heisenberg_chain_10.yaml";
config const kInputDataPath = "data/matvec/heisenberg_chain_10.h5";
config const kOutputDataPath = "output.h5";
config const kDebugLatticeSymmetries = false;
config const kDistributionChunkSize = 10000;
config const kPrintOutput = false;

proc main() {
  initRuntime(false);
  var timer = new Timer();
  var basis = new DistributedBasis(kInputBasisPath);
  var states = makeStates(basis);
  var matrix = new DistributedOperator(kInputBasisPath);

  timer.start();
  var mask = distributionMask(states, chunkSize = kDistributionChunkSize);
  timer.stop();
  writeln("[Chapel] distributionMask took ", timer.elapsed());

  timer.clear();
  timer.start();
  var X = distributeVectors(loadVectors(kInputDataPath, "/x")[0 ..# 1, ..], mask);
  timer.stop();
  writeln("[Chapel] loadVectors+distributeVectors took ", timer.elapsed());

  timer.clear();
  timer.start();
  var Y : [OnePerLocale] [X[0].domain] real(64);
  writeln("[Chapel] initializing Y took ", timer.elapsed());
  timer.stop();

  matvecSerial(matrix, states, X, Y);
  if (kPrintOutput) {
    writeln(Y);
  }

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

  timer.clear();
  timer.start();
  var YExpected = distributeVectors(loadVectors(kInputDataPath, "/y")[0 ..# 1, ..], mask);
  timer.stop();
  writeln("[Chapel] loadVectors+distributeVectors took ", timer.elapsed());
  if kPrintOutput {
    writeln(YExpected);
  }

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
