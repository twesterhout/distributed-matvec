module DistributedMatrixVector {

use CTypes;
use RangeChunk;
use AllLocalesBarriers;
use Time;
use DynamicIters;

use FFI;
use ForeignTypes;
use ConcurrentAccessor;
use BatchedOperator;
use CommunicationQueue;

/* 
 */
private proc localDiagonalBatch(indices : range(int, BoundedRangeType.bounded, false),
                                matrix : Operator, const ref x : [] ?eltType, ref y : [] eltType,
                                const ref representatives : [] uint(64)) {
  const batchSize = indices.size;
  // assert(workspace.size >= batchSize);
  // if workspace.size < batchSize then
  //   workspace.domain = {0 ..# batchSize};
  ls_internal_operator_apply_diag_x1(
    matrix.payload, batchSize, c_const_ptrTo(representatives[indices.low]),
    c_ptrTo(y[indices.low]), c_const_ptrTo(x[indices.low]));
  // ls_hs_operator_apply_diag_kernel(
  //   matrix.payload, batchSize,
  //   c_const_ptrTo(representatives[indices.low]), 1,
  //   c_ptrTo(workspace));
  // foreach i in indices {
  //   y[i] = x[i] * workspace[i - indices.low]:eltType;
  // }
}

config const matrixVectorDiagonalNumChunks : int = 10 * here.maxTaskPar;
config const matrixVectorOffDiagonalNumChunks : int = 150 * here.maxTaskPar;
config const matrixVectorMainLoopNumTasks : int = here.maxTaskPar;

private proc localDiagonal(matrix : Operator, const ref x : [] ?eltType, ref y : [] eltType,
                           const ref representatives : [] uint(64),
                           numChunks : int = min(matrixVectorDiagonalNumChunks,
                                                 representatives.size)) {
  const totalSize = representatives.size;
  const batchSize = (totalSize + numChunks - 1) / numChunks;
  var ranges : [0 ..# numChunks] range(int, BoundedRangeType.bounded, false) =
    chunks(0 ..# totalSize, numChunks);
  // var workspace : [0 ..# batchSize] complex(128) = noinit;
  forall r in ranges {
    localDiagonalBatch(r, matrix, x, y, representatives);
  }
}

var globalPtrStore : [LocaleSpace] (c_ptr(Basis), c_ptr(ConcurrentAccessor(real(64))));

private proc localOffDiagonal(matrix : Operator, const ref x : [] ?eltType, ref y : [] eltType,
                              const ref representatives : [] uint(64),
                              numChunks : int = min(matrixVectorOffDiagonalNumChunks,
                                                    representatives.size)) {
  var timer = new Timer();
  timer.start();
  var initTime : real;
  var computeOffDiagTime : atomic real;
  var stagingAddTime : atomic real;
  var stagingFlushTime : atomic real;
  var queueDrainTime : real;
  var initTimer = new Timer();
  initTimer.start();

  const chunkSize = (representatives.size + numChunks - 1) / numChunks;
  logDebug("Local dimension=", representatives.size, ", chunkSize=", chunkSize);
  // Used to calculate the action of matrix on a chunk of basis elements.
  // batchedOperator stores some task-local buffers and a pointer to `matrix`.
  // var batchedOperator = new BatchedOperator(matrix, chunkSize);
  // Used to then process the matrix elements.
  // logDebug("new ConcurrentAccessor: " + y.locale:string);
  var accessor = new ConcurrentAccessor(y);
  globalPtrStore[here.id] = (c_const_ptrTo(matrix.basis), c_ptrTo(accessor));
  allLocalesBarrier.barrier();

  var queue = new CommunicationQueue(eltType, globalPtrStore);

  const ranges : [0 ..# numChunks] range(int) =
    chunks(0 ..# representatives.size, numChunks);
  initTimer.stop();
  initTime += initTimer.elapsed();
  forall rangeIdx in dynamic(0 ..# numChunks, chunkSize=1,
                             numTasks=matrixVectorMainLoopNumTasks)
              // ranges
                     with (ref queue,
                           var batchedOperator = new BatchedOperator(matrix, chunkSize),
                           var staging = new StagingBuffers(queue)) {
    const r : range(int) = ranges[rangeIdx];
    // logDebug("Processing ", r, " ...");

    var timer = new Timer();
    timer.start();
    const (basisStatesPtr, coeffsPtr, offsetsPtr) = batchedOperator.computeOffDiag(
        r.size, c_const_ptrTo(representatives[r.low]), c_const_ptrTo(x[r.low]));
    const n = offsetsPtr[r.size];
    timer.stop();
    computeOffDiagTime.add(timer.elapsed(), memoryOrder.relaxed);

    timer.clear();
    timer.start();
    staging.add(n, basisStatesPtr, coeffsPtr);
    timer.stop();
    stagingAddTime.add(timer.elapsed(), memoryOrder.relaxed);

    // timer.clear();
    // timer.start();
    // staging.flush();
    // timer.stop();
    // stagingFlushTime.add(timer.elapsed(), memoryOrder.relaxed);
  }

  // logDebug("Draining ...");
  var queueTimer = new Timer();
  queueTimer.start();
  queue.drain();
  queueTimer.stop();
  queueDrainTime += queueTimer.elapsed();

  allLocalesBarrier.barrier();
  timer.stop();
  logDebug("localOffDiagonal took ", timer.elapsed(), "\n",
           "  ├─ ", initTime, " in initialization\n",
           "  ├─ ", computeOffDiagTime, " in computeOffDiag (total)\n",
           "  ├─ ", stagingAddTime, " in staging.add (total)\n",
           "  │   └─ ", queue.enqueueTime, " in queue.enqueue\n",
           "  │       ├─ ", queue.localProcessTimeHere, " in localProcess on here\n",
           "  │       ├─ ", queue.lockWaitingTimeHere, " waiting for queue.lock\n",
           "  │       └─ ", queue.enqueueUnsafeTime, " in queue._enqueueUnsafe\n",
           "  │           └─ ", queue.flushBufferTime, " in queue._flushBuffer\n",
           "  │               ├─ ", queue.lockWaitingTimeRemote, " waiting for remoteBuffer.lock\n",
           "  │               ├─ ", queue.flushBufferPutTime, " in remote PUTs\n",
           "  │               └─ ", queue.flushBufferRemoteTime, " in remote tasks\n",
           "  │                   └─ ", queue.localProcessTimeRemote, " in localProcess on remote\n",
           "  └─ ", queueDrainTime, " in queue.drain\n");
  // logDebug("Spent ", queue.flushBufferRemoteTime, " in queue._flushBuffer (remote part)");
  // logDebug("Spent ", queue.enqueueUnsafeTime, " in queue._enqueueUnsafe");
  // logDebug("Spent ", queue.enqueueTime, " in queue.enqueue");
  // logDebug("Spent ", queue.localProcessTime, " in localProcess");
  // logDebug("Spent ", globalLocalProcessTime[here.id][0].read(), " in localProcess (all, indexing)");
  // logDebug("Spent ", globalLocalProcessTime[here.id][1].read(), " in localProcess (all, allocate)");
  // logDebug("Spent ", globalLocalProcessTime[here.id][2].read(), " in localProcess (all, accessing)");
  // logDebug("Spent ", queue._remoteBuffers.localProcessTime, " in remoteLocalProcess");
  // logDebug("Frequency ", frequency.read() / norm.read());
  // logDebug("Processed ", queue!.numRemoteCalls.read(), " remote jobs");
}

private proc localMatrixVector(matrix : Operator, const ref x : [] ?eltType, ref y : [] eltType,
                               const ref representatives : [] uint(64)) {
  assert(matrix.locale == here);
  assert(x.locale == here);
  assert(y.locale == here);
  assert(representatives.locale == here);
  localDiagonal(matrix, x, y, representatives);
  localOffDiagonal(matrix, x, y, representatives);
}

proc matrixVectorProduct(matrixFilename : string,
                         const ref x,
                         ref y,
                         const ref representatives) {
  coforall loc in Locales with (ref y) do on loc {
    var myMatrix = loadHamiltonianFromYaml(matrixFilename);
    const ref myX = x.getBlock(loc.id)[0, ..];
    ref myY = y.getBlock(loc.id)[0, ..];
    const ref myBasisStates = representatives.getBlock(loc.id);
    myMatrix.basis.uncheckedSetRepresentatives(myBasisStates);
    // myY += 1;
    var timer = new Timer();
    timer.start();
    localMatrixVector(myMatrix, myX, myY, myBasisStates);
    timer.stop();
    logDebug("Spent ", timer.elapsed(), " in localMatrixVector");
  }
}

/*
proc localMatrixVectorPart(
      H : Operator,
      const ref x : [] ?eltType,
      ref y : [] eltType,
      const ref representatives : [] uint(64)) {
  y = 0;

  // Diagonal contribution
  forall (sigma, i) in zip(representatives, 0..) {

  }


  // Off-diagonal contribution

  coforall loc in Locales do on loc {
    var accessor = new unmanaged ConcurrentAccessor();

    for i,  in zip
    
  }


  for i in 0 ..# x.size {
    // representatives is a distributed vector of basis states. It is pre-computed
    const |σᵢ⟩ = representatives[i];
    for (cⱼ, |σⱼ⟩) in H|σᵢ⟩ {
      const t = x[i] * cⱼ;
      const j = indexOf(|σⱼ⟩);
      y[j] += t;
    }
  }

  y = 0; // initialize y with zeros
  for i in 0 ..# x.size {
    // representatives is a distributed vector of basis states. It is pre-computed
    const |σᵢ⟩ = representatives[i];
    for (cⱼ, |σⱼ⟩) in H|σᵢ⟩ {
      const t = x[i] * cⱼ;
      const j = indexOf(|σⱼ⟩);
      y[j] += t;
    }
  }
}
*/

}
