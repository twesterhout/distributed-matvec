module DistributedMatrixVector {

use CTypes;
use RangeChunk;
use AllLocalesBarriers;
use Time;
use DynamicIters;
use CommDiagnostics;

use FFI;
use ForeignTypes;
use ConcurrentAccessor;
use BatchedOperator;
use CommunicationQueue;

config const kVerboseComm : bool = false; 
config const kVerboseGetTiming : bool = false; 

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
  }

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
}

record PartitionInfo {
    var _countOrOffset : int;
    var nextOffset : int;

    inline proc count ref { return _countOrOffset; }
    inline proc offset ref { return _countOrOffset; }
};

proc partitionBy(in first : c_ptr(?eltType), last : c_ptr(eltType), predicate) {
    while true {
      if first == last then return last;
      if !predicate(first.deref()) then break;
      first += 1;
    }

    var it = first + 1;
    while it != last {
      if predicate(it.deref()) {
        first.deref() <=> it.deref();
        first += 1;
      }
      it += 1;
    }
    return first;
}

inline proc swapElements(a : int, b : int, arr : c_ptr(?t1)) {
  arr[a] <=> arr[b];
}
inline proc swapElements(a : int, b : int, arr1 : c_ptr(?t1), arr2 : c_ptr(?t2)) {
  swapElements(a, b, arr1);
  swapElements(a, b, arr2);
}

proc radixOneStep(numKeys : int, keys : c_ptr(uint(8)), offsets : c_array(int, 257), arrs...?numArrs)
{
  var partitions : c_array(PartitionInfo, 256);
  foreach i in 0 ..# numKeys {
    partitions[keys[i]:int].count += 1;
  }

  var remainingPartitions : c_array(uint(8), 256);
  var numPartitions : int;
  var total : int;
  for i in 0 ..# 256 {
    const count = partitions[i].count;
    if count > 0 {
      partitions[i].offset = total;
      total += count;
      remainingPartitions[numPartitions] = i:uint(8);
      numPartitions += 1;
    }
    partitions[i].nextOffset = total;
  }

  var lastRemaining = remainingPartitions:c_ptr(uint(8)) + numPartitions;
  var endPartition = remainingPartitions:c_ptr(uint(8)) + 1;
  while lastRemaining - endPartition > 0 {
    record Func {
      inline proc this(partitionIdx : uint(8)) {
        ref beginOffset = partitions[partitionIdx:int].offset;
        ref endOffset = partitions[partitionIdx:int].nextOffset;
        if beginOffset == endOffset then return false;

        for i in beginOffset .. endOffset - 1 {
          ref offset = partitions[keys[i]:int].offset;
          keys[i] <=> keys[offset];
          swapElements(i, offset, (...arrs));
          offset += 1;
        }
        return beginOffset != endOffset;
      }
    }
    lastRemaining = partitionBy(remainingPartitions:c_ptr(uint(8)), lastRemaining, new Func());
  }

  offsets[0] = 0;
  foreach i in 1 ..# 256 {
    offsets[i] = partitions[i - 1].nextOffset;
  }
}


var globalPtrStoreNoQueue : [LocaleSpace] (c_ptr(Basis), c_ptr(ConcurrentAccessor(real(64))));

private proc localOffDiagonalNoQueue(matrix : Operator, const ref x : [] ?eltType, ref y : [] eltType,
                                     const ref representatives : [] uint(64),
                                     numChunks : int = min(matrixVectorOffDiagonalNumChunks,
                                                           representatives.size)) {
  var timer = new Timer();
  timer.start();
  var initTime : real;
  var computeOffDiagTime : atomic real;
  var radixTime : atomic real;
  var remoteTime : atomic real;
  var allocateTime : atomic real;
  var getTime : atomic real;
  var processTime : atomic real;
  var initTimer = new Timer();
  initTimer.start();

  const chunkSize = (representatives.size + numChunks - 1) / numChunks;
  // logDebug("Local dimension=", representatives.size, ", chunkSize=", chunkSize);
  var accessor = new ConcurrentAccessor(y);
  globalPtrStoreNoQueue[here.id] = (c_const_ptrTo(matrix.basis), c_ptrTo(accessor));
  allLocalesBarrier.barrier();

  const ptrStore : [0 ..# numLocales] (c_ptr(Basis), c_ptr(ConcurrentAccessor(real(64)))) =
    globalPtrStoreNoQueue;
  // var queue = new CommunicationQueue(eltType, globalPtrStore);

  const ranges : [0 ..# numChunks] range(int) = chunks(0 ..# representatives.size, numChunks);
  initTimer.stop();
  initTime += initTimer.elapsed();
  forall rangeIdx in dynamic(0 ..# numChunks, chunkSize=1,
                             numTasks=0) // matrixVectorMainLoopNumTasks)
      with (var batchedOperator = new BatchedOperator(matrix, chunkSize)) {
    const r : range(int) = ranges[rangeIdx];

    var timer = new Timer();
    timer.start();
    const (count, basisStatesPtr, coeffsPtr, keysPtr) = batchedOperator.computeOffDiag(
        r.size, c_const_ptrTo(representatives[r.low]), c_const_ptrTo(x[r.low]));
    // for i in 0 ..# count {
    //   logDebug(rangeIdx, ": ", i, ": ", basisStatesPtr[i], ", ", coeffsPtr[i]);
    // }
    // const n = offsetsPtr[r.size];
    timer.stop();
    computeOffDiagTime.add(timer.elapsed(), memoryOrder.relaxed);
    for i in 0 ..# numLocales {
      chpl_task_yield();
    }

    timer.clear();
    timer.start();
    var offsets : c_array(int, 257);
    radixOneStep(count, keysPtr, offsets, basisStatesPtr, coeffsPtr);
    timer.stop();
    radixTime.add(timer.elapsed(), memoryOrder.relaxed);
    for i in 0 ..# numLocales {
      chpl_task_yield();
    }

    for localeIdx in 0 ..# numLocales {
      const myCount = offsets[localeIdx + 1] - offsets[localeIdx];
      if myCount != 0 {
        const mainBasisStatesPtr = basisStatesPtr + offsets[localeIdx];
        const mainCoeffsPtr = coeffsPtr + offsets[localeIdx];
        const mainLocaleIdx = here.id;
        const (myBasisPtr, myAccessorPtr) = ptrStore[localeIdx];
        var mainTimers : c_array(real, 3);
        const mainTimersPtr = mainTimers:c_ptr(real);
        if localeIdx == mainLocaleIdx {
          localProcess(myBasisPtr, myAccessorPtr,
                       mainBasisStatesPtr,
                       mainCoeffsPtr,
                       myCount);
          chpl_task_yield();
        }
        else {
          if kVerboseComm then startVerboseCommHere();
          timer.clear();
          timer.start();
          on Locales[localeIdx] {
            var myAllocateTimer = new Timer();
            var myGetTimer = new Timer();
            var myProcessTimer = new Timer();
            // const (myBasisPtr, myAccessorPtr) = globalPtrStoreNoQueue[localeIdx];
            myAllocateTimer.start();
            var myBasisStates : [0 ..# myCount] uint(64);
            var myCoeffs : [0 ..# myCount] complex(128);
            myAllocateTimer.stop();

            myGetTimer.start();
            GET(c_ptrTo(myBasisStates[0]), mainLocaleIdx, mainBasisStatesPtr,
                myCount:c_size_t * c_sizeof(uint(64)));
            myGetTimer.stop();
            if kVerboseGetTiming then logDebug("GET from locale ", mainLocaleIdx, " of ", myCount:c_size_t * c_sizeof(uint(64)),
                                               " bytes took ", myGetTimer.elapsed(), ", i.e. ",
                                               (myCount:c_size_t * c_sizeof(uint(64))):real / myGetTimer.elapsed()
                                                 / 1024.0 / 1024.0 / 1024.0, " GB/s");
            myGetTimer.start();
            GET(c_ptrTo(myCoeffs[0]), mainLocaleIdx, mainCoeffsPtr,
                myCount:c_size_t * c_sizeof(myCoeffs.eltType));
            myGetTimer.stop();
            // logDebug("myBasisStates=", myBasisStates);
            // logDebug("myBasisStates=", myCoeffs);
            myProcessTimer.start();
            localProcess(myBasisPtr, myAccessorPtr,
                         c_ptrTo(myBasisStates[0]),
                         c_ptrTo(myCoeffs[0]),
                         myCount);
            myProcessTimer.stop();

            var myTimers : c_array(real, 3);
            myTimers[0] = myAllocateTimer.elapsed();
            myTimers[1] = myGetTimer.elapsed();
            myTimers[2] = myProcessTimer.elapsed();
            PUT(myTimers:c_ptr(real), mainLocaleIdx, mainTimersPtr,
                myTimers.size:c_size_t * c_sizeof(real));
          }
          timer.stop();
          remoteTime.add(timer.elapsed(), memoryOrder.relaxed);
          allocateTime.add(mainTimers[0], memoryOrder.relaxed);
          getTime.add(mainTimers[1], memoryOrder.relaxed);
          processTime.add(mainTimers[2], memoryOrder.relaxed);
          if kVerboseComm then stopVerboseCommHere();
        }
      }
    }
    // timer.clear();
    // timer.start();
    // staging.add(n, basisStatesPtr, coeffsPtr);
    // timer.stop();
    // stagingAddTime.add(timer.elapsed(), memoryOrder.relaxed);
  }

  // var queueTimer = new Timer();
  // queueTimer.start();
  // queue.drain();
  // queueTimer.stop();
  // queueDrainTime += queueTimer.elapsed();

  allLocalesBarrier.barrier();
  timer.stop();
  logDebug("localOffDiagonalNoQueue took ", timer.elapsed(), "\n",
           "  ├─ ", initTime, " in initialization\n",
           "  ├─ ", computeOffDiagTime, " in computeOffDiag (total)\n",
           "  ├─ ", radixTime, " in radixOneStep (total)\n",
           "  └─ ", remoteTime, " in remote on clauses (total)\n",
           "      ├─ ", allocateTime, " in allocations\n",
           "      ├─ ", getTime, " in remote GETs\n",
           "      └─ ", processTime, " in localProcess");
}

private proc localMatrixVector(matrix : Operator, const ref x : [] ?eltType, ref y : [] eltType,
                               const ref representatives : [] uint(64)) {
  assert(matrix.locale == here);
  assert(x.locale == here);
  assert(y.locale == here);
  assert(representatives.locale == here);
  localDiagonal(matrix, x, y, representatives);
  // localOffDiagonal(matrix, x, y, representatives);
  localOffDiagonalNoQueue(matrix, x, y, representatives);
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
