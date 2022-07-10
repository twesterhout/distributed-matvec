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
config const kUseQueue : bool = true;

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
    const (n, basisStatesPtr, coeffsPtr, keysPtr) = batchedOperator.computeOffDiag(
        r.size, c_const_ptrTo(representatives[r.low]), c_const_ptrTo(x[r.low]));
    chpl_task_yield();
    var radixOffsets : c_array(int, 257);
    radixOneStep(n, keysPtr, radixOffsets, basisStatesPtr, coeffsPtr);
    chpl_task_yield();
    // const n = offsetsPtr[r.size];
    timer.stop();
    computeOffDiagTime.add(timer.elapsed(), memoryOrder.relaxed);

    timer.clear();
    timer.start();
    var numberToDo = numLocales;
    var isComplete : [0 ..# numLocales] bool = false;
    while numberToDo > 0 {
      const prev = numberToDo;
      for localeIdx in 0 ..# numLocales {
        if isComplete[localeIdx] then continue;

        const k = radixOffsets[localeIdx];
        const n = radixOffsets[localeIdx + 1] - k;
        const isDone = if n != 0 then queue.tryEnqueue(localeIdx, n, basisStatesPtr + k, coeffsPtr + k)
                                 else true;
        if isDone {
          isComplete[localeIdx] = true;
          numberToDo -= 1;
        }
        // else
        //   chpl_task_yield();
      }
      if numberToDo == prev then chpl_task_yield();
    }
    // staging.add(n, basisStatesPtr, coeffsPtr);
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

config const remoteBufferSize = 1000000;

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
  var putTime : atomic real;
  var onTime : atomic real;
  var processTime : atomic real;
  var processTimeHere : atomic real;
  var barrierTime1 : real;
  var barrierTime2 : real;
  var totalSendSize : atomic int;
  var initTimer = new Timer();
  initTimer.start();

  // const chunkSize = (representatives.size + numChunks - 1) / numChunks;
  // logDebug("Local dimension=", representatives.size, ", chunkSize=", chunkSize);
  var accessor = new ConcurrentAccessor(y);
  globalPtrStoreNoQueue[here.id] = (c_const_ptrTo(matrix.basis), c_ptrTo(accessor));
  allLocalesBarrier.barrier();

  // const remoteBufferSize = 1000000;
  const numTasks = here.maxTaskPar;
  const numChunks = max(numTasks * representatives.size * matrix.numberOffDiagTerms()
                          / remoteBufferSize,
                        1);
  const ranges : [0 ..# numChunks] range(int) = chunks(0 ..# representatives.size, numChunks);
  const ptrStore : [0 ..# numLocales] (c_ptr(Basis), c_ptr(ConcurrentAccessor(real(64)))) =
    globalPtrStoreNoQueue;

  var localBasisStatesBuffers : [0 ..# numLocales, 0 ..# remoteBufferSize] uint(64);
  var localCoeffsBuffers : [0 ..# numLocales, 0 ..# remoteBufferSize] complex(128);
  var remoteBuffers : [0 ..# numLocales] RemoteBuffer(complex(128)) =
    [i in 0 ..# numLocales] new RemoteBuffer(complex(128), remoteBufferSize, i);

  const _chunkSize = (representatives.size + (numChunks - 1)) / numChunks;
  var batchedOperators : [0 ..# numTasks] BatchedOperator =
    [i in 0 ..# numTasks] new BatchedOperator(matrix, _chunkSize);
  var radixOffsets : [0 ..# numTasks] c_array(int, 257);
  var taskPtrs : [0 ..# numTasks] (c_ptr(uint(64)), c_ptr(complex(128)));
  var totalCounts : [0 ..# numLocales] int;

  // var queue = new CommunicationQueue(eltType, globalPtrStore);

  // const ranges : [0 ..# numChunks] range(int) = chunks(0 ..# representatives.size, numChunks);
  initTimer.stop();
  initTime += initTimer.elapsed();

  const numIters = (numChunks + numTasks - 1) / numTasks;
  logDebug("numIters = ", numIters);
  for i in 0 ..# numIters {
    const firstChunkIdx = i * numTasks;
    const lastChunkIdx = min(firstChunkIdx + numTasks, numChunks) - 1;

    var computeTimer = new Timer();
    computeTimer.start();
    coforall chunkIdx in firstChunkIdx .. lastChunkIdx {
      const r = ranges[chunkIdx];
      const taskIdx = chunkIdx - firstChunkIdx;
      ref batchedOperator = batchedOperators[taskIdx];
      const (count, basisStatesPtr, coeffsPtr, keysPtr) = batchedOperator.computeOffDiag(
          r.size, c_const_ptrTo(representatives[r.low]), c_const_ptrTo(x[r.low]));
      radixOneStep(count, keysPtr, radixOffsets[taskIdx], basisStatesPtr, coeffsPtr);
      taskPtrs[taskIdx] = (basisStatesPtr, coeffsPtr);
    }
    computeTimer.stop();
    computeOffDiagTime.add(computeTimer.elapsed());

    var remoteTimer = new Timer();
    remoteTimer.start();
    coforall loc in Locales {
      if loc == here {
        var myTimer = new Timer();
        myTimer.start();
        const (basisPtr, remoteAccessorPtr) = ptrStore[loc.id];
        for chunkIdx in firstChunkIdx .. lastChunkIdx {
          const taskIdx = chunkIdx - firstChunkIdx;
          const k = radixOffsets[taskIdx][loc.id];
          const n = radixOffsets[taskIdx][loc.id + 1] - k;
          const (basisStatesPtr, coeffsPtr) = taskPtrs[taskIdx];

          localProcess(basisPtr, remoteAccessorPtr, basisStatesPtr + k,
                       coeffsPtr + k, n);
          chpl_task_yield();
        }
        myTimer.stop();
        processTimeHere.add(myTimer.elapsed());
      }
      else {
        var total = 0;
        for chunkIdx in firstChunkIdx .. lastChunkIdx {
          const taskIdx = chunkIdx - firstChunkIdx;
          const k = radixOffsets[taskIdx][loc.id];
          const n = radixOffsets[taskIdx][loc.id + 1] - k;
          assert(total + n <= remoteBufferSize);
          const (basisStatesPtr, coeffsPtr) = taskPtrs[taskIdx];
          c_memcpy(c_ptrTo(localBasisStatesBuffers[loc.id, total]),
                   basisStatesPtr + k, n:c_size_t * c_sizeof(uint(64)));
          c_memcpy(c_ptrTo(localCoeffsBuffers[loc.id, total]),
                   coeffsPtr + k, n:c_size_t * c_sizeof(complex(128)));
          total += n;
        }
        totalCounts[loc.id] = total;
        // chpl_task_yield();
        // logDebug("Sending ", localBasisStatesBuffers[loc.id, 0 ..# total], " to ", loc);

        var putTimer = new Timer();
        putTimer.start();
        remoteBuffers[loc.id].put(
          c_ptrTo(localBasisStatesBuffers[loc.id, 0]),
          c_ptrTo(localCoeffsBuffers[loc.id, 0]),
          total);
        putTimer.stop();
        putTime.add(putTimer.elapsed());
        totalSendSize.add(total);
        // chpl_task_yield();
      }
        // const (remoteBasis, remoteAccessor) = ptrStore[loc.id];
        // const remoteBasisStates = remoteBuffers[loc.id].basisStates;
        // const remoteCoeffs = remoteBuffers[loc.id].coeffs;
        // const remoteSize = total;
        // var onTimer = new Timer();
        // var _processTime : real;
        // var _processTimePtr = c_ptrTo(_processTime);
        // const mainLocaleIdx = here.id;
        // onTimer.start();
        // on loc {
        //   var myTimer = new Timer();
        //   myTimer.start();
        //   localProcess(remoteBasis, remoteAccessor, remoteBasisStates,
        //                remoteCoeffs, remoteSize);
        //   myTimer.stop();
        //   var myTime = myTimer.elapsed();
        //   PUT(c_ptrTo(myTime), mainLocaleIdx, _processTimePtr, c_sizeof(real));
        // }
        // onTimer.stop();
        // onTime.add(onTimer.elapsed());
        // processTime.add(_processTime);
    }
    var barrierTimer1 = new Timer();
    barrierTimer1.start();
    allLocalesBarrier.barrier();
    barrierTimer1.stop();
    coforall loc in Locales {
      const (remoteBasis, remoteAccessor) = ptrStore[loc.id];
      const remoteBasisStates = remoteBuffers[loc.id].basisStates;
      const remoteCoeffs = remoteBuffers[loc.id].coeffs;
      const remoteSize = totalCounts[loc.id];
      var onTimer = new Timer();
      var _processTime : real;
      var _processTimePtr = c_ptrTo(_processTime);
      const mainLocaleIdx = here.id;
      onTimer.start();
      on loc {
        var myTimer = new Timer();
        myTimer.start();
        localProcess(remoteBasis, remoteAccessor, remoteBasisStates,
                     remoteCoeffs, remoteSize);
        myTimer.stop();
        var myTime = myTimer.elapsed();
        PUT(c_ptrTo(myTime), mainLocaleIdx, _processTimePtr, c_sizeof(real));
      }
      onTimer.stop();
      onTime.add(onTimer.elapsed());
      processTime.add(_processTime);
    }
    var barrierTimer2 = new Timer();
    barrierTimer2.start();
    allLocalesBarrier.barrier();
    barrierTimer2.stop();

    remoteTimer.stop();
    remoteTime.add(remoteTimer.elapsed());
  }

  allLocalesBarrier.barrier();
  timer.stop();
  logDebug("localOffDiagonalNoQueue took ", timer.elapsed(), "\n",
           "  ├─ ", initTime, " in initialization\n",
           "  ├─ ", computeOffDiagTime, " in computeOffDiag (total)\n",
           "  ├─ ", radixTime, " in radixOneStep (total)\n",
           "  └─ ", remoteTime, " in remote stuff\n",
           "      ├─ ", processTimeHere, " in localProcess on here\n",
           "      ├─ ", putTime, " in remote PUTs\n",
           "      └─ ", onTime, " in remote on clauses\n",
           "          └─ ", processTime, " in localProcess on remote\n",
           " totalSendSize=", totalSendSize);
}

private proc localMatrixVector(matrix : Operator, const ref x : [] ?eltType, ref y : [] eltType,
                               const ref representatives : [] uint(64)) {
  assert(matrix.locale == here);
  assert(x.locale == here);
  assert(y.locale == here);
  assert(representatives.locale == here);
  localDiagonal(matrix, x, y, representatives);
  if kUseQueue then
    localOffDiagonal(matrix, x, y, representatives);
  else
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
