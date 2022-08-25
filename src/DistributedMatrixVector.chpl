module DistributedMatrixVector {

use CTypes;
use RangeChunk;
use AllLocalesBarriers;
use Time;
use DynamicIters;
use CommDiagnostics;
use ChapelLocks;

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
config const matrixVectorOffDiagonalNumChunks : int = 32 * here.maxTaskPar;
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

record _LocalBuffer {
  type coeffType;
  var destLocaleIdx : int;
  var srcLocaleIdx : int;
  var capacity : int;
  var size : int;
  var basisStates : c_ptr(uint(64));
  var coeffs : c_ptr(coeffType);
  var isFull : c_ptr(chpl__processorAtomicType(bool));
  var isEmpty : chpl__processorAtomicType(bool);
  var isEOF : chpl__processorAtomicType(bool);

  proc postinit() {
    isEmpty.write(true);
    isEOF.write(false);
    assert(destLocaleIdx == here.id);
    // isFull.store(false);
    // const rvf_size = size;
    // on Locales[destLocaleIdx] {
      basisStates = c_malloc(uint(64), capacity);
      coeffs = c_malloc(coeffType, capacity);
      // isEmpty.deref().store(true);
    // }
  }

  proc deinit() {
    // assert(!isFull.read());
    // if basisStates != nil {
      // assert(coeffs != nil);
      // const rvf_basisStates = basisStates;
      // const rvf_coeffs = coeffs;
      // on Locales[destLocaleIdx] {
        c_free(basisStates);
        c_free(coeffs);
      // }
    // }
  }
}

record _RemoteBuffer {
  type coeffType;
  var destLocaleIdx : int;
  var srcLocaleIdx : int;
  var capacity : int;
  var size : c_ptr(int);
  var basisStates : c_ptr(uint(64));
  var coeffs : c_ptr(coeffType);
  var isFull : chpl__processorAtomicType(bool);
  var isEmpty : c_ptr(chpl__processorAtomicType(bool));
  var isEOF : c_ptr(chpl__processorAtomicType(bool));

  var _putTime : real;
  var _waitTime : real;
  var _submitTime : real;

  proc postinit() {
    isFull.write(false);
  }

  proc put(localBasisStates : c_ptr(uint(64)),
           localCoeffs : c_ptr(coeffType),
           count : int) {
    var timer = new Timer();
    timer.start();

    assert(here.id == srcLocaleIdx);
    assert(count <= capacity);
    assert(!isFull.read());
    PUT(localBasisStates, destLocaleIdx, basisStates, count:c_size_t * c_sizeof(uint(64)));
    PUT(localCoeffs, destLocaleIdx, coeffs, count:c_size_t * c_sizeof(coeffType));
    // PUT(c_const_ptrTo(count), destLocaleIdx, size, c_sizeof(int));

    timer.stop();
    _putTime += timer.elapsed();
  }

  inline proc put(localBasisStates : [] uint(64),
                  localCoeffs : [] coeffType,
                  size : int) {
    put(c_const_ptrTo(localBasisStates[0]), c_const_ptrTo(localCoeffs[0]), size);
  }

  proc submit(basisStatesPtr : c_ptr(uint(64)),
              coeffsPtr : c_ptr(coeffType),
              count : int,
              inTheMeantime) {
    var timer = new Timer();
    timer.start();

    var waitTimer = new Timer();
    waitTimer.start();
    // Wait for the buffer to become empty
    while isFull.read() {
      inTheMeantime();
      chpl_task_yield();
      // tryProcessLocal(taskIdx, newLocalBuffers, matrix, accessor);
    }
    waitTimer.stop();
    _waitTime += waitTimer.elapsed();

    put(basisStatesPtr, coeffsPtr, count);
    isFull.write(true);
    const atomicPtr = isEmpty;
    const sizePtr = size;
    on Locales[destLocaleIdx] {
      sizePtr.deref() = count;
      atomicStoreBool(atomicPtr, false);
    }

    timer.stop();
    _submitTime += timer.elapsed();
  }

  proc submit(basisStatesPtr : c_ptr(uint(64)),
              coeffsPtr : c_ptr(coeffType),
              count : int) {
    record InTheMeantime {
      inline proc this() {
        return;
      }
    }
    submit(basisStatesPtr, coeffsPtr, count, new InTheMeantime());
  }

  proc finish() {
    // startVerboseCommHere();
    const atomicPtr = isEOF;
    on Locales[destLocaleIdx] {
      atomicStoreBool(atomicPtr, true);
      // pragma "local fn" pragma "fast-on safe extern function"
      // extern proc atomic_store_bool(ref obj : chpl__processorAtomicType(bool), value : bool) : void;

      // atomic_store_bool(atomicPtr.deref(), true);
      // isEOF.deref().write(true);
    }
    // stopVerboseCommHere();
  }
}






var globalPtrStoreNoQueue : [LocaleSpace] (c_ptr(Basis),
                                           c_ptr(ConcurrentAccessor(real(64))),
                                           c_ptr(_RemoteBuffer(complex(128))),
                                           c_ptr(_LocalBuffer(complex(128))),
                                           int);

config const kRemoteBufferSize = 150000;
config const kNumTasks = here.maxTaskPar;
config const kNumConsumerTasks = 1;
config const kVerbose = false;
config const specialCase = false;

extern proc chpl_task_getId(): chpl_taskID_t;

inline proc atomicStoreBool(p : c_ptr(chpl__processorAtomicType(bool)), value : bool) {
  pragma "local fn" pragma "fast-on safe extern function"
  extern proc atomic_store_bool(ref _obj : chpl__processorAtomicType(bool),
                                _value : bool) : void;

  atomic_store_bool(p.deref(), value);
}

proc tryProcessLocal(taskIdx : int, srcLocaleIdx : int, newLocalBuffers, matrix, accessor) {
  ref localBuffer = newLocalBuffers[srcLocaleIdx, taskIdx];
  if !localBuffer.isEmpty.read() {
    // if kVerbose then
    //   logDebug("335: Calling localProcess for newLocalBuffers[", srcLocaleIdx, ", ", taskIdx, "]...");
    localProcess(c_const_ptrTo(matrix.basis), c_const_ptrTo(accessor),
                 localBuffer.basisStates,
                 localBuffer.coeffs,
                 localBuffer.size);
    localBuffer.isEmpty.write(true);
    const atomicPtr = localBuffer.isFull;
    on Locales[localBuffer.srcLocaleIdx] {
      atomicStoreBool(atomicPtr, false);
    }
    // if kVerbose then
    //   logDebug("335: done with localProcess for newLocalBuffers[", srcLocaleIdx, ", ", taskIdx, "]");
    return true;
  }
  return false;
}
proc tryProcessLocal(taskIdx : int, newLocalBuffers, matrix, accessor) {
  var hasDoneWork : bool = false;
  for srcLocaleIdx in 0 ..# numLocales {
    hasDoneWork = hasDoneWork ||
      tryProcessLocal(taskIdx, srcLocaleIdx, newLocalBuffers, matrix, accessor);
  }
  return hasDoneWork;
}

record LocalOffDiagonalTimers {
  var total : Timer;
  var initialization : Timer;
  var computeOffDiag : atomic real;
  var submit : real;
  var wait : real;
  var put : real;
}

proc _offDiagMakeLocalBuffers(numTasks : int, remoteBufferSize : int) {
  var localBuffers : [0 ..# numLocales, 0 ..# numTasks] _LocalBuffer(complex(128))
    = [(srcLocaleIdx, taskIdx) in {0 ..# numLocales, 0 ..# numTasks}]
        new _LocalBuffer(complex(128), here.id, srcLocaleIdx, remoteBufferSize);
  return localBuffers;
}

proc _offDiagMakeRemoteBuffers(numTasks : int, remoteBufferSize : int) {
  var remoteBuffers : [0 ..# numLocales, 0 ..# numTasks] _RemoteBuffer(complex(128))
    = [(destLocaleIdx, taskIdx) in {0 ..# numLocales, 0 ..# numTasks}]
        new _RemoteBuffer(complex(128), destLocaleIdx, here.id, remoteBufferSize);
  return remoteBuffers;
}

proc _offDiagInitLocalBuffers(numTasks : int, ref localBuffers, const ref ptrStore) {
  const destLocaleIdx = localBuffers.locale.id;
  coforall loc in Locales do on loc {
    const srcLocaleIdx = loc.id;
    const (_basis, _accessor, _remoteBufferPtr, _localBufferPtr, _numChunks) = ptrStore[srcLocaleIdx];
    for taskIdx in 0 ..# numTasks {
      ref remoteBuffer = (_remoteBufferPtr + destLocaleIdx * numTasks + taskIdx).deref();
      localBuffers[srcLocaleIdx, taskIdx].isFull = c_ptrTo(remoteBuffer.isFull);
    }
  }
}

proc _offDiagInitRemoteBuffers(numTasks : int, ref remoteBuffers, const ref ptrStore) {
  const srcLocaleIdx = remoteBuffers.locale.id;
  coforall loc in Locales do on loc {
    const destLocaleIdx = loc.id;
    const (_basis, _accessor, _remoteBufferPtr, _localBufferPtr, _numChunks) = ptrStore[destLocaleIdx];
    for taskIdx in 0 ..# numTasks {
      ref myLocalBuffer = (_localBufferPtr + srcLocaleIdx * numTasks + taskIdx).deref();
      ref remoteBuffer = remoteBuffers[destLocaleIdx, taskIdx];
      remoteBuffer.basisStates = myLocalBuffer.basisStates;
      remoteBuffer.coeffs = myLocalBuffer.coeffs;
      remoteBuffer.size = c_ptrTo(myLocalBuffer.size);
      remoteBuffer.isEmpty = c_ptrTo(myLocalBuffer.isEmpty);
      remoteBuffer.isEOF = c_ptrTo(myLocalBuffer.isEOF);
    }
  }
}

record Producer {
  type eltType;
  var _taskIdx : int;
  var numChunks : int;
  var numProducerTasks : int;

  var batchedOperator : BatchedOperator;
  var basisPtr : c_ptr(Basis);
  var accessorPtr : c_ptr(ConcurrentAccessor(eltType));
  var representativesPtr : c_ptr(uint(64));
  var xPtr : c_ptr(eltType);

  var rangesPtr : c_ptr(range(int));
  var moreWorkPtr : c_ptr(atomic bool);
  var currentChunkIdxPtr : c_ptr(atomic int);

  var remoteBuffersPtr : c_ptr(_RemoteBuffer(complex(128)));

  proc init(taskIdx : int, numChunks : int, in batchedOperator : BatchedOperator,
            ref accessor : ConcurrentAccessor(?eltType),
            const ref representatives : [] uint(64),
            const ref x : [] eltType,
            const ref ranges : [] range(int),
            ref moreWork : atomic bool,
            ref currentChunkIdx : atomic int,
            ref remoteBuffers : [] _RemoteBuffer(complex(128))) {
    // logDebug("Creating Producer(", taskIdx, ")...");
    this.eltType = eltType;
    this._taskIdx = taskIdx;
    this.numChunks = numChunks;
    this.numProducerTasks = remoteBuffers.domain.dim(1).size;
    this.batchedOperator = batchedOperator;
    this.basisPtr = c_const_ptrTo(this.batchedOperator._matrixPtr.deref().basis);
    this.accessorPtr = c_ptrTo(accessor);
    this.representativesPtr = c_const_ptrTo(representatives);
    this.xPtr = c_const_ptrTo(x);
    this.rangesPtr = c_const_ptrTo(ranges);
    this.moreWorkPtr = c_ptrTo(moreWork);
    this.currentChunkIdxPtr = c_ptrTo(currentChunkIdx);
    this.remoteBuffersPtr = c_ptrTo(remoteBuffers);
    // logDebug("Done creating Producer(", taskIdx, ")...");
  }

  proc remoteBuffers(localeIdx : int, taskIdx : int) ref {
    assert(0 <= localeIdx && localeIdx < numLocales);
    assert(0 <= taskIdx && taskIdx < numProducerTasks);
    return remoteBuffersPtr[localeIdx * numProducerTasks + taskIdx];
  }

  proc run() {
    // logDebug("Running Producer(", taskIdx, ")...");
    while moreWorkPtr.deref().read() {
      const rangeIdx = currentChunkIdxPtr.deref().fetchAdd(1);
      // Multiple threads passed moreWork.read() at once.
      // All whose fetchAdd() was after the one
      // that grabbed the final chunk just break.
      if rangeIdx >= numChunks then
        break;
      // Final rangeIdx -- signal that to everybody
      if rangeIdx == numChunks - 1 then
        moreWorkPtr.deref().write(false);

      // logDebug("Doing work in Producer(", taskIdx, ")...");
      const r : range(int) = rangesPtr[rangeIdx];
      // Compute a r.size rows of the matrix
      var timer = new Timer();
      timer.start();
      // logDebug("Calling computeOffDiag in Producer(", taskIdx, ")...");
      const (n, basisStatesPtr, coeffsPtr, keysPtr) = batchedOperator.computeOffDiag(
          r.size, representativesPtr + r.low, xPtr + r.low);
      var radixOffsets : c_array(int, 257);
      radixOneStep(n, keysPtr, radixOffsets, basisStatesPtr, coeffsPtr);
      timer.stop();
      // timers.computeOffDiag.add(timer.elapsed(), memoryOrder.relaxed);

      for destLocaleIdx in 0 ..# numLocales {
        const k = radixOffsets[destLocaleIdx];
        const n = radixOffsets[destLocaleIdx + 1] - k;
        ref remoteBuffer = remoteBuffers[destLocaleIdx, _taskIdx];
        if destLocaleIdx == here.id {
          assert(!remoteBuffer.isFull.read());
          // logDebug("Calling localProcess in Producer(", taskIdx, ")...");
          localProcess(basisPtr,
                       accessorPtr,
                       basisStatesPtr + k,
                       coeffsPtr + k,
                       n);
        }
        else {
          // logDebug("Calling submit in Producer(", _taskIdx, ")...");
          remoteBuffer.submit(basisStatesPtr + k, coeffsPtr + k, n);
        }
      }
      // logDebug("Done doing work in Producer(", taskIdx, ")");
    }
    // logDebug("Done running Producer(", taskIdx, ")");
  }
}

record Consumer {
  type eltType;

  var _taskIdx : int;
  var numConsumerTasks : int;
  var numProducerTasks : int;

  var slots;

  var basisPtr : c_ptr(Basis);
  var accessorPtr : c_ptr(ConcurrentAccessor(eltType));

  var totalNumberChunks : int;
  var numProcessedPtr : c_ptr(atomic int);

  var localBuffersPtr : c_ptr(_LocalBuffer(complex(128)));

  proc init(taskIdx : int, numConsumerTasks : int, numProducerTasks : int,
            const ref basis : Basis, ref accessor : ConcurrentAccessor(?eltType),
            totalNumberChunks : int, ref numProcessedChunks : atomic int,
            ref localBuffers : [] _LocalBuffer(complex(128))) {
    // logDebug("Creating Consumer(", taskIdx, ")...");
    this.eltType = eltType;
    this._taskIdx = taskIdx;
    this.numConsumerTasks = numConsumerTasks;
    this.numProducerTasks = numProducerTasks;
   
    // logDebug(0 ..# (numLocales - 1) * numProducerTasks);
    const everything = 0 ..# (numLocales - 1) * numProducerTasks;
    if taskIdx < everything.size {
      const r = chunk(everything, numConsumerTasks, taskIdx);
      var pairs : [0 ..# r.size] (int, int);
      var size = 0;
      var offset = 0;
      for localeIdx in 0 ..# numLocales {
        if localeIdx == here.id then continue;
        for otherTaskIdx in 0 ..# numProducerTasks {
          if r.first <= offset && offset <= r.last {
            pairs[size] = (localeIdx, otherTaskIdx);
            size += 1;
          }
          offset += 1;
        }
      }
      this.slots = pairs;
      logDebug("Slots: ", slots);
    }
    else {
      var pairs : [0 ..# 0] (int, int);
      this.slots = pairs;
    }

    this.basisPtr = c_const_ptrTo(basis);
    this.accessorPtr = c_ptrTo(accessor);
    this.totalNumberChunks = totalNumberChunks;
    this.numProcessedPtr = c_ptrTo(numProcessedChunks);
    this.localBuffersPtr = c_ptrTo(localBuffers[0, 0]);
  }

  proc localBuffers(localeIdx : int, taskIdx : int) ref {
    assert(0 <= localeIdx && localeIdx < numLocales);
    assert(0 <= taskIdx && taskIdx < numProducerTasks);
    return localBuffersPtr[localeIdx * numProducerTasks + taskIdx];
  }

  proc run() {
    // logDebug("numProcessed = ", numProcessedPtr.deref().read(),
    //          ", totalNumberChunks = ", totalNumberChunks);
    while numProcessedPtr.deref().read() < totalNumberChunks {
      for (localeIdx, otherTaskIdx) in slots {
      // for localeIdx in 0 ..# numLocales {
      //   for otherTaskIdx in 0 ..# numProducerTasks {
          ref localBuffer = localBuffers[localeIdx, otherTaskIdx];
          if !localBuffer.isEmpty.read() {
            // logDebug("Calling localProcess in Consumer(", _taskIdx, ")...");
            local {
              numProcessedPtr.deref().add(1);

              localProcess(basisPtr, accessorPtr,
                           localBuffer.basisStates,
                           localBuffer.coeffs,
                           localBuffer.size);
              localBuffer.isEmpty.write(true);
            }
            const atomicPtr = localBuffer.isFull;


            // logDebug("Calling on Locales[...] in Consumer(", _taskIdx, ")...");
            on Locales[localBuffer.srcLocaleIdx] {
              atomicStoreBool(atomicPtr, false);
            }

            // logDebug("Done with localProcess in Consumer(", _taskIdx, ")...");
          }
      //  }
      }
      chpl_task_yield();
    }
  }
}

config const kUseConsumer : bool = false;

private proc localOffDiagonalNoQueue(matrix : Operator, const ref x : [] ?eltType, ref y : [] eltType,
                                     const ref representatives : [] uint(64),
                                     numChunks : int = min(matrixVectorOffDiagonalNumChunks,
                                                           representatives.size)) {
  var timers = new LocalOffDiagonalTimers();
  timers.total.start();
  timers.initialization.start();

  const numTasks = kNumTasks;
  const numConsumerTasks = min(kNumConsumerTasks, numTasks - 1);
  const numProducerTasks = numTasks - numConsumerTasks;

  const remoteBufferSize = max(kRemoteBufferSize, matrix.numberOffDiagTerms());
  // Other locales need access to our newLocalBuffers, newRemoteBuffers, etc.
  // so we create and store them in a common location (on locale 0) before a
  // barrier.
  var newLocalBuffers = _offDiagMakeLocalBuffers(numProducerTasks, remoteBufferSize);
  var newRemoteBuffers = _offDiagMakeRemoteBuffers(numProducerTasks, remoteBufferSize);
  var accessor = new ConcurrentAccessor(y);
  const numChunks = max(representatives.size * matrix.numberOffDiagTerms() / remoteBufferSize, 1);
  globalPtrStoreNoQueue[here.id] = (c_const_ptrTo(matrix.basis),
                                    c_ptrTo(accessor),
                                    c_ptrTo(newRemoteBuffers[0, 0]),
                                    c_ptrTo(newLocalBuffers[0, 0]),
                                    numChunks);
  allLocalesBarrier.barrier();

  const ptrStore : [0 ..# numLocales] globalPtrStoreNoQueue.eltType = globalPtrStoreNoQueue;
  _offDiagInitLocalBuffers(numProducerTasks, newLocalBuffers, ptrStore);
  _offDiagInitRemoteBuffers(numProducerTasks, newRemoteBuffers, ptrStore);
  allLocalesBarrier.barrier();
  timers.initialization.stop();

  const ranges : [0 ..# numChunks] range(int) = chunks(0 ..# representatives.size, numChunks);
  const batchedOperatorChunkSize = (representatives.size + numChunks - 1) / numChunks;
  var moreWork : atomic bool = true;
  var curChunkIdx : atomic int = 0;

  var totalNumberChunks = 0;
  for localeIdx in 0 ..# numLocales {
    if localeIdx != here.id {
      const (_basis, _accessor, _remoteBufferPtr, _localBufferPtr, _numChunks) = ptrStore[localeIdx];
      totalNumberChunks += _numChunks;
    }
  }
  
  var numProcessedChunks : atomic int = 0;

  logDebug("numChunks = ", numChunks, ", chunkSize = ", batchedOperatorChunkSize, ", numTasks = ", numTasks);

  // allLocalesBarrier.barrier();

  coforall taskIdx in 0 ..# numTasks with (ref timers, ref accessor) {
    if taskIdx < numProducerTasks {
      // I'm a producer
      var producer = new Producer(
        taskIdx,
        numChunks,
        new BatchedOperator(matrix, batchedOperatorChunkSize),
        accessor,
        representatives,
        x,
        ranges,
        moreWork,
        curChunkIdx,
        newRemoteBuffers);
      // logDebug("ptr(matrix.basis) = ", c_const_ptrTo(matrix.basis),
      //          ", producer.basisPtr = ", producer.basisPtr);
      producer.run();
    }
    else {

      if kUseConsumer {
        var consumer = new Consumer(
          taskIdx - numProducerTasks,
          numConsumerTasks,
          numProducerTasks,
          matrix.basis,
          accessor,
          totalNumberChunks,
          numProcessedChunks,
          newLocalBuffers);

        consumer.run();
      }
      else {
        // startVerboseCommHere();
        // I'm a consumer
        // logDebug("totalNumberChunks = ", totalNumberChunks);
        var consumer = new Consumer(
          taskIdx - numProducerTasks,
          numConsumerTasks,
          numProducerTasks,
          matrix.basis,
          accessor,
          totalNumberChunks,
          numProcessedChunks,
          newLocalBuffers);

        // logDebug("Starting Consumer(", taskIdx, ")...");
        while consumer.numProcessedPtr.deref().read() < totalNumberChunks {
          for localeIdx in 0 ..# numLocales {
            for otherTaskIdx in 0 ..# consumer.numProducerTasks {
              ref localBuffer = newLocalBuffers[localeIdx, otherTaskIdx];
              if !localBuffer.isEmpty.read() {
                // logDebug("Calling localProcess in Consumer(", taskIdx, ")...");
                localProcess(consumer.basisPtr,
                             consumer.accessorPtr,
                             // c_const_ptrTo(matrix.basis), c_const_ptrTo(accessor),
                             localBuffer.basisStates,
                             localBuffer.coeffs,
                             localBuffer.size);
                localBuffer.isEmpty.write(true);
                const atomicPtr = localBuffer.isFull;

                // logDebug("Calling on Locales[...] in Consumer(", taskIdx, ")...");
                on Locales[localBuffer.srcLocaleIdx] {
                  atomicStoreBool(atomicPtr, false);
                }

                consumer.numProcessedPtr.deref().add(1);
                // logDebug("Done with localProcess in Consumer(", taskIdx, ")...");
              }
            }
          }
          chpl_task_yield();
        }
        // logDebug("Done with Consumer(", taskIdx, ")...");
        // stopVerboseCommHere();

      }

      
    }

    /*
    while moreWork.read() {
      record InTheMeantime {
        inline proc this() {
          return tryProcessLocal(taskIdx, newLocalBuffers, matrix, accessor);
        }
      }
      const inTheMeantime = new InTheMeantime();

      const rangeIdx = curChunkIdx.fetchAdd(1);
      // Multiple threads passed moreWork.read() at once.
      // All whose fetchAdd() was after the one
      // that grabbed the final chunk just break.
      if rangeIdx >= numChunks then
        break;
      // Final rangeIdx -- signal that to everybody
      if rangeIdx == numChunks - 1 then
        moreWork.write(false);
      const r : range(int) = ranges[rangeIdx];

      // Compute a r.size rows of the matrix
      var timer = new Timer();
      timer.start();
      const (n, basisStatesPtr, coeffsPtr, keysPtr) = batchedOperator.computeOffDiag(
          r.size, c_const_ptrTo(representatives[r.low]), c_const_ptrTo(x[r.low]));
      var radixOffsets : c_array(int, 257);
      radixOneStep(n, keysPtr, radixOffsets, basisStatesPtr, coeffsPtr);
      timer.stop();
      timers.computeOffDiag.add(timer.elapsed(), memoryOrder.relaxed);

      for destLocaleIdx in 0 ..# numLocales {
        const k = radixOffsets[destLocaleIdx];
        const n = radixOffsets[destLocaleIdx + 1] - k;
        ref remoteBuffer = newRemoteBuffers[destLocaleIdx, taskIdx];
        // if n == 0 then continue;
        if destLocaleIdx == here.id {
          assert(!remoteBuffer.isFull.read());
          localProcess(c_const_ptrTo(matrix.basis), c_const_ptrTo(accessor),
                       basisStatesPtr + k,
                       coeffsPtr + k,
                       n);
        }
        else {
          remoteBuffer.submit(basisStatesPtr + k, coeffsPtr + k, n, inTheMeantime);
          inTheMeantime();
        }
      }
    }
    */

    /*
    coforall localeIdx in 0 ..# numLocales {
      ref remoteBuffer = newRemoteBuffers[localeIdx, taskIdx];
      while remoteBuffer.isFull.read() {
        if !tryProcessLocal(taskIdx, localeIdx, newLocalBuffers, matrix, accessor) then
          chpl_task_yield();
      }
      remoteBuffer.finish();

      ref localBuffer = newLocalBuffers[localeIdx, taskIdx];
      while !localBuffer.isEOF.read() {
        if !tryProcessLocal(taskIdx, localeIdx, newLocalBuffers, matrix, accessor) then
          chpl_task_yield();
      }
    }
    */
  }

  allLocalesBarrier.barrier();

  for srcLocaleIdx in 0 ..# numLocales {
    for taskIdx in 0 ..# numProducerTasks {
      // assert(newLocalBuffers[srcLocaleIdx, taskIdx].isEOF.read());
      assert(newLocalBuffers[srcLocaleIdx, taskIdx].isEmpty.read());
    }
  }
  for srcLocaleIdx in 0 ..# numLocales {
    for taskIdx in 0 ..# numProducerTasks {
      assert(!newRemoteBuffers[srcLocaleIdx, taskIdx].isFull.read());
    }
  }

  timers.put = + reduce newRemoteBuffers._putTime;
  timers.wait = + reduce newRemoteBuffers._waitTime;
  timers.submit = + reduce newRemoteBuffers._submitTime;

  timers.total.stop();
  logDebug("localOffDiagonalNoQueue took ", timers.total.elapsed(), "\n",
           "  ├─ ", timers.initialization.elapsed(), " in initialization\n",
           "  ├─ ", timers.computeOffDiag, " in computeOffDiag\n",
           "  ├─ ", timers.submit, " in remoteBuffer.submit\n",
           "      ├─ ", timers.wait, " in remoteBuffer.isFull.waitFor(false)\n",
           "      └─ ", timers.put, " in remote PUTs\n");
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
