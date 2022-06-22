module DistributedMatrixVector {

use CTypes;
use RangeChunk;
use AllLocalesBarriers;

use FFI;
use ForeignTypes;
use ConcurrentAccessor;
use BatchedOperator;
use CommunicationQueue;

/* 
 */
private proc localDiagonalBatch(indices : range(int, BoundedRangeType.bounded, false),
                                ref workspace : [] complex(128), matrix : Operator,
                                const ref x : [] ?eltType, ref y : [] eltType,
                                const ref representatives : [] uint(64)) {
  const batchSize = indices.size;
  assert(workspace.size >= batchSize);
  // if workspace.size < batchSize then
  //   workspace.domain = {0 ..# batchSize};
  ls_hs_operator_apply_diag_kernel(
    matrix.payload, batchSize,
    c_const_ptrTo(representatives[indices.low]), 1,
    c_ptrTo(workspace));
  foreach i in indices {
    y[i] = x[i] * workspace[i - indices.low]:eltType;
  }
}

config const matrixVectorDiagonalNumChunks : int = here.maxTaskPar;
config const matrixVectorOffDiagonalNumChunks : int = 10 * here.maxTaskPar;

private proc localDiagonal(matrix : Operator, const ref x : [] ?eltType, ref y : [] eltType,
                           const ref representatives : [] uint(64),
                           numChunks : int = min(matrixVectorDiagonalNumChunks,
                                                 representatives.size)) {
  const totalSize = representatives.size;
  const batchSize = (totalSize + numChunks - 1) / numChunks;
  var ranges : [0 ..# numChunks] range(int, BoundedRangeType.bounded, false) =
    chunks(0 ..# totalSize, numChunks);
  var workspace : [0 ..# batchSize] complex(128) = noinit;
  forall r in ranges with (in workspace) {
    localDiagonalBatch(r, workspace, matrix, x, y, representatives);
  }
}

private proc localOffDiagonal(matrix : Operator, const ref x : [] ?eltType, ref y : [] eltType,
                              const ref representatives : [] uint(64),
                              numChunks : int = min(matrixVectorOffDiagonalNumChunks,
                                                    representatives.size)) {
  const chunkSize = (representatives.size + numChunks - 1) / numChunks;
  logDebug("Using chunkSize=", chunkSize);
  // Used to calculate the action of matrix on a chunk of basis elements.
  // batchedOperator stores some task-local buffers and a pointer to `matrix`.
  var batchedOperator = new BatchedOperator(matrix, chunkSize);
  // Used to then process the matrix elements.
  // logDebug("new ConcurrentAccessor: " + y.locale:string);
  var accessor = new ConcurrentAccessor(y);
  globalAllQueues[here.id] = new CommunicationQueue(matrix.basis, accessor);
  ref queue = globalAllQueues[here.id];
  var staging = new StagingBuffers();
  allLocalesBarrier.barrier();
  // writeln(globalAllQueues);

  // logDebug(representatives.size:string + " vs. " + numChunks:string);
  var ranges : [0 ..# numChunks] range(int) =
    chunks(0 ..# representatives.size, numChunks);
  forall r in ranges
      with (in batchedOperator, in staging) {
    logDebug("Processing ", r, " ...");
    const (basisStatesPtr, coeffsPtr, offsetsPtr) = batchedOperator.computeOffDiag(
        r.size, c_const_ptrTo(representatives[r.low]), c_const_ptrTo(x[r.low]));
    const n = offsetsPtr[r.size];
    staging.add(n, basisStatesPtr, coeffsPtr);
    staging.flush();
  }
  queue!.drain();
  allLocalesBarrier.barrier();
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
    localMatrixVector(myMatrix, myX, myY, myBasisStates);
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
