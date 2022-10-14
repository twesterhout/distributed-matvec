module Diagonalize {
  use AllLocalesBarriers;
  use CTypes;
  use CommDiagnostics;
  import Random;

  use PRIMME;
  use LatticeSymmetries;

  /*
  // A buffer located on locale 0 to help with the broadcast                      
  var tmpBuffDom = {0..0:c_int};
  var tmpBuff: [tmpBuffDom] real;

  proc broadcastReal(buffer: c_ptr(real), count: c_ptr(c_int)) {
    const n = count.deref(),
          inds = 0..<n;

    if here.id == 0 {
      // grow the temp buff if it's not big enough                                
      if n > tmpBuffDom.size then
        tmpBuffDom = {inds};

      // copy locale 0's data into the buffer                                     
      forall i in inds do
        tmpBuff[i] = buffer[i];
    }

    // wait until locale 0's got tmpBuff set up before proceeding                 
    allLocalesBarrier.barrier();

    // Locale 0 already has the data so doesn't need to do anything               
    if (here.id != 0) then
      forall i in inds do
        buffer[i] = tmpBuff[i];
  }

  export proc ls_chpl_broadcast_real(sendBuf : c_void_ptr, count : c_ptr(c_int),
                                     primme : c_ptr(primme_params), ierr : c_ptr(c_int)) {
    if primme.deref().broadcastReal_type != primme_op_double then
      halt("broadcastReal is implemented for double precision only");
    broadcastReal(sendBuf : c_ptr(real), count);
    ierr.deref() = 0;
  }


  record AtomicBuf {
    var dom : domain(1);
    var arr : [dom] atomic real;
  }

  // A buffer of atomics on locale 0 for computing the reduction                  
  // var atomicBuffDom: domain(1); // = {0..0:c_int};
  var atomicBuff = new AtomicBuf();
  // : [atomicBuffDom] atomic real;

  // var D : domain(1);
  // var A : [D] real;

  proc globalSumReal(sendBuf: c_ptr(real), recvBuf: c_ptr(real),
                     count: c_ptr(c_int)) {
    // startVerboseCommHere();
    const n = count.deref();
    const inds = 0 ..# n:int;
    // writeln("calling globalSumReal(", sendBuf, ", ", recvBuf, ", ", n, ") from ", here);

    const sendArr = [i in inds] sendBuf[i];
    writeln(here, ": sendBuf=", ([i in inds] sendBuf[i]):string);
    writeln(here, ": sendArr=", sendArr);
    writeln(here, ": recvBuf=", ([i in inds] recvBuf[i]):string);

    // ref atomicBuffRef = atomicBuff;

    // grow the temp buff if it's not big enough                                  
    if here.id == 0 then
      if n > atomicBuff.dom.size {
        // atomicBuffDom = {inds};
        // writeln(atomicBuffDom);
        // writeln(atomicBuff);
        atomicBuff.dom = {inds};

        // D = {0 ..# n:int};
      }

    // Make sure locale 0 has had the chance to resize before proceeding          
    allLocalesBarrier.barrier();
    assert(atomicBuff.dom.size >= n);
    assert(atomicBuff.arr.size >= n);
    // writeln(here, ": D = ", D, ", A = ", A);

    // allLocalesBarrier.barrier();
    // assert(D.size >= n);
    // assert(A.domain.size >= n);
    // assert(atomicBuffRef.size >= n);

    // have all locales atomically add their results to the atomicBuff            
    for i in inds {
      writeln("adding " + sendBuf[i]:string + " to atomicBuff.arr[" + i:string + "]");
      const elt = sendBuf[i];
      assert(elt == sendArr[i]);
      atomicBuff.arr[i].add(sendBuf[i]);
    }

    // Make sure all locales have accumulated their contributions                 
    allLocalesBarrier.barrier();

    // Have each locale copy the results out into its buffer                      
    for i in inds do
      recvBuf[i] = atomicBuff.arr[i].read();

    writeln(here, ": recvBuf=", ([i in inds] recvBuf[i]):string);
    // stopVerboseCommHere();
  }

  export proc ls_chpl_global_sum_real(sendBuf : c_void_ptr, recvBuf : c_void_ptr,
                                      count : c_ptr(c_int), primme : c_ptr(primme_params),
                                      ierr : c_ptr(c_int)) {
    if primme.deref().globalSumReal_type != primme_op_double then
      halt("globalSum is implemented for double precision only");
    globalSumReal(sendBuf : c_ptr(real), recvBuf : c_ptr(real), count);
    ierr.deref() = 0;
  }

  */

  inline proc borrowOperator(primme : c_ptr(primme_params)) {
    const p = primme.deref().matrix : c_ptr(ls_hs_operator);
    return new Operator(p, owning=false);
  }

  export proc ls_chpl_primme_matvec(_x : c_void_ptr, _ldx : c_ptr(int(64)),
                                    _y : c_void_ptr, _ldy : c_ptr(int(64)),
                                    _blockSize : c_ptr(c_int), primme : c_ptr(primme_params),
                                    _ierr : c_ptr(c_int)) {
    logDebug("Calling ls_chpl_primme_matvec ...");
    const blockSize = _blockSize.deref():int;
    // if blockSize != 1 then
    //   halt("currently only blockSize=1 is supported");
    const ldx = _ldx.deref();
    const ldy = _ldy.deref();
    const n = primme.deref().nLocal;
    assert(ldx >= n);
    assert(ldy >= n);

    // const precision = primme.deref().matrixMatvec_type;
    // const dtype = 
    type eltType = real;
    const matrix = borrowOperator(primme);

    for k in 0 ..# blockSize {
      ref X = makeArrayFromPtr(_x : c_ptr(eltType) + ldx * k, (n,));
      ref Y = makeArrayFromPtr(_y : c_ptr(eltType) + ldy * k, (n,));
      localMatrixVector(matrix, X, Y, matrix.basis.representatives());
    }

    _ierr.deref() = 0;
    logDebug("Done with ls_chpl_primme_matvec ...");
  }

  proc main() {
    initRuntime();
    defer deinitRuntime();
    writeln("Hello world!");

    const filename = "data/heisenberg_chain_10.yaml";
    const matrix = loadHamiltonianFromYaml(filename);

    const masks;
    const basisStates = enumerateStates(matrix.basis, masks);

    const numEvals = 2;
    var evecs = new BlockVector(real(64), numEvals, basisStates.counts);

    // var matrix : [0 ..# 4, 0 ..# 4] real;
    // Random.fillRandom(matrix, seed=53);

    // Symmetrize
    // for i in matrix.dim(0) do
    //   for j in matrix.dim(1) do
    //     if i > j then
    //       matrix[i, j] = matrix[j, i];

    // writeln(matrix);


    coforall loc in Locales with (ref evecs)
                            do on loc {

      // var myEvecs : [0 ..# (4 / numLocales)] real;
      // var myEvals : [0 ..# 1] real;
      // var myResNorms : [0 ..# 1] real;
      // const myMatrix = matrix[loc.id * (4 / numLocales) .. (loc.id + 1) * (4 / numLocales), ..];
      // logDebug("initially: ", globalPtrStoreNoQueue.arr[here.id]);


      var myMatrix = loadHamiltonianFromYaml(filename);
      const ref myBasisStates = basisStates.getBlock(loc.id);
      myMatrix.basis.uncheckedSetRepresentatives(myBasisStates);

      ref myEvecs = evecs.getBlock(loc.id);
      var myEvals : [0 ..# numEvals] real(64);
      var myResNorms : [0 ..# numEvals] real(64);

      var params : primme_params;
      primme_initialize(params);
      // params.matrix = c_const_ptrTo(myMatrix[0, 0]);
      // params.matrixMatvec = c_ptrTo(primmeMatrixMatvec);

      params.matrix = myMatrix.payload;
      params.matrixMatvec = c_ptrTo(ls_chpl_primme_matvec);
      params.maxBlockSize = 1;
      params.n = + reduce basisStates.counts;
      params.numEvals = numEvals:c_int;
      params.eps = 1e-6;
      params.target = primme_smallest;

      params.numProcs = numLocales:int(32);
      params.procID = loc.id:int(32);
      params.nLocal = myBasisStates.size; // myEvecs.dim(0).size;
      params.globalSumReal = c_ptrTo(primmeGlobalSumReal);
      params.globalSumReal_type = primme_op_double;
      params.broadcastReal = c_ptrTo(primmeBroadcastReal);
      params.broadcastReal_type = primme_op_double;

      params.printLevel = 5;

      primme_set_method(PRIMME_DEFAULT_METHOD, params);
      primme_display_params(params);

      allLocalesBarrier.barrier();

      writeln(here, ": calling dprimme ...");
      allLocalesBarrier.barrier();

      const ierr = dprimme(c_ptrTo(myEvals[0]),
                           c_ptrTo(myEvecs[myEvecs.domain.low]),
                           c_ptrTo(myResNorms[0]),
                           params);

      writeln(ierr, ": ", myEvals, ", ", myResNorms);

      primme_free(params);
    }
  }
}
