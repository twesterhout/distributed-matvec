module MyHDF5 {
  use CTypes;
  use BlockDist;
  import HDF5;
  import HDF5.C_HDF5;

  use LatticeSymmetries.FFI;

  /*
  proc datasetRank(filename : string, dataset : string) {
    return ls_hs_hdf5_get_dataset_rank(
      filename.localize().c_str(), dataset.localize().c_str()):int;
  }
  */

  /* Get shape of a dataset in a HDF5 file.
   *
   * :arg filename: Path to HDF5 file.
   * :arg dataset:  Path to dataset within HDF5 file.
   * 
   * :returns: shape of the dataset as a one-dimensional array.
   * :rtype: [] int
   */
  proc datasetShape(filename : string, dataset : string, param rank : int) : rank * int {
    const file_id = C_HDF5.H5Fopen(filename.localize().c_str(),
                                   C_HDF5.H5F_ACC_RDONLY, C_HDF5.H5P_DEFAULT);
    defer C_HDF5.H5Fclose(file_id);

    const dset_id = C_HDF5.H5Dopen(file_id, dataset.localize().c_str(), C_HDF5.H5P_DEFAULT);
    defer C_HDF5.H5Dclose(dset_id);

    const dspace_id = C_HDF5.H5Dget_space(dset_id);
    defer C_HDF5.H5Sclose(dspace_id);

    const dset_rank = C_HDF5.H5Sget_simple_extent_ndims(dspace_id);
    if dset_rank != rank then
      halt("expected '" + dataset + "' to be " + rank:string +
           "-dimensional, but it has rank " + dset_rank:string);
    var c_shape : [0 ..# rank] C_HDF5.hsize_t;
    C_HDF5.H5Sget_simple_extent_dims(dspace_id, c_ptrTo(c_shape), nil);

    return _asTuple(c_shape, rank):(rank * int);
  }

  inline proc _domainFromTuple((a,) : 1 * range) : domain(1) { return { a }; }
  inline proc _domainFromTuple((a, b) : 2 * range) : domain(2) { return { a, b }; }
  inline proc _domainFromTuple((a, b, c) : 3 * range) : domain(3) { return { a, b, c }; }
  inline proc _domainFromTuple((a, b, c, d) : 4 * range) : domain(4) { return { a, b, c, d }; }

  private inline proc _rangesFromShape(shape) {
    param rank = shape.size;
    var ranges : rank * range;
    for (r, n) in zip(ranges, shape) do
      r = 0 ..# n;
    return ranges;
  }

  private inline proc _makeDomain(shape) where isTuple(shape) {
    return _domainFromTuple(_rangesFromShape(shape));
  }

  private inline proc _asTuple(xs : [] ?eltType, param rank : int) {
    var ys : rank * eltType;
    for (y, x) in zip(ys, xs) do
      y = x;
    return ys;
  }

  proc readDataset(filename : string, dataset : string, type eltType, param rank : int) {
    const file_id = C_HDF5.H5Fopen(filename.localize().c_str(),
                                   C_HDF5.H5F_ACC_RDONLY, C_HDF5.H5P_DEFAULT);
    defer C_HDF5.H5Fclose(file_id);

    const dset_id = C_HDF5.H5Dopen(file_id, dataset.localize().c_str(), C_HDF5.H5P_DEFAULT);
    defer C_HDF5.H5Dclose(dset_id);

    const dspace_id = C_HDF5.H5Dget_space(dset_id);
    defer C_HDF5.H5Sclose(dspace_id);

    const datasetRank = C_HDF5.H5Sget_simple_extent_ndims(dspace_id);
    if rank != datasetRank then
      halt("rank mismatch in file: '" + filename + "' dataset: '" + dataset +
           "'  " + rank:string + " != " + datasetRank:string);

    const dtype_id = C_HDF5.H5Dget_type(dset_id);
    defer C_HDF5.H5Tclose(dtype_id);
    if C_HDF5.H5Tequal(HDF5.getHDF5Type(eltType), dtype_id) <= 0 then
      halt("type mismatch in file: '" + filename + "' dataset: '" + dataset +
           "'  " + HDF5.getHDF5Type(eltType):string + " != " + dtype_id:string);
  
    var c_shape : [0 ..# rank] C_HDF5.hsize_t;
    C_HDF5.H5Sget_simple_extent_dims(dspace_id, c_ptrTo(c_shape), nil);

    const mspace_id = C_HDF5.H5Screate_simple(rank, c_ptrTo(c_shape), nil);
    defer C_HDF5.H5Sclose(mspace_id);

    var arr : [_makeDomain(c_shape:int)] eltType;
    C_HDF5.H5Dread(dset_id, dtype_id, mspace_id, dspace_id,
                   C_HDF5.H5P_DEFAULT, c_ptrTo(arr[arr.domain.low]));
    return arr;
  }

  proc readDatasetChunk(filename : string, dataset : string, offset, ref arr : [] ?eltType)
      where isTuple(offset) && offset.size == arr.rank {
    const file_id = C_HDF5.H5Fopen(filename.localize().c_str(),
                                   C_HDF5.H5F_ACC_RDONLY, C_HDF5.H5P_DEFAULT);
    defer C_HDF5.H5Fclose(file_id);

    const dset_id = C_HDF5.H5Dopen(file_id, dataset.localize().c_str(), C_HDF5.H5P_DEFAULT);
    defer C_HDF5.H5Dclose(dset_id);

    const dspace_id = C_HDF5.H5Dget_space(dset_id);
    defer C_HDF5.H5Sclose(dspace_id);

    const datasetRank = C_HDF5.H5Sget_simple_extent_ndims(dspace_id);
    if arr.rank != datasetRank then
      halt("rank mismatch in file: '" + filename + "' dataset: '" + dataset +
           "'  " + arr.rank:string + " != " + datasetRank:string);

    const dtype_id = C_HDF5.H5Dget_type(dset_id);
    defer C_HDF5.H5Tclose(dtype_id);
    if C_HDF5.H5Tequal(HDF5.getHDF5Type(eltType), dtype_id) <= 0 then
      halt("type mismatch in file: '" + filename + "' dataset: '" + dataset +
           "'  " + HDF5.getHDF5Type(eltType):string + " != " + dtype_id:string);
    
    var c_offset : [0 ..# arr.rank] C_HDF5.hsize_t;
    var c_shape : [0 ..# arr.rank] C_HDF5.hsize_t;
    for i in 0 ..# arr.rank {
      c_offset[i] = offset[i]:C_HDF5.hsize_t;
      c_shape[i] = arr.dim(i).size:C_HDF5.hsize_t;
    }

    const mspace_id = C_HDF5.H5Screate_simple(arr.rank, c_ptrTo(c_shape), nil);
    defer C_HDF5.H5Sclose(mspace_id);

    C_HDF5.H5Sselect_hyperslab(dspace_id, C_HDF5.H5S_SELECT_SET,
                               c_ptrTo(c_offset), nil,
                               c_ptrTo(c_shape), nil);

    C_HDF5.H5Dread(dset_id, dtype_id, mspace_id, dspace_id,
                   C_HDF5.H5P_DEFAULT, c_ptrTo(arr[arr.domain.low]));
  }

  /* Read part of a dataset from a HDF5 file.
   *
   * :arg filename: Path to HDF5 file.
   * :arg dataset:  Path to dataset within HDF5 file.
   * :arg eltType:  Array datatype.
   * :arg offset:   A tuple of offsets along each dimension.
   * :arg shape:    Array shape.
   *
   * :returns: part of the dataset which is read from file.
   * :rtype: [] eltType
   */
  proc readDatasetChunk(filename : string, dataset : string, type eltType, offset, shape) {
    const dom = _makeDomain(shape);
    var array : [dom] eltType;
    readDatasetChunk(filename, dataset, offset, array);
    return array;
  }

  /* Read part of a dataset from a HDF5 file. This function modifies `array` inplace.
   * 
   */
  /*
  proc readDatasetChunk(filename : string, dataset : string, offset, ref array : [] ?eltType)
      where isTuple(offset) && offset.size == array.rank {
    const rank = offset.size;
    var c_offset : [0 .. rank - 1] uint(64) = noinit;
    var c_shape : [0 .. rank - 1] uint(64) = noinit;
    for i in 0 .. rank - 1 {
      c_offset[i] = offset[i]:uint;
      c_shape[i] = array.dim(i).size:uint;
    }
    var args = (filename.localize().c_str(), dataset.localize().c_str(), rank:c_uint, c_ptrTo(c_offset), c_ptrTo(c_shape), c_ptrTo(array));
    if (eltType == uint(64)) {
      ls_hs_hdf5_read_chunk_u64((...args));
    }
    else if (eltType == real(64)) {
      ls_hs_hdf5_read_chunk_f64(filename.localize().c_str(), dataset.localize().c_str(),
        rank:c_uint, c_ptrTo(c_offset), c_ptrTo(c_shape), c_ptrTo(array));
    }
    else {
      halt("readDatasetChunk does not support " + eltType:string);
    }
  }
  */

  /* Create an HDF5 dataset of given shape and data type.
   *
   */
  /*
  proc createHDF5Dataset(filename : string, dataset : string, type eltType, shape) {
    assert(filename.locale == here && dataset.locale == here);
    var c_shape : [0 .. shape.size - 1] uint(64) = noinit;
    for i in 0 .. shape.size - 1 { c_shape[i] = shape[i]:uint; }
    if (eltType == uint(64)) {
      ls_hs_hdf5_create_dataset_u64(filename.localize().c_str(), dataset.localize().c_str(),
        c_shape.size:c_uint, c_ptrTo(c_shape));
    }
    else if (eltType == real(64)) {
      ls_hs_hdf5_create_dataset_f64(filename.localize().c_str(), dataset.localize().c_str(),
        c_shape.size:c_uint, c_ptrTo(c_shape));
    }
    else {
      assert(false);
    }
  }
  */

  /* Write array to a HDF5 dataset.
   */
  /*
  proc writeHDF5Chunk(filename : string, dataset : string, offset, array : [?D] ?eltType)
      where isTuple(offset) && offset.size == D.rank {
    assert(D.rank == offset.size);
    var c_offset : [0 .. D.rank - 1] uint(64) = noinit;
    var c_shape  : [0 .. D.rank - 1] uint(64) = noinit;
    for i in c_offset.domain { c_offset[i] = offset[i]:uint; }
    for i in c_shape.domain { c_shape[i] = array.dim(i).size:uint; }

    if (eltType == uint(64)) {
      ls_hs_hdf5_write_chunk_u64(filename.localize().c_str(), dataset.localize().c_str(),
        D.rank:c_uint, c_ptrTo(c_offset), c_ptrTo(c_shape), c_ptrTo(array));
    } else if (eltType == real(64)) {
      ls_hs_hdf5_write_chunk_f64(filename.localize().c_str(), dataset.localize().c_str(),
        D.rank:c_uint, c_ptrTo(c_offset), c_ptrTo(c_shape), c_ptrTo(array));
    }
    else {
      assert(false);
    }
  }
  */

  /*
  proc readBasisStatesAsBlocks(filename : string, dataset : string) {
    const shape = datasetShape(filename, dataset);
    if shape.size != 1 then
      halt("expected '" + dataset + "' to be one-dimensional, but it has shape " + shape:string);
    const totalNumberStates = shape[0];
  
    const box = {0 ..# totalNumberStates};
    const dom = box dmapped Block(box, Locales);
    var states : [dom] uint(64);
    coforall loc in Locales do on loc {
      const indices = states.localSubdomain();
      readDatasetChunk(filename, dataset, (indices.low,), states[indices]);
    }
    return states;
  }
  */

  proc readDatasetAsBlocks(filename : string, dataset : string, param rank = 2, type eltType = real(64)) {
    const shape = datasetShape(filename, dataset, rank);
  
    const boundingBox = _makeDomain(shape);
    const targetLocales = if rank == 1 then Locales
                                       else reshape(Locales, {0 ..# 1, 0 ..# numLocales});
    const dom = boundingBox dmapped Block(boundingBox, targetLocales);
    var vectors : [dom] eltType;
    coforall loc in Locales do on loc {
      const indices = vectors.localSubdomain();
      const offset = if rank == 1 then (indices.low,) else indices.low;
      readDatasetChunk(filename, dataset, offset, vectors[indices]);
    }
    return vectors;
  }

}
