module MyHDF5 {
  use CTypes;
  use BlockDist;
  import HDF5;
  import HDF5.C_HDF5;

  use LatticeSymmetries.FFI;

  proc datasetRank(filename : string, dataset : string) {
    return ls_hs_hdf5_get_dataset_rank(
      filename.localize().c_str(), dataset.localize().c_str()):int;
  }

  /* Get shape of a dataset in a HDF5 file.
   *
   * :arg filename: Path to HDF5 file.
   * :arg dataset:  Path to dataset within HDF5 file.
   * 
   * :returns: shape of the dataset as a one-dimensional array.
   * :rtype: [] int
   */
  proc datasetShape(filename : string, dataset : string) {
    const file_id = C_HDF5.H5Fopen(filename.localize().c_str(),
                                   C_HDF5.H5F_ACC_RDONLY, C_HDF5.H5P_DEFAULT);
    defer C_HDF5.H5Fclose(file_id);

    const dset_id = C_HDF5.H5Dopen(file_id, dataset.localize().c_str(), C_HDF5.H5P_DEFAULT);
    defer C_HDF5.H5Dclose(dset_id);

    const dspace_id = C_HDF5.H5Dget_space(dset_id);
    defer C_HDF5.H5Sclose(dspace_id);

    const rank = C_HDF5.H5Sget_simple_extent_ndims(dspace_id);
    var c_shape : [0 ..# rank] C_HDF5.hsize_t;
    H5Sget_simple_extent_dims(dspace_id, c_ptrTo(c_shape), nil);

    return c_shape:int;
  }

  private inline proc _makeDomain(shape : 1 * int) : domain(1) {
    return {0 ..# shape[0]};
  }
  private inline proc _makeDomain(shape : 2 * int) : domain(2) {
    return {0 ..# shape[0], 0 ..# shape[1]};
  }
  private inline proc _makeDomain(shape : 3 * int) : domain(3) {
    return {0 ..# shape[0], 0 ..# shape[1], 0 ..# shape[2]};
  }
  private inline proc _makeDomain(shape : 4 * int) : domain(3) {
    return {0 ..# shape[0], 0 ..# shape[1], 0 ..# shape[2], 0 ..# shape[3]};
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
    if arr.rank != datasetRank then
      halt("rank mismatch in file: '" + filename + "' dataset: '" + dataset +
           "'  " + arr.rank:string + " != " + datasetRank:string);

    const dtype_id = C_HDF5.H5Dget_type(dset_id);
    defer C_HDF5.H5Tclose(dtype_id);
    if C_HDF5.H5Tequal(HDF5.getHDF5Type(eltType), dtype_id) <= 0 then
      halt("type mismatch in file: '" + filename + "' dataset: '" + dataset +
           "'  " + HDF5.getHDF5Type(eltType):string + " != " + dtype_id:string);
  
    var c_shape : [0 ..# arr.rank] C_HDF5.hsize_t;
    H5Sget_simple_extent_dims(dspace_id, c_ptrTo(c_shape), nil);

    const mspace_id = C_HDF5.H5Screate_simple(arr.rank, c_ptrTo(c_shape), nil);
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

  /* Write array to a HDF5 dataset.
   */
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

  proc readVectorsAsBlocks(filename : string, dataset : string, type eltType = real(64)) {
    const shape = datasetShape(filename, dataset);
    if shape.size != 2 then
      halt("expected '" + dataset + "' to be two-dimensional, but it has shape " + shape:string);
    const numberVectors = shape[0];
    const totalNumberStates = shape[1];
  
    const boundingBox = {0 ..# numberVectors, 0 ..# totalNumberStates};
    const targetLocales = reshape(Locales, {0 ..# 1, 0 ..# numLocales});
    const dom = boundingBox dmapped Block(boundingBox, targetLocales);
    var vectors : [dom] eltType;
    coforall loc in Locales do on loc {
      const indices = vectors.localSubdomain();
      readDatasetChunk(filename, dataset, indices.low, vectors[indices]);
    }
    return vectors;
  }

}
