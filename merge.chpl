module Merge {
  use List;
  use Search;
  use SysCTypes;
  use profiling;

  record Node {
    var chunkIndex : int;
    var indexInChunk : int;
    var outerParentIndex : int;

    inline proc isEmpty() { return chunkIndex == -1; }
    inline proc isInfinity() { return indexInChunk == -1; }
  }

  record TournamentTree {
    var numberChunks : int;
    var nodes : [0 .. numberChunks - 1] Node;

    proc init(numberChunks : int)
    {
      this.numberChunks = numberChunks;
      this.nodes = [i in 0..#numberChunks] new Node(-1, -1, -1);
    }

    inline proc innerParentIndex(i : int) {
      assert(i > 0);
      return i / 2;
    }
  }

  private inline proc indexAccess(const ref chunks, i : int, j) {
    return chunks[i][j];
  }
  private inline proc indexAccess(const ref chunks : [?D] ?eltType, i : int, j)
      where (D.rank == 2) {
    return chunks[i, j];
  }

  private inline proc _isLess(const ref a : Node, const ref b : Node,
                              const ref chunks) : bool {
    assert(!a.isEmpty() && !b.isEmpty());
    if a.isInfinity() { return false; }
    if b.isInfinity() { return true; }
    return indexAccess(chunks, a.chunkIndex, a.indexInChunk)
             < indexAccess(chunks, b.chunkIndex, b.indexInChunk);
  }

  private proc _insert(ref tree: TournamentTree, i : int, const ref node : Node,
                       const ref chunks) : void {
    assert(i != -1);
    assert(!node.isEmpty());
    if (tree.nodes[i].isEmpty()) { tree.nodes[i] = node; }
    else {
      const parent : int = tree.innerParentIndex(i);
      if (_isLess(node, tree.nodes[i], chunks)) {
        _insert(tree, parent, node, chunks);
      }
      else {
        const old = tree.nodes[i];
        tree.nodes[i] = node;
        _insert(tree, parent, old, chunks);
      }
    }
  }
  private inline proc _insert(ref tree: TournamentTree, const ref node : Node,
                              const ref chunks) : void {
    _insert(tree, node.outerParentIndex, node, chunks);
  }

  private inline proc _makeNode(chunkIndex : int, in indexInChunk : int,
                                outerParentIndex : int, const ref chunks,
                                const ref counts) {
    if chunkIndex >= chunks.dim(0).size || indexInChunk >= counts[chunkIndex] {
      indexInChunk = -1;
    }
    return new Node(chunkIndex, indexInChunk, outerParentIndex);
  }

  private proc _buildTree(const ref chunks, const ref counts) {
    const numberChunks = chunks.dim(0).size + (chunks.dim(0).size % 2);
    var tree = new TournamentTree(numberChunks);
    for i in 0 .. numberChunks / 2 - 1 {
      const outerParentIndex = numberChunks - 1 - i;
      _insert(tree,
          _makeNode(numberChunks - 1 - 2 * i, 0, outerParentIndex, chunks, counts),
          chunks);
      _insert(tree,
          _makeNode(numberChunks - 2 - 2 * i, 0, outerParentIndex, chunks, counts),
          chunks);
    }
    return tree;
  }

  private proc _processOne(ref tree : TournamentTree, const ref chunks, const ref counts) {
    const result = tree.nodes[0];
    if (result.isInfinity()) { return result; }
    tree.nodes[0] = new Node(-1, -1, -1);
    const node = _makeNode(result.chunkIndex, result.indexInChunk + 1,
                           result.outerParentIndex, chunks, counts);
    _insert(tree, node, chunks);
    return result;
  }

  iter kMergeIndices(const ref chunks, const ref counts) {
    var tree = _buildTree(chunks, counts);
    var root = _processOne(tree, chunks, counts);
    while (!root.isInfinity()) {
      yield (root.chunkIndex, root.indexInChunk);
      root = _processOne(tree, chunks, counts);
    }
  }

  iter kMerge(const ref chunks, const ref counts) {
    for (chunkIndex, indexInChunk) in kMergeIndices(chunks, counts) {
      yield indexAccess(chunks, chunkIndex, indexInChunk);
    }
  }

  private proc chunkBounds(const ref arrays, const ref counts, chunkSize : int) {
    assert(chunkSize >= 1);
    assert(&& reduce [i in counts.domain] counts[i] > 0);
    type eltType = indexAccess(arrays, 0, 0).type;
    const lower = indexAccess(arrays, 0, 0);
    const upper = indexAccess(arrays, 0, counts[0] - 1);
    // const lower = min reduce [i in counts.domain] indexAccess(arrays, i, 0);
    // const upper = max reduce [i in counts.domain] indexAccess(arrays, i, counts[i] - 1);

    var bounds : list(eltType);
    bounds.append(lower);
    var offset = chunkSize;
    while (offset < counts[0]) {
      bounds.append(indexAccess(arrays, 0, offset));
      offset += chunkSize;
    }
    if (bounds.last() != upper || bounds.size < 2) { bounds.append(upper); }
    // writeln("chunkBounds: ", bounds);
    return bounds;
  }

  var _chunkedOffsetsTime = new MeasurementTable("Merge.chunkOffsets");

  private proc chunkOffsets(const ref arrays, const ref counts, chunkSize : int) {
    var __time = getTimerFor(_chunkedOffsetsTime);
    const bounds = chunkBounds(arrays, counts, chunkSize);

    var offsets : [0 ..# counts.size, 0 ..# bounds.size] int;
    forall i in counts.domain {
      offsets[i, 0] = 0;
      offsets[i, bounds.size - 1] = counts[i];
      for j in 1 .. bounds.size - 2 {
        const (_found, location) = binarySearch(
          indexAccess(arrays, i, 0 ..# counts[i]), bounds[j]);
        offsets[i, j] = location;
      }
    }
    return offsets;
  }

  private inline proc linearizeOffsets(const ref offsets) {
    return [j in offsets.dim(1)] (+ reduce offsets[.., j]);
  }

  var _chunkedIterBodyTime = new MeasurementTable("Merge.chunkedIterBody");
  var _kMergeIndicesChunked = new MeasurementTable("Merge.kMergeIndicesChunked");

  private proc chunkedIterBody(chunkId : int, const ref arrays, const ref counts,
                               const ref offsets, numberArrays, type eltType) {
    var __time = getTimerFor(_chunkedIterBodyTime);
    // Gather chunk
    const localCounts = [i in 0 ..# numberArrays] offsets[i, chunkId + 1] - offsets[i, chunkId];
    const totalLocalCount = + reduce localCounts;
    const maxLocalCount = max reduce localCounts;
    var localArrays : [0 ..# numberArrays, 0 ..# maxLocalCount] eltType;
    for i in 0 ..# numberArrays {
      localArrays[i, 0 ..# localCounts[i]] =
        indexAccess(arrays, i, offsets[i, chunkId] .. offsets[i, chunkId + 1] - 1);
    }
    // writeln("chunkId=", chunkId, " localArrays:");
    // writeln(localArrays);
    var localMergeIndices : [0 ..# totalLocalCount] (int, int);
    for ((arrayIndex, indexInArray), j) in zip(kMergeIndices(localArrays, localCounts),
                                               localMergeIndices) {
      j = (arrayIndex, offsets[arrayIndex, chunkId] + indexInArray);
    }
    return localMergeIndices;
  }

  iter kMergeIndicesChunked(const ref arrays, const ref counts, chunkSize : int) {
    var __time = getTimerFor(_kMergeIndicesChunked);
    const offsets = chunkOffsets(arrays, counts, chunkSize);
    const linear = linearizeOffsets(offsets);
    const numberArrays = counts.size;
    const numberChunks = offsets.dim(1).size - 1;
    type eltType = indexAccess(arrays, 0, 0).type;
    for chunkId in 0 ..# numberChunks {
      const indices = chunkedIterBody(chunkId, arrays, counts, offsets, numberArrays, eltType);
      yield (linear[chunkId], indices);
    }
  }
  iter kMergeIndicesChunked(param tag : iterKind,
                            const ref arrays, const ref counts, chunkSize : int)
      where tag == iterKind.standalone {
    var __time = getTimerFor(_kMergeIndicesChunked);
    const offsets = chunkOffsets(arrays, counts, chunkSize);
    const linear = linearizeOffsets(offsets);
    const numberArrays = counts.size;
    const numberChunks = offsets.dim(1).size - 1;
    type eltType = indexAccess(arrays, 0, 0).type;
    forall chunkId in 0 ..# numberChunks {
      const indices = chunkedIterBody(chunkId, arrays, counts, offsets, numberArrays, eltType);
      yield (linear[chunkId], indices);
    }
  }

  proc merge_test()
  {
    var chunks: [0 .. 4, 0 .. 3] int;
    chunks[0, ..] = [45, 47,  0,  0];
    chunks[1, ..] = [0,   0,  0,  0];
    chunks[2, ..] = [4,  15, 16, 19];
    chunks[3, ..] = [17, 82,  0,  0];
    chunks[4, ..] = [0,   0,  0,  0];
    var counts = [2, 0, 4, 2, 1];
    for x in kmerge(chunks, counts) do
      writeln(x);
  }
} // end module
