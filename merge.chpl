module Merge {
  use List;
  use Search;

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

  inline proc _isLess(const ref a : Node, const ref b : Node,
                      const ref chunks) : bool {
    assert(!a.isEmpty() && !b.isEmpty());
    if a.isInfinity() { return false; }
    if b.isInfinity() { return true; }
    return chunks[a.chunkIndex][a.indexInChunk] < chunks[b.chunkIndex][b.indexInChunk];
  }

  proc _insert(ref tree: TournamentTree, i : int, const ref node : Node,
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
  inline proc _insert(ref tree: TournamentTree, const ref node : Node,
                      const ref chunks) : void {
    _insert(tree, node.outerParentIndex, node, chunks);
  }

  inline proc _makeNode(chunkIndex : int, in indexInChunk : int,
                        outerParentIndex : int, const ref chunks,
                        const ref counts) {
    if chunkIndex >= chunks.dim(0).size || indexInChunk >= counts[chunkIndex] {
      indexInChunk = -1;
    }
    return new Node(chunkIndex, indexInChunk, outerParentIndex);
  }

  proc _buildTree(const ref chunks, const ref counts) {
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

  proc _processOne(ref tree : TournamentTree, const ref chunks, const ref counts) {
    const result = tree.nodes[0];
    if (result.isInfinity()) { return result; }
    tree.nodes[0] = new Node(-1, -1, -1);
    const node = _makeNode(result.chunkIndex, result.indexInChunk + 1,
                           result.outerParentIndex, chunks, counts);
    _insert(tree, node, chunks);
    return result;
  }

  iter kmerge(const ref chunks, const ref counts) {
    var tree = _buildTree(chunks, counts);
    var root = _processOne(tree, chunks, counts);
    while (!root.isInfinity()) {
      yield chunks[root.chunkIndex][root.indexInChunk];
      root = _processOne(tree, chunks, counts);
    }
  }

  proc indicesForBounds(states : [] ?eltType, bounds: list(?eltType2)) {
    assert(eltType == eltType2);
    assert(states.size > 0 && bounds.size > 1);
    assert(bounds[0] <= states[0]);
    assert(bounds[bounds.size - 1] >= states[states.size - 1]);

    var edges : [0 .. bounds.size - 1] int;
    edges[0] = 0;
    edges[edges.size - 1] = states.size;
    for i in 1 .. bounds.size - 2 {
      // TODO: can be improved by setting stricted bounds where to search
      const (found, location) = binarySearch(states, bounds[i]);
      edges[i] = location;
      // writeln(i, ", ", bounds[i], ", ", found, ", ", location);
    }
    // writeln(states, ", ", bounds, " -> ", edges);
    return edges;
  }

  proc computeRanges(const ref states, const ref counts : [] int) {
    assert(states.domain.dim(0) == counts.domain.dim(0));
    const D1 = states.domain.dim(0);
    const lower = min reduce [i in D1] states[i][0];
    const upper = max reduce [i in D1] states[i][counts[i] - 1];
    const maxCount = max reduce counts;
    const chunkSize = max(counts[0] / (5 * numLocales), 1);

    var bounds : list(uint(64));
    bounds.append(lower);
    var offset = chunkSize;
    while (offset < counts[0]) {
      bounds.append(states[0][offset]);
      offset += chunkSize;
    }
    bounds.append(upper);
    writeln("bounds: ", bounds);

    var edges : [0 .. D1.size - 1, 0 .. bounds.size - 1] int;
    for i in D1 do on counts[i] {
      edges[i, ..] = indicesForBounds(states[i][0 .. counts[i] - 1], bounds);
    }

    return edges;
  }

  proc mergeStates(const ref states, const ref counts : [] int,
                   filename: string, dataset: string) {
    // Prepare the output file
    const totalCount = (+ reduce counts);
    ls_hs_hdf5_create_dataset_u64(
      filename.c_str(), dataset.c_str(), 1, c_ptrTo([totalCount:uint]));

    // Split all states into smaller chunks
    const indexEdges = computeRanges(states, counts);
    const numberRanges = indexEdges.domain.dim(1).size - 1;
    // Compute sizes of those chunks
    var offsets : [0 .. numberRanges] uint(64);
    offsets[0] = 0;
    for j in 1 .. numberRanges {
      const n = (+ reduce [i in indexEdges.dim(0)] indexEdges[i, j] - indexEdges[i, j - 1]);
      offsets[j] = offsets[j - 1] + n:uint;
    }
    // Iterate over chunks
    // TODO: parallelize
    for k in 0 .. numberRanges - 1 {
      // Prepare for kmerge
      var localCounts = [i in indexEdges.dim(0)] indexEdges[i, k + 1] - indexEdges[i, k];
      var maxLocalCount = max reduce localCounts;
      // Copy chunks from different locales to here
      var localStates : [0 .. states.dim(0).size - 1] [0 .. maxLocalCount - 1] uint(64);
      for i in 0 .. states.dim(0).size - 1 {
        localStates[i][0 .. localCounts[i] - 1] =
          states[i][indexEdges[i, k] .. indexEdges[i, k + 1] - 1];
      }
      // Run kmerge locally
      const localCount = offsets[k + 1] - offsets[k];
      assert(localCount:int == + reduce localCounts);
      var combinedLocalStates : [0 .. localCount - 1] uint(64) =
        kmerge(localStates, localCounts);

      // Save to file
      // TODO: make async
      ls_hs_hdf5_write_1d_chunk_u64(filename.c_str(), dataset.c_str(),
        offsets[k], localCount, c_ptrTo(combinedLocalStates));
    }
  }

  proc merge_test()
  {
    var chunks: [0 .. 4, 0 .. 3] int = [
      [45, 47,  0,  0],
      [0,   0,  0,  0],
      [4,  15, 16, 19],
      [17, 82,  0,  0],
      [0,   0,  0,  0],
    ];
    var counts = [2, 0, 4, 2, 1];
    for x in kmerge(chunks, counts) do
      writeln(x);
  }

} // end module
