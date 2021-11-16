/* k-way merge
 */
use List;

record Node {
  var chunkIndex : int;
  var indexInChunk : int;
  var outerParentIndex : int;

  inline proc isEmpty() { return chunkIndex == -1; }
  inline proc isInfinity() { return indexInChunk == -1; }
}

record TournamentTree {
  var numberChunks : int;
  var nodes : [0..#numberChunks] Node;

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

inline proc _isLess(const ref a : Node, const ref b : Node, chunks...) : bool {
  assert(!a.isEmpty() && !b.isEmpty());
  if a.isInfinity() { return false; }
  if b.isInfinity() { return true; }
  return chunks[a.chunkIndex][a.indexInChunk] < chunks[b.chunkIndex][b.indexInChunk];
}

proc _insert(ref tree: TournamentTree, i : int, const ref node : Node, chunks...) : void {
  assert(i != -1);
  assert(!node.isEmpty());
  if (tree.nodes[i].isEmpty()) { tree.nodes[i] = node; }
  else {
    const parent : int = tree.innerParentIndex(i);
    if (_isLess(node, tree.nodes[i], (...chunks))) {
      _insert(tree, parent, node, (...chunks));
    }
    else {
      const old = tree.nodes[i];
      tree.nodes[i] = node;
      _insert(tree, parent, old, (...chunks));
    }
  }
}
inline proc _insert(ref tree: TournamentTree, const ref node : Node, chunks...) : void {
  _insert(tree, node.outerParentIndex, node, (...chunks));
}

inline proc _makeNode(chunkIndex : int, in indexInChunk : int,
                      outerParentIndex : int, chunks...) {
  if chunkIndex >= chunks.size || indexInChunk >= chunks[chunkIndex].size {
    indexInChunk = -1;
  }
  return new Node(chunkIndex, indexInChunk, outerParentIndex);
}

proc _buildTree(chunks...) {
  const numberChunks = chunks.size + (chunks.size % 2);
  var tree = new TournamentTree(numberChunks);
  for i in 0 .. #numberChunks / 2 {
    const outerParentIndex = numberChunks - 1 - i;
    _insert(tree,
        _makeNode(numberChunks - 1 - 2 * i, 0, outerParentIndex, (...chunks)),
        (...chunks));
    _insert(tree,
        _makeNode(numberChunks - 2 - 2 * i, 0, outerParentIndex, (...chunks)),
        (...chunks));
  }
  return tree;
}

proc _processOne(ref tree : TournamentTree, chunks...) {
  const result = tree.nodes[0];
  if (result.isInfinity()) { return result; }
  tree.nodes[0] = new Node(-1, -1, -1);
  const node = _makeNode(result.chunkIndex, result.indexInChunk + 1,
                         result.outerParentIndex, (...chunks));
  _insert(tree, node, (...chunks));
  return result;
}

iter kmerge(chunks...) {
  var tree = _buildTree((...chunks));
  var root = _processOne(tree, (...chunks));
  while (!root.isInfinity()) {
    yield chunks[root.chunkIndex][root.indexInChunk];
    root = _processOne(tree, (...chunks));
  }
}

proc main()
{
  var chunks = (
    new list([45, 47]),
    new list(int),
    new list([4, 15, 16, 19]),
    new list([17, 82]),
    new list([0, 16]),
    // new list([1, 1, 90]),
  );
  for x in kmerge((...chunks)) do
    writeln(x);
}
