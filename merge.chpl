/* k-way merge
 */
use List;


record Chunk {
  type eltType;
  var D: domain(1);
  var data: borrowed [D] eltType;

  proc init(type eltType) {
    this.eltType = eltType;
  }

  proc init(array: [?X] ?t) {
    this.eltType = t;
    this.D = X;
    this.data = array;
  }
}

record TournamentTree {
  record Node {
    var chunkIndex : int;
    var indexInChunk : int;
    var outerParentIndex : int;
    // var innerParentIndex : int;
  }

  type eltType;
  var numberChunks : int;
  var tree : [0..#numberChunks] Node;
  var chunks : [0..#numberChunks] Chunk(eltType);

  proc init(chunks: [?D] ?t)
  {
    this.eltType = t;
    this.numberChunks = D.size;
    this.tree = [i in 0..#numberChunks] new Node(-1, -1, -1);
    this.chunks = [i in 0..#numberChunks] chunks[i];
    this.complete();

    for outerParentIndex in numberChunks / 2 .. #numberChunks by -1 {
      for chunkIndex in [2 * outerParentIndex + 1, 2 * outerParentIndex] {
        _insert(outerParentIndex, new Node(chunkIndex, 0, outerParentIndex));
      }
    }
  }

  proc _innerParentIndex(i : int) { return if i == 0 then -1 else i / 2; }
  proc _isEmpty(i : int) { return this.tree[i].chunkIndex == -1; }

  proc _insert(i : int, node : Node) {
    assert(i != -1);
    if (_isEmpty(i)) {
      var innerParentIndex = if i == 0 then -1 else i / 2;
      this.tree[i] = node;
    }
    else {
      var currentKey : eltType = this.chunks[this.tree[i].chunkIndex].data[this.tree[i].indexInChunk];
      var newKey : eltType = this.chunks[node.chunkIndex].data[node.indexInChunk];
      var parent : int = _innerParentIndex(i);
      assert(parent != -1);
      // if (currentKey > newKey) {
      //   return _insert(parent, node);
      // }
      // else {
      //   var old = this.tree[i];
      //   this.tree[i] = node;
      //   return _insert(parent, old);
      // }
    }
  }
}

record Int {
  var x: int;
}

proc main()
{
  writeln("Hello world!");
  // for i in 0 .. 10 by -1 {
  //   writeln(i);
  // }
  for i in 1 .. 10 {
    writeln("upper(log2(", i, ")) = ", log2(2 * i - 1));
  }

  // var x : list(int) = list([1, 2, 3]);
  var chunks: [0..#5] shared Chunk(int) = [
    new Chunk([1, 2]),
    new Chunk([3, 4]),
    new Chunk([5, 6, 3]),
    new Chunk([7, 8]),
    new Chunk([9, 10])
  ];
  var tree = new TournamentTree(chunks);
  // writeln(chunks.type:string);
  // for x in chunks do
  //   writeln(x);

}
