module Vector {

// TODO: probably shouldn't use it...
use ArrayViewSlice;

class Vector {
  type eltType;
  var _dom : domain(1);
  var _arr : [_dom] eltType;
  var _size : int;

  proc init(type t) {
    this.eltType = t;
    this._size = 0;
  }
  proc init(arr : [] ?t) {
    this.eltType = t;
    this._dom = arr.dom;
    this._arr = arr;
    this._size = arr.size;
  }

  forwarding _arr only this;
  // proc this(i : int) {
  //   if i >= _size then
  //     halt("index " + i:string + " is out of bounds for domain {0 ..# " + _size:string + "}");
  //   return _arr[i];
  // }

  proc reserve(capacity : int) {
    if capacity > _dom.size then
      _dom = {0 ..# capacity};
  }

  inline proc size { return _size; }

  proc defaultGrow(factor : real = 1.5) {
    const currentCapacity = _dom.size;
    const newCapacity =
      max(currentCapacity + 1, round(factor * currentCapacity):int);
    reserve(newCapacity);
  }

  proc pushBack(x : eltType) {
    if _size == _dom.size then
      defaultGrow();
    _arr[_size] = x;
    _size += 1;
  }

  proc append(xs : [] eltType) {
    if _size + xs.size > _dom.size then
      reserve(_size + xs.size);
    _arr[_size ..# xs.size] = xs;
    _size += xs.size;
  }
  proc append(const ref xs : Vector(eltType)) {
    append(xs._arr[0 ..# xs._size]);
  }

  proc shrink() {
    if _size < _dom.size then
      _dom = {0 ..# _size};
  }

  pragma "reference to const when const this"
  pragma "fn returns aliasing array"
  proc toArray() {
    pragma "no auto destroy" var d = {0 ..# _size};
    d._value._free_when_no_arrs = true;
    d._value.definedConst = true;
    var a = new unmanaged ArrayViewSliceArr(
        eltType=_arr.eltType,
        _DomPid=d._pid, dom=d._instance,
        _ArrPid=_arr._pid, _ArrInstance=_arr._value);
    d._value.add_arr(a, locking=false, addToList=false);
    return _newArray(a);
  }
}

proc isVector(type x : owned Vector) param { return true; }
proc isVector(type x : shared Vector) param { return true; }
proc isVector(type x : borrowed Vector) param { return true; }
proc isVector(type x : unmanaged Vector) param { return true; }
proc isVector(type x) param { return false; }

class BlockVector {
  type eltType;
  param innerRank : int;
  var _outerDom;
  var _counts : [_outerDom] int;
  var _innerDom : domain(innerRank);
  var _data : [_outerDom] [_innerDom] eltType;

  proc init(type eltType, counts : [] int, outerDom) {
    this.eltType = eltType;
    this.innerRank = 1;
    this._outerDom = outerDom;
    this._counts = counts;
    const maxNumElts = max reduce _counts;
    this._innerDom = {0 ..# maxNumElts};
  }
  proc init(type eltType, counts : [] int) {
    init(eltType, counts, counts.domain);
  }

  proc init(type eltType, batchSize : int, counts : [] int, outerDom) {
    this.eltType = eltType;
    this.innerRank = 2;
    this._outerDom = outerDom;
    this._counts = counts;
    const maxNumElts = max reduce _counts;
    this._innerDom = {0 ..# batchSize, 0 ..# maxNumElts};
  }
  proc init(type eltType, batchSize : int, counts : [] int) {
    init(eltType, batchSize, counts, counts.domain);
  }

  proc init(chunks : [] ?t)
      where chunks.domain.rank == 1 && isVector(t) {
    this.eltType = chunks[0].eltType;
    this.innerRank = 1;
    this._outerDom = chunks.domain;
    this._counts = [v in chunks] v.size;
    const maxNumElts = max reduce _counts;
    this._innerDom = {0 ..# maxNumElts};
    this.complete();
    forall (arr, count, chunk) in zip(this._data, this._counts, chunks) {
      arr[0 ..# count] = chunk.toArray();
    }
  }

  proc init(chunks : [] ?eltType, counts : [] int)
      where !isArray(eltType) {
    this.eltType = eltType;
    this.innerRank = chunks.rank - 1;
    this._outerDom = counts.domain;
    this._counts = counts;
    const maxNumElts = max reduce _counts;
    this._innerDom = {0 ..# maxNumElts};
    this.complete();
    forall (arr, count, i) in zip(this._data, this._counts, 0 ..) {
      arr[0 ..# count] = chunks[i, 0 ..# count];
    }
  }

  inline proc numBlocks { return _outerDom.size; }
  inline proc count(i) { return _counts[i]; }

  inline proc getBlockDomain(blockIdx : int)
      where innerRank == 1 {
    return {0 ..# _counts[blockIdx]};
  }
  inline proc getBlockDomain(blockIdx : int)
      where innerRank == 2 {
    return {0 ..# _innerDom.dim(0).size, 0 ..# _counts[blockIdx]};
  }

  pragma "reference to const when const this"
  pragma "fn returns aliasing array"
  proc getBlock(blockIdx : int) {
    pragma "no auto destroy" var d = getBlockDomain(blockIdx);
    d._value._free_when_no_arrs = true;
    d._value.definedConst = true;
    var a = new unmanaged ArrayViewSliceArr(
        eltType=eltType,
        _DomPid=d._pid, dom=d._instance,
        _ArrPid=_data[blockIdx]._pid, _ArrInstance=_data[blockIdx]._value);
    d._value.add_arr(a, locking=false, addToList=false);
    return _newArray(a);
  }

  pragma "reference to const when const this"
  pragma "fn returns aliasing array"
  inline proc this(loc : locale) { return getBlock(loc.id); }

  inline proc this(i : int, j) ref where innerRank == 1 {
    return _data[i][j];
  }
  inline proc this(i : int, j, k) ref where innerRank == 2 {
    return _data[i][j, k];
  }
}

} // end module Vector
