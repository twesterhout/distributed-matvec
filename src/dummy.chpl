
class C {
  var dom : domain(1);
  var arr : [dom] int;

  proc init(const ref xs) {
    this.dom = xs.domain;
    this.arr = xs;
  }
}

proc main() {
  var xs = [1, 2, 3, 4];
  var c = new C(xs);
  writeln(c.arr);
  c.dom = {0 ..# 6};
  writeln(c.arr);

  c.dom = {0 ..# 2};
  writeln(c.arr);
  var cs : [0 ..# 2] owned C = [new C(xs), new C([1, 5, 7])];
}
