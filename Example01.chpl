use ApplyOperator;
use StatesEnumeration;
use Time;

proc bas() {
}

proc testConstruction() {
  // var basis = SpinBasis(4);
  // writeln(localEnumerateRepresentatives(basis));

  // basis = SpinlessFermionicBasis(4, 1);
  // writeln(localEnumerateRepresentatives(basis));

  // const basis = SpinfulFermionicBasis(numberSites=3, numberUp=2, numberDown=1);
  // const representatives = localEnumerateRepresentatives(basis);

  const basis = SpinChain10();
  writeln(localEnumerateRepresentatives(basis));

  var alphas = [
    31,
    47,
    55,
    28,
    30
  ]:uint(64);
  var (areRepresentatives, norms) = isRepresentative(basis, alphas);
  writeln(areRepresentatives);
  writeln(norms);
}



proc main() {
  ls_hs_init();

  testConstruction();
  return 0;

  var basis = SpinBasis(4); // ls_hs_create_spin_basis(10, -1);
  var tuples = reshape(
      [0, 1,
       1, 2,
       2, 3,
       3, 0], {0 ..# 4, 0 ..# 2});
  var a = new Operator(basis, "σᶻ₀ σᶻ₁", tuples); 
  var b = new Operator(basis, "σˣ₀ σˣ₁", tuples);
  var c = new Operator(basis, "σʸ₀ σʸ₁", tuples);
  var timer = new Timer();
  timer.start();
  var op = a + b + c;
  timer.stop();
  writeln(timer.elapsed());

  writeln(localEnumerateRepresentatives(basis, 0, 16));

  var states : [0 ..# 16] uint(64) = [i in 0 ..# 16] i:uint;
  var x : [0 ..# 16] real(64) = [0.1030783014650637, 0.9576374437392635, 0.4373286486833481, 0.8093162494152863, 0.9922428768677899, 0.26718521544439966, 0.24723682304851302, 0.38044071869638985, 0.053766497354555076, 0.1469359882375053, 0.5787416396568855, 0.951776938627626, 0.5199039424498677, 0.7953057411693661, 0.9261507114991209, 0.700102655754242];
  var y : [0 ..# 16] real(64);

  localMatVec(op, states, x, y);
  writeln(y);

  writeln("Hello world!");
}

proc deinit() {
  ls_hs_exit();
}
