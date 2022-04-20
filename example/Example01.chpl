use LatticeSymmetries;
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
  basis.build();
  writeln(basis.representatives().size);
  return;

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


config const kPhysicalSystem : string = "kagome12";

proc benchmarkBuild(model : string) {
  var timer = new Timer();

  var json : string;
  if (model == "chain10") {
    json = "{\"number_spins\": 10, \"hamming_weight\": 5, \"spin_inversion\": -1,\
             \"symmetries\": [\
               {\"permutation\": [1, 2, 3, 4, 5, 6, 7, 8, 9, 0], \"sector\": 5},\
               {\"permutation\": [9, 8, 7, 6, 5, 4, 3, 2, 1, 0], \"sector\": 1}\
             ] }";
  } else if (model == "kagome12") {
    json = "{\"number_spins\": 12, \"hamming_weight\": 6}";
  }
  else if (model == "kagome16") {
    json = "{\"number_spins\": 16, \"hamming_weight\": 8}";
  }
  else if (model == "square4x4") {
    json = "{\"number_spins\": 16, \"hamming_weight\": 8, \"spin_inversion\": 1,\
             \"symmetries\": [\
               {\"permutation\": [1, 2, 3, 0, 5, 6, 7, 4, 9, 10, 11, 8, 13, 14, 15, 12], \"sector\": 0},\
               {\"permutation\": [4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 0, 1, 2, 3], \"sector\": 0},\
               {\"permutation\": [3, 2, 1, 0, 7, 6, 5, 4, 11, 10, 9, 8, 15, 14, 13, 12], \"sector\": 0},\
               {\"permutation\": [12, 13, 14, 15, 8, 9, 10, 11, 4, 5, 6, 7, 0, 1, 2, 3], \"sector\": 0},\
               {\"permutation\": [3, 7, 11, 15, 2, 6, 10, 14, 1, 5, 9, 13, 0, 4, 8, 12], \"sector\": 0}\
             ] }";
  }
  else {
    halt("invalid model: " + model);
  }

  timer.start();
  var basis = SpinBasis(json);
  timer.stop();
  writeln("Basis construction took: ", timer.elapsed(), " seconds");

  timer.clear();
  timer.start();
  basis.build();
  timer.stop();
  writeln("Building representatives took: ", timer.elapsed(), " seconds");
}


proc main() {
  ls_hs_init();

  benchmarkBuild(kPhysicalSystem);

  // testConstruction();
  // return 0;
  var basis = SpinBasis(4, 2);
  writeln(basis.numberSites());
  writeln(basis.isHammingWeightFixed());

  basis = SpinBasis("{\"number_spins\": 10, \"hamming_weight\": 5, \"spin_inversion\": -1,\
                      \"symmetries\": [\
                        {\"permutation\": [1, 2, 3, 4, 5, 6, 7, 8, 9, 0], \"sector\": 5},\
                        {\"permutation\": [9, 8, 7, 6, 5, 4, 3, 2, 1, 0], \"sector\": 1}\
                      ] }");
  var tuples = reshape(
      [0, 1,
       1, 2,
       2, 3,
       3, 4,
       4, 5,
       5, 6,
       6, 7,
       7, 8,
       8, 9,
       9, 0], {0 ..# 10, 0 ..# 2});
  var a = new Operator(basis, "sigma sigma", tuples); 
  var b = new Operator(basis, "σˣ₀ σˣ₁", tuples);
  var c = new Operator(basis, "σʸ₀ σʸ₁", tuples);
  var timer = new Timer();
  timer.start();
  var op = a + b + c;
  timer.stop();
  writeln(timer.elapsed());

  timer.clear();
  timer.start();
  basis.build();
  timer.stop();
  writeln("building representatives took ", timer.elapsed());
  // writeln(localEnumerateRepresentatives(basis, 0, 16));

  ref states = basis.representatives();
  // var x : [0 ..# 16] real(64) = [0.1030783014650637, 0.9576374437392635, 0.4373286486833481, 0.8093162494152863, 0.9922428768677899, 0.26718521544439966, 0.24723682304851302, 0.38044071869638985, 0.053766497354555076, 0.1469359882375053, 0.5787416396568855, 0.951776938627626, 0.5199039424498677, 0.7953057411693661, 0.9261507114991209, 0.700102655754242];
  // var y : [0 ..# 16] real(64);
  // var x : [0 ..# 6] real(64) = [0.1030783014650637, 0.9576374437392635, 0.4373286486833481, 0.8093162494152863, 0.9922428768677899, 0.26718521544439966];
  // var y : [0 ..# 6] real(64);
  var x : [0 ..# 16] real(64) = [0.1030783014650637, 0.9576374437392635, 0.4373286486833481, 0.8093162494152863, 0.9922428768677899, 0.26718521544439966, 0.24723682304851302, 0.38044071869638985, 0.053766497354555076, 0.1469359882375053, 0.5787416396568855, 0.951776938627626, 0.5199039424498677, 0.7953057411693661, 0.9261507114991209, 0.700102655754242];
  var y : [0 ..# 16] real(64);

  timer.clear();
  timer.start();
  localMatVec(op, states, x, y);
  timer.stop();
  writeln(y);
  writeln(timer.elapsed());
}

proc deinit() {
  ls_hs_exit();
}
