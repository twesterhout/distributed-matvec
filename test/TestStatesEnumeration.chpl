use LatticeSymmetries;
use Time;

proc localLoadRepresentatives(filename : string, dataset : string = "/representatives") {
  return readDataset(filename, dataset, uint(64), rank = 1);
}

config const kHamiltonian = "data/heisenberg_kagome_12.yaml";
config const kRepresentatives = "data/heisenberg_kagome_12.h5";

proc main() {
  initRuntime();
  defer deinitRuntime();

  const reference = localLoadRepresentatives(kRepresentatives);
  var matrix = loadHamiltonianFromYaml(kHamiltonian);
  ref basis = matrix.basis;

  var masks;
  const states = enumerateStates(basis, masks);
  for loc in Locales do
    writeln(states[loc]);

  if numLocales > 1 then return 0;

  var timer = new Timer();
  timer.start();
  basis.build();
  timer.stop();
  const ref predicted = basis.representatives();
  writeln(predicted.size);
  const theSame = && reduce [i in reference.domain] reference[i] == predicted[i];
  if !theSame {
    ref predicted = basis.representatives();
    var maxErrorCount = 10;
    for i in reference.domain {
      if reference[i] != predicted[i] && maxErrorCount > 0 {
        writeln("at index ", i, ": ", reference[i], " != ", predicted[i]);
        maxErrorCount -= 1;
      }
    }
  }
  writeln(timer.elapsed());
  return 0;
}
