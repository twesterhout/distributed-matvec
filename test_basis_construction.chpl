use basis;
use MatVec.IO;
use wrapper only initRuntime, deinitRuntime;
use Spawn;

config const kInputBasisPath = "/home/tom/src/spin-ed/example/heisenberg_kagome_12.yaml";
  // "/home/tom/src/spin-ed/example/heisenberg_square_4x4.yaml";
  // "/home/tom/src/spin-ed/example/heisenberg_chain_4.yaml";
config const kInputDataPath = "/home/tom/src/spin-ed/example/data/heisenberg_kagome_12.h5";
// config const representativesDataset = "/basis/representatives";

config const kOutputDataPath = "output.h5";
config const kDebugLatticeSymmetries = false;

// The idea is that we load YAML data from inputBasisPath. This data is fed to
// lattice_symmetries which then constructs ls_spin_basis (we do it for each locale).
// We then build basis representatives in Chapel using DistributedBasis.
// These representatives are compared to representatives loaded from inputDataPath.
// Then dump our DistributedBasis representatives to outputDataPath. Finally we
// compare representativesDataset in inputDataPath and outputDataPath using
// h5diff.
proc main() {
  initRuntime(kDebugLatticeSymmetries);

  var basis = new DistributedBasis(kInputBasisPath);
  var states = makeStates(basis);

  var expected_representatives = loadRepresentatives(kInputDataPath, "/basis/representatives");
  for loc in Locales {
    const equal = && reduce (states.representatives[loc.id] == expected_representatives[loc.id]);
    if !equal {
      writeln("[Chapel] either loadRepresentatives or makeStates is broken");
      return 1;
    }
  }

  mergeAndWriteStates(states.representatives, states.counts, kOutputDataPath, "/representatives");
  var h5diff = spawn(
    ["h5diff", kOutputDataPath, kInputDataPath, "/representatives", "/basis/representatives"]);
  h5diff.wait();
  if h5diff.exitCode != 0 {
      writeln("[Chapel] mergeAndWriteStates is broken");
      return 1;
  }
  return 0;
}

proc deinit() {
  deinitRuntime();
}
