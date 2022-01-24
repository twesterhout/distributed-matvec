use basis;
use MatVec.IO;
use wrapper only initRuntime, deinitRuntime;
use Spawn;

config const kInputDataPath = "/home/tom/src/spin-ed/example/data/heisenberg_kagome_12.h5";
config const kOutputDataPath = "output.h5";
config const kDebugLatticeSymmetries = false;

proc main() {
  initRuntime(kDebugLatticeSymmetries);

  const countPerLocale =
    calculateNumberStatesPerLocale(kInputDataPath, "/basis/representatives");
  var states = loadRepresentatives(kInputDataPath, "/basis/representatives", countPerLocale);
  var vectors = loadVectors(kInputDataPath,
      "/basis/representatives", "/hamiltonian/eigenvectors", real(64), countPerLocale);
  mergeAndWriteVectors(states, vectors, countPerLocale, kOutputDataPath, "/y");

  var h5diff = spawn(
    ["h5diff", kOutputDataPath, kInputDataPath, "/y", "/hamiltonian/eigenvectors"]);
  h5diff.wait();
  if h5diff.exitCode != 0 {
      writeln("[Chapel] mergeAndWriteVectors is broken");
      return 1;
  }
  return 0;
}

proc deinit() {
  deinitRuntime();
}

