use LatticeSymmetries;
use Time;
use RangeChunk;

proc main() {
  initRuntime();
  defer deinitRuntime();

  var basis0 = new Basis("{\"number_spins\": 10, \"hamming_weight\": 5}");
  var buckets0 = enumerateStates(basis0, 5);

  for loc in Locales {
    writeln(buckets0[loc]);
  }
  return 0;
}
