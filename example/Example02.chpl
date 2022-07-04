use LatticeSymmetries;
use Time;
use RangeChunk;

config const kNumSpins = 10;
config const kVerbose = false;

proc main() {
  initRuntime();
  defer deinitRuntime();

  var basis0 = new Basis("{\"number_spins\": " + kNumSpins:string
                        + ", \"hamming_weight\": " + (kNumSpins / 2):string
                        + "}");
  var r = enumerateStatesNew(basis0);
  const ref buckets0 = r[0];
  const ref masks0 = r[1];
  // var (buckets0, masks0) = ;

  if kVerbose then
    for loc in Locales do
      on loc do
        logDebug(buckets0[loc]);
  logDebug("total number of states: ", + reduce buckets0.counts);
  if kVerbose then
    logDebug("masks: ", masks0);
  return 0;
}
