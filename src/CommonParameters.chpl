module CommonParameters {
  config const kDisplayTimings : bool = false;
  config const kHashedToBlockNumChunks = 2 * here.maxTaskPar;
  config const kBlockToHashedNumChunks = 2 * numLocales * here.maxTaskPar;
}
