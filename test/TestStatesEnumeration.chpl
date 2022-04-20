use LatticeSymmetries;

proc localLoadRepresentatives(filename : string, dataset : string = "/basis/representatives") {
  return readDataset(filename, dataset, uint(64), rank = 1);
}

config const kBasis = "data/heisenberg_kagome_12.yaml";
config const kRepresentatives = "data/heisenberg_kagome_12.h5";

proc main() {
  ls_hs_init();

  const reference = localLoadRepresentatives(kRepresentatives);
  var basis = loadBasisFromYaml(kBasis);
  writeln("Calling basis.build() ...");
  basis.build();
  writeln("Calling basis.representatives() ...");
  assert(reference == basis.representatives());
  writeln(global_count);
}

proc deinit() {
  ls_hs_exit();
}
