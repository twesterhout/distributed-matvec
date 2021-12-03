#!/bin/bash

set -eu

PREFIX="/home/tom/src/spin-ed/example"

for numLocales in 1 2 3 4; do
  for yamlPath in \
    "heisenberg_chain_10.yaml" \
    "heisenberg_kagome_12.yaml" \
    "heisenberg_kagome_16.yaml" \
    "heisenberg_square_4x4.yaml" \
    "heisenberg_square_5x5.yaml" \
    "heisenberg_triangular_19.yaml"; do

    referenceOutput="$PREFIX/data/${yamlPath%.yaml}.h5"
    ./basis -nl $numLocales --yamlPath "$PREFIX/$yamlPath" --hdf5Path "output.h5"
    h5diff "output.h5" "$referenceOutput" \
      "/representatives" "/basis/representatives"
  done
done
