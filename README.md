# distributed-matvec

This repository contains Chapel implementation of a few key algorithms for
lattice-symmetries package.


## Obtaining the dependencies

[`lattice-symmetries-haskell`](https://github.com/twesterhout/lattice-symmetries-haskell)
is a required dependency. You can either compile it from source (not recommended
unless you plan to contribute patches) or download a pre-built version:

```bash
LS_HS_VERSION=abe34c4
# On Linux
wget -q https://github.com/twesterhout/lattice-symmetries-haskell/releases/download/continuous/lattice-symmetries-haskell-${LS_HS_VERSION}-ubuntu-18.04.tar.bz2
tar -xf lattice-symmetries-haskell-${LS_HS_VERSION}-ubuntu-18.04.tar.bz2
# On OS X
# wget -q https://github.com/twesterhout/lattice-symmetries-haskell/releases/download/continuous/lattice-symmetries-haskell-${LS_HS_VERSION}-macos-latest.tar.bz2
# tar -xf lattice-symmetries-haskell-${LS_HS_VERSION}-macos-latest.tar.bz2
mv lattice-symmetries-haskell-${LS_HS_VERSION} third_party
```

See also the [Github Actions workflow](.github/workflows/ci.yml) for up-to-date
installation instructions (for instance, the latest tested `LS_HS_VERSION`).


## Building

Use the [Makefile](Makefile) in this repository to build the code:

```bash
make
```

This assumes that you haven an existing installation of
[Chapel](https://chapel-lang.org/).


## Testing

To run tests, type

```bash
make check
```

This will automatically download some HDF5 input files which are required to run
the tests.


## Benchmarking

There are some benchmarks which you can run using

```bash
make benchmark-states-enumeration
make benchmark-matrix-vector-product
```

Timings should be compared with the plots from the
[`lattice-symmetries`](https://github.com/twesterhout/lattice-symmetries#bicyclist-performance)
repo.

Please, note, that running benchmarks requires some extra HDF5 input files which
will use 1-2GB of storage.
