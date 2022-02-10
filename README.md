## Obtaining dependencies

### Using pre-built libraries

In the project root run the following:

```sh
wget https://github.com/twesterhout/distributed-matvec/releases/download/v0.0.1-pre/third_party.tar.bz2
tar -xf third_party.tar.bz2
```

Shared libraries in `third_party/lib` depend on the following system libraries:

  * libc (GLIBC >= 2.12)
  * libgomp (GOMP >= 4.5)
  * libnuma
  * libz
  * libutil *(you probably already have it)*
  * libm *(you probably already have it)*
  * libdl *(you probably already have it)*

### Building from source

Alternatively, you can generate `third_party.tar.bz2` manually. For this you will need [Singularity](https://sylabs.io/).

1) Download GHC and HDF5 from the releases.

   ~~~~~~~~sh
   mkdir pre
   cd pre/
   wget https://github.com/twesterhout/distributed-matvec/releases/download/v0.0.0/ghc-8.10.7-x86_64-linux-ubuntu-18.04-2021-12-16.tar.xz
   wget https://github.com/twesterhout/distributed-matvec/releases/download/v0.0.0/hdf5-1.12.1.tar.bz2
   ~~~~~~~~
   (SHA256 checksums for the binaries are available on the [Releases](https://github.com/twesterhout/distributed-matvec/releases) page)

   Source distribution of HDF5 can alternatively be downloaded from the [official website](https://www.hdfgroup.org/downloads/hdf5/source-code/).

   Generating the bindist of GHC yourself is more tricky. We will be compiling Haskell into a standalone shared library and hence need GHC RTS compiled with `-fPIC`. This requires building GHC from source. So just stick to the pre-compiled GHC distributed with this repo :wink:

2) Build the containers. It might take a while.

   ~~~~~~~~sh
   make dependencies
   ~~~~~~~~


## Building

Use the [Makefile](Makefile) in this repository to build the code:

```sh
make
```


## Running

To run the code you will also need input files. They can be obtained from
[here](https://surfdrive.surf.nl/files/index.php/s/OK5527Awfgl1hT2).

Place all HDF5 files from `data/matvec` into `data/matvec` in this repository.
Files correspond to a few physical systems of different sizes.

As a small test, try running the following from the project's root folder:

```sh
./test_matvec -nl1
```

this will run matrix-vector product on a 10-spin chain which has a dimension of
13. I.e. you are multiplying a 13x13 matrix by a vector.

If you want to run a larger test, try the 5x5 square. Matrix dimension is now
5200300.

```sh
./test_matvec -nl4 \
  --kInputBasisPath data/heisenberg_square_5x5.yaml \
  --kInputDataPath data/matvec/heisenberg_square_5x5.h5 \
  --kBatchSize 10000 \
  --kProfilingThreads=false
```

For comparison, `liblattice_symmetries.so` (which is written in (mostly) C++ and
uses OpenMP for parallelism takes about 0.2 seconds on a 64-core Epyc to perform
such a matrix-vector).

Finally, you could also try running the 6x6 square although with the current
version of the code it takes way too long:

```sh
./test_matvec -nl4 \
  --kInputBasisPath data/heisenberg_square_6x6.yaml \
  --kInputDataPath data/matvec/heisenberg_square_6x6.h5 \
  --kBatchSize 10000 \
  --kProfilingThreads=false
```

(C++ + OpenMP performs this in about 12.8 seconds on a 64-core Epyc)

There is also `heisenberg_kagome_16` which is somewhere between
`heisenberg_chain_10` and `heisenberg_square_5x5` and is quite good for testing.
