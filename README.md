## Obtaining the dependencies

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

   ```sh
   mkdir pre
   cd pre/
   wget https://github.com/twesterhout/distributed-matvec/releases/download/v0.0.0/ghc-8.10.7-x86_64-linux-ubuntu-18.04-2021-12-16.tar.xz
   wget https://github.com/twesterhout/distributed-matvec/releases/download/v0.0.0/hdf5-1.12.1.tar.bz2
   ```
   (SHA256 checksums for the binaries are available on the [Releases](https://github.com/twesterhout/distributed-matvec/releases) page)

   Source distribution of HDF5 can alternatively be downloaded from the [official website](https://www.hdfgroup.org/downloads/hdf5/source-code/).

   Generating the bindist of GHC yourself is more tricky. We will be compiling Haskell into a standalone shared library and hence need GHC RTS compiled with `-fPIC`. This requires building GHC from source. So just stick to the pre-compiled GHC distributed with this repo :wink:

2) Build the containers. It might take a while.

   ```sh
   make dependencies
   ```
