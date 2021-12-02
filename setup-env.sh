export PKG_CONFIG_PATH=/home/tom/src/lattice-symmetries/prefix/lib/pkgconfig
export LD_LIBRARY_PATH=/home/tom/src/lattice-symmetries/prefix/lib:/home/tom/src/lattice-symmetries-haskell/build
export CHPL_COMM_SUBSTRATE=smp CHPL_COMM=gasnet
# export CHPL_COMM=none
export CHPL_LLVM=system CHPL_UNWIND=system CHPL_LIB_PIC=pic
source "${HOME}/src/chapel-1.25.0/util/setchplenv.bash"
export CHPL_RT_OVERSUBSCRIBED=yes CHPL_RT_NUM_THREADS_PER_LOCALE=1
export OMP_NUM_THREADS=1
