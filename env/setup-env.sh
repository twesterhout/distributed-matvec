if [ -z "${CHPL_HOME+x}" ]; then
  export CHPL_HOME="$HOME/src/chapel"
fi
export GASNET_BACKTRACE=1
export CHPL_COMM_SUBSTRATE=smp CHPL_COMM=gasnet
export CHPL_LLVM=system
export CHPL_UNWIND=bundled
export CHPL_TARGET_CPU=native
export CHPL_HWLOC=bundled
export CHPL_RE2=none
export CHPL_GMP=none
. "${HOME}/src/chapel/util/setchplenv.bash"
export CHPL_RT_OVERSUBSCRIBED=yes
# export CHPL_RT_NUM_THREADS_PER_LOCALE=1
# export OMP_NUM_THREADS=1
export HDF5_USE_FILE_LOCKING=FALSE
export LD_LIBRARY_PATH="$PWD/third_party/lib:$LD_LIBRARY_PATH"
