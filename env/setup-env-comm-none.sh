if [ -z "${CHPL_HOME+x}" ]; then
  export CHPL_HOME="$HOME/src/chapel"
fi
export GASNET_BACKTRACE=1
export CHPL_COMM=none
export CHPL_LLVM=system
export CHPL_TARGET_CPU=none
export CHPL_TASKS=qthreads
export CHPL_MEM=jemalloc
export CHPL_ATOMICS=cstdlib
export CHPL_UNWIND=bundled
export CHPL_HWLOC=bundled
export CHPL_RE2=none
export CHPL_GMP=none
export CHPL_LIB_PIC=pic
. "${CHPL_HOME}/util/setchplenv.bash"
export LD_LIBRARY_PATH="$PWD/third_party/lib:$LD_LIBRARY_PATH"
