export GASNET_BACKTRACE=1
export CHPL_COMM=none
export CHPL_LLVM=system
export CHPL_UNWIND=bundled
# export CHPL_TARGET_CPU=native
export CHPL_HWLOC=bundled
export CHPL_RE2=none
export CHPL_GMP=none
export CHPL_LIB_PIC=pic
. "${HOME}/src/chapel/util/setchplenv.bash"
export LD_LIBRARY_PATH="$PWD/third_party/lib:$LD_LIBRARY_PATH"
