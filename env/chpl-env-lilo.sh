export CHPL_COMM_SUBSTRATE=mpi CHPL_COMM=gasnet
# export CHPL_COMM=none
export CHPL_LLVM=bundled CHPL_UNWIND=bundled
export CHPL_HOME="/vol/tcm28/westerhout_tom/chapel"
export MPI_CC="$CHPL_HOME/third-party/llvm/install/linux64-x86_64-gnu/bin/clang"
# export MPI_CFLAGS="-I/vol/opt/intelcompilers/intel-2019u5/compilers_and_libraries_2019.5.281/linux/mpi/intel64/include"
# export MPI_LIBS="-L/vol/opt/intelcompilers/intel-2019u5/compilers_and_libraries_2019.5.281/linux/mpi/intel64/lib/release -L/vol/opt/intelcompilers/intel-2019u5/compilers_and_libraries_2019.5.281/linux/mpi/intel64/lib -Xlinker --enable-new-dtags -Xlinker -rpath -Xlinker /vol/opt/intelcompilers/intel-2019u5/compilers_and_libraries_2019.5.281/linux/mpi/intel64/lib/release -Xlinker -rpath -Xlinker /vol/opt/intelcompilers/intel-2019u5/compilers_and_libraries_2019.5.281/linux/mpi/intel64/lib -lmpifort -lmpi -lrt -lpthread -Wl,-z,now -Wl,-z,relro -Wl,-z,noexecstack -Xlinker --enable-new-dtags -ldl"
export MPI_CFLAGS="-I/usr/lib/x86_64-linux-gnu/openmpi/include/openmpi -I/usr/lib/x86_64-linux-gnu/openmpi/include/openmpi/opal/mca/event/libevent2022/libevent -I/usr/lib/x86_64-linux-gnu/openmpi/include/openmpi/opal/mca/event/libevent2022/libevent/include -I/usr/lib/x86_64-linux-gnu/openmpi/include"
export MPI_LIBS="-L/usr/lib -L/usr/lib/x86_64-linux-gnu/openmpi/lib -lmpi"
