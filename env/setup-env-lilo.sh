# . env/chpl-env-lilo.sh
. env/chpl-env-lilo-ofi.sh
export LD_LIBRARY_PATH="$PWD/third_party/lib:$LD_LIBRARY_PATH"
source "${CHPL_HOME}/util/setchplenv.bash"
# export CHPL_RT_OVERSUBSCRIBED=yes
# export OMP_NUM_THREADS=1
export HDF5_USE_FILE_LOCKING=FALSE

run_with_mpi() {
    APP="$1"
    shift

    case "$1" in
        -nl)
            NUM_LOCALES=$2
            shift
            shift
            ;;
        -nl*)
            NUM_LOCALES="${1#-nl}"
            shift
            ;;
        *)
            echo "Expected -nl argument specifying the number of locales"
            return 1
    esac
    mpirun -np $NUM_LOCALES --bind-to none -mca btl_base_warn_component_unused 0 -- \
        "$APP" -nl$NUM_LOCALES "$@"
}

example05() {
    run_with_mpi "bin/Example05_real" "$@"
}
