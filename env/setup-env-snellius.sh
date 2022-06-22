. env/chpl-env-snellius.sh
export LD_LIBRARY_PATH="$PWD/third_party/lib:$LD_LIBRARY_PATH"
source "${CHPL_HOME}/util/setchplenv.bash"
# export CHPL_RT_OVERSUBSCRIBED=yes
# export OMP_NUM_THREADS=1
export GASNET_PHYSMEM_MAX='167 GB'
export HDF5_USE_FILE_LOCKING=FALSE
