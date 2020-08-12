#!/bin/bash
#
module load icc/2019
module load ifort/2019
module load impi/2019
#module load openmpi/3.1.2
#
# we compiled curl with netcdf, etc. but this module might work as well:
# module load system curl/7.54.0
echo "Modules: "
module list
#
export MIMA_LOADER="MIMA_READY"
#
# Eventually these parts should be handled in module scripts, but for now we'll do it here:
# base directories:
#SW_DIR=$SCRATCH/.local/intel/19.1.0.166
SW_DIR=.local/intel/19.1.0.166
export NETCDF_DIR=${SW_DIR}/netcdfc/4.7.3
export NETCDF_FORTRAN_DIR=${SW_DIR}/netcdff/4.5.2
export HDF5_DIR=${SW_DIR}/hdf5/1.10.6
export ZLIB_DIR=${SW_DIR}/zlib/1.2.11
export CURL_DIR=${SW_DIR}/curl/7.68.0
export MPI_DIR=$I_MPI_ROOT
#
# prepend PATH variables:
export PATH=$HDF5_DIR/bin:$CURL_DIR/bin:$ZLIB_DIR/bin:$NETCDF_DIR/bin:$NETCDF_FORTRAN_DIR/bin:MPI_DIR/bin:$PATH
export INCLUDE=$HDF5_DIR/include:$CURL_DIR/include:$ZLIB_DIR/include:$NETCDF_DIR/include:$NETCDF_FORTRAN_DIR/include:$MPI_DIR/include:$INCLUDE
export LD_LIBRARY_PATH=$HDF5_DIR/lib:$CURL_DIR/lib:$ZLIB_DIR/lib:$NETCDF_DIR/lib:$NETCDF_FORTRAN_DIR/lib:$MPI_DIR/lib:$LD_LIBRARY_PATH
#
# and add runtime path to PATH as well:
export PATH=$SCRATCH/MiMA_x/exp/exec.Sherlock:$PATH

