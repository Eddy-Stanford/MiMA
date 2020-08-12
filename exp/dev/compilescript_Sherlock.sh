#!/bin/bash
#
#Minimal runscript for atmospheric dynamical cores
#
# Compile script for Stanford SRCC Sherlock HPC, copied from Stanford Earth SE3Mazama template, copied from nyu template.
# NOTE: Sherlock does not support csh, so this script has been converted to bash.
#
# get proper compiler (intel), mpi environment:
module purge
#
# SRCC Sherlock:
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
# Eventually these parts should be handled in module scripts, but for now we'll do it here:
# base directories:
SW_DIR=$SCRATCH/.local/intel/19.1.0.166
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
export LIBRARY_PATH=$LD_LIBRARY_PATH
#
export CC=icc
export FC=ifort
export LD=ifort
export CXX=icpc
#
# MPI compilers:
export MPICC=mpiicc
export MPIFC=mpiifort
export MPILD=mpiifort
export MPICXX=mpiicpc
#
export FC=$MPICC
export LD=$MPIFC
export CC=$MPICC
export CXX=$MPICXX

#
######################################################################
######################################################################
#
# shorthand to comply with earlier csh syntax:
cwd=`pwd`
#
# some convenience variables (can be defined at our discretion):
#HDF5_INC=$HDF5_DIR/include
#HDF5_LIB=$HDF5_DIR/lib
#
# openmpi3:
# pnetcdf-config --fflags; pnetcdf-config --fcflags;
# NOTE: something, I think, in the `nc-config`, etc. calls renderes these strings toxic to the end process, so (at least right now)
#   hand-code these in the template file. eventually though, we want to make these work.
export MIMA_CONFIG_FFLAGS="`nc-config --fflags; nf-config --fflags` -I${HDF5_DIR}/include -I${HDF5_DIR}/lib -I${MPI_DIR}/include  -I${MPI_DIR}/lib -I${NETCDF_FORTRAN_DIR}/include -I${NETCDF_FORTRAN_DIR}/lib"
#$ pnetcdf-config --cflags ompi; pnetcdf-config --cflags; pkg-config --cflags ompi-fort;
export MIMA_CONFIG_CFLAGS="`nc-config --cflags; nf-config --cflags` -I${MPI_DIR}/include -I${HDF5_DIR}/include"
# x; ncxx4-config --libs; pnetcdf-config --ldflags; pnetcdf-config --libs
export MIMA_CONFIG_LDFLAGS=" -shared-intel -L/usr/local/lib `nc-config --libs; nc-config --flibs; nf-config --flibs` -L${MPI_DIR}/lib -L${HDF5_DIR}/lib"
#
# TODO: in the template, LDPATH should be getting set to MIMA_CONFIG_LDFLAGS, but it appears that it is not. It is, in fact, the very final compile step
#  that appears to be failing because the linking paths are not provided... and after all of this, it looks like it may just need the -L and -l prams for
#
#set echo
#set -v
#
echo "** cwd: " $cwd
echo "MIMA_CONFIG_FFLAGS: ${MIMA_CONFIG_FFLAGS}"
echo "MIMA_CONFIG_CFLAGS: ${MIMA_CONFIG_CFLAGS}"
echo "MIMA_CONFIG_LDFLAGS: ${MIMA_CONFIG_LDFLAGS}"
#
if [ -z "$MIMA_nPROCS" ]; then
	MIMA_nPROCS=1
fi
#
echo "Compile on N=${MIMA_nPROCS} process"
#
#
#--------------------------------------------------------------------------------------------------------
# define variables
#set platform=nyu                                     # A unique identifier for your platform
export platform=Sherlock
export npes=$MIMA_nPROCS                                       # number of processors
export template=$cwd/../bin/mkmf.template.$platform   # path to template for your platform
export mkmf=$cwd/../bin/mkmf                      # path to executable mkmf
export sourcedir=$cwd/../src                           # path to directory containing model source code
export mppnccombine=$cwd/../bin/mppnccombine.$platform # path to executable mppnccombine
#--------------------------------------------------------------------------------------------------------
export execdir=$cwd/exec.$platform       # where code is compiled and executable is created
export workdir=$cwd/workdir              # where model is run and model output is produced
export pathnames=$cwd/path_names           # path to file containing list of source paths
export namelist=$cwd/namelists            # path to namelist file
export diagtable=$cwd/diag_table           # path to diagnositics table
export fieldtable=$cwd/field_table         # path to field table (specifies tracers)
#--------------------------------------------------------------------------------------------------------
#
echo "** **"
echo "*** compile step..."
# compile mppnccombine.c, will be used only if $npes > 1
if [ ! -f "$mppnccombine" ]; then
  #icc -O -o $mppnccombine -I$NETCDF_INC -L$NETCDF_LIB $cwd/../postprocessing/mppnccombine.c -lnetcdf
  #  icc -O -o $mppnccombine -I$NETCDF_INC -I$NETCDF_FORTRAN_INC  -L$NETCDF_LIB -L$NETCDF_FORTRAN_LIB -lnetcdf -lnetcdff $cwd/../postprocessing/mppnccombine.c
  icc -O -o $mppnccombine -I${NETCDF_DIR}/include -I${NETCDF_FORTRAN_DIR}/include  -L${NETCDF_DIR}/lib -L${NETCDF_FORTRAN_DIR}/lib -lnetcdf -lnetcdff $cwd/../postprocessing/mppnccombine.c
fi
#--------------------------------------------------------------------------------------------------------
#
echo "** **"
echo "*** export up directory structure..."
# setup directory structure
if [ ! -d "$execdir" ]; then
    mkdir $execdir
fi
#
if [ -e "$workdir" ]; then
  # TODO: how 'bout we just remove it?
  echo "removing old workdir: $workdir"
  rm -rf $workdir
  #echo "ERROR: Existing workdir may contaminate run. Move or remove $workdir and try again."
  #exit 1
fi
#--------------------------------------------------------------------------------------------------------
echo "** **"
echo "*** compile the model code and create executable"
# compile the model code and create executable
cd $execdir
#export cppDefs="-Duse_libMPI -Duse_netCDF"
export cppDefs="-Duse_libMPI -Duse_netCDF -DgFortran"
#

#$mkmf -p mima.x -t $template -c "$cppDefs" -a $sourcedir $pathnames /usr/local/include $NETCDF_INC $sourcedir/shared/mpp/include $sourcedir/shared/include
$mkmf -p mima.x -t $template -c "$cppDefs" -a $sourcedir $pathnames /usr/local/include ${NETCDF_DIR}/include ${NETCDF_FORTRAN_DIR}/include ${HDF5_DIR}/include ${MPI_DIR}/include $sourcedir/shared/mpp/include $sourcedir/shared/include
#
# yoder:
# this might work here -- setting LDFLAGS after all other things have passed, but it's breaking -- probably on some sort of c-shell + mkmf combined
#  syntax elsewhere. I think that the pkg-config, nc-config, nf-config, etc. calls produce some content that mkmf+cshell cannot handle, resulting in an
#  error like, "variable name must start with a letter" (so I think it's trying to define a variable in the middle of the variable definition string).
#export LDFLAGS="${MIMA_CONFIG_LDFLAGS}"

echo "Makefile created; now execute:: `pwd`"

make clean
make -f Makefile -j $npes
#make
