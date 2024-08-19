set -e
CC="mpicc"
MPIFC="mpif90"
FC="mpifort"
CXX="mpicxx"

CORE_LDFLAGS="-L/Library/Developer/CommandLineTools/SDKs/MacOSX.sdk/usr/lib -L/opt/homebrew/lib -F/Library/Developer/CommandLineTools/SDKs/MacOSX.sdk/System/Library/Frameworks"

MPI_FFLAGS="$(pkg-config --cflags mpich) -I/opt/homebrew/lib -I/Library/Developer/CommandLineTools/SDKs/MacOSX.sdk/usr/include"
MPI_CFLAGS="$(pkg-config --cflags mpich) -I/Library/Developer/CommandLineTools/SDKs/MacOSX.sdk/usr/include"
MPI_LDFLAGS="${CORE_LDFLAGS} $(pkg-config --libs mpich) -lmpifort"
# Config Flags
MIMA_CONFIG_FFLAGS="-fallow-argument-mismatch -fdefault-real-8 -fdefault-double-8 -g ${MPI_FFLAGS} $(nf-config --fflags) $(nc-config --fflags)  $(nc-config --cflags) -cpp -fcray-pointer -fallow-invalid-boz -fno-range-check -ffree-line-length-none"
MIMA_CONFIG_CFLAGS=" -g -pthread  ${MPI_CFLAGS} $(nc-config --cflags) $(nf-config --cflags) "
MIMA_CONFIG_LDFLAGS="${MPI_LDFLAGS} -pthread $(nf-config --flibs) $(nc-config --libs) $(python3-config --ldflags)"
DEBUG="-g"
OPT="-O2"

export FFLAGS="${DEBUG} ${OPT} ${MIMA_CONFIG_FFLAGS} "
export CFLAGS=${MIMA_CONFIG_CFLAGS}
export LDFLAGS=${MIMA_CONFIG_LDFLAGS}
cwd=$(pwd)
MIMA_ROOT_PATH=$(cd ..;pwd)
#--------------------------------------------------------------------------------------------------------
# define variables
platform="apple"
template="mkmf.template.$platform"    # path to template for your platform
mkmf="${MIMA_ROOT_PATH}/bin/mkmf"                           # path to executable mkmf
sourcedir="${MIMA_ROOT_PATH}/src"                           # path to directory containing model source code
mppnccombine="${MIMA_ROOT_PATH}/bin/mppnccombine.$platform" # path to executable mppnccombine
#--------------------------------------------------------------------------------------------------------
execdir="${cwd}/exec.$platform"       # where code is compiled and executable is created
workdir="${cwd}/workdir"              # where model is run and model output is produced
pathnames="${cwd}/path_names"           # path to file containing list of source paths
namelist="${cwd}/namelists"            # path to namelist file
diagtable="${cwd}/diag_table"           # path to diagnositics table
fieldtable="${cwd}/field_table"         # path to field table (specifies tracers)
#--------------------------------------------------------------------------------------------------------
#
if [[ ! -f ${template} ]] ; then touch ${template}; fi

NETCDF_INC=$(dirname $(brew ls netcdf | grep -m1 "include"))
NETCDF_FORTRAN_INC=$(dirname $(brew ls netcdf-fortran | grep -m1 "include"))
HDF5_INC=$(dirname $(brew ls hdf5 | grep -m1 "include"))
echo $NETCDF_INC
echo $NETCDF_FORTRAN_INC
echo $HDF5_INC

NETCDF_LIB="/opt/homebrew/Cellar/netcdf/4.9.2_1/lib"
NETCDF_FORTRAN_LIB="/opt/homebrew/Cellar/netcdf-fortran/4.6.1/lib"
HDF5_LIB="/opt/homebrew/Cellar/hdf5/1.14.3_1/lib"

echo "*** compile step..."
# compile mppnccombine.c, will be used only if $npes > 1
if [[ -f "${mppnccombine}" ]]; then 
  rm ${mppnccombine}
fi
if [[ ! -f "${mppnccombine}" ]]; then
  #icc -O -o $mppnccombine -I$NETCDF_INC -L$NETCDF_LIB ${cwd}/../postprocessing/mppnccombine.c -lnetcdf
  # NOTE: this can be problematic if the SPP and MPI CC compilers get mixed up. this program often requires the spp compiler.
   ${CC} -g -o ${mppnccombine} -I/opt/homebrew/include -I${NETCDF_INC} -I${NETCDF_FORTRAN_INC} -I${HDF5_INC} -L${NETCDF_LIB} -L${NETCDF_FORTRAN_LIB} -L${HDF5_LIB}  -lnetcdf -lnetcdff ${cwd}/../postprocessing/mppnccombine.c
  echo "${mppncombine} compiled"
else
    echo "${mppnccombine} exists?"
fi

if [[ ! $? = 0 ]]; then
    echo "Something Broke! after mppnccombine..."
fi

# --------------------------------------------------------------------------------------------------------

echo "*** set up directory structure..."
# note though, we really have no busines doing anything with $workdir here, but we'll leave it to be consistent with
#  documentation.
# setup directory structure
# yoder: just brute force these. If the files/directories, exist, nuke them...
if [[ -d ${execdir} ]]; then rm -rf ${execdir}; fi
if [[ ! -d "${execdir}" ]]; then mkdir -p ${execdir}; fi
#
if [[ -e "${workdir}" ]]; then
  #echo "ERROR: Existing workdir may contaminate run. Move or remove $workdir and try again."
  #exit 1
  rm -rf ${workdir}
  mkdir -p ${workdir}
fi
#--------------------------------------------------------------------------------------------------------
echo "**"
echo "*** compile the model code and create executable"

# compile the model code and create executable
cd ${execdir}
if [[ ! $? = 0 ]]; then
    echo "Something Broke! after execdir..."
fi
#
export cppDefs="-Duse_libMPI -Duse_netCDF -DgFortran"
#
# NOTE: not sure how much of this we still need for mkmf, but this does work...
MPI_INC="/opt/homebrew/Cellar/mpich/4.2.2/include"
MPI_LIB="/opt/homebrew/Cellar/mpich/4.2.2/lib"


MPICC="mpicc"
MPICCX="mpicxx"

${mkmf} -p mima.x -t "../${template}" -c "${cppDefs}" -a $sourcedir $pathnames ${NETCDF_INC} ${NETCDF_LIB} ${NETCDF_FORTRAN_INC} ${NETCDF_FORTRAN_LIB} ${HDF5_INC} ${HDF5_LIB} ${MPI_INC} ${MPI_LIB} $sourcedir/shared/mpp/include $sourcedir/shared/include

echo "Compilers: "

echo "** CC: $CC"
echo "** CXX: $CXX"
echo "** FC: $FC"
echo "** MPICC: $MPICC"
echo "** MPICXX: $MPICXX"
echo "** MPIFC: $MPIFC"

echo $FFLAGS



make -f Makefile clean
make -f Makefile 

if [[ ! $? -eq 0 ]]; then
  echo "*** Error during make."
fi
#
echo "Exiting intentionally after make "

