#!/bin/csh
#SBATCH -n 1
#SBATCH -o mima_compile.out
#SBATCH -e mima_compile.err
#
#Minimal runscript for atmospheric dynamical cores
# Compile script for MiMA... or a template for it at least. For more information, see docs, etc. on GitHub:
# JeffersonScientific fork (may include some site-specific tweaks):
# https://github.com/jeffersonscientific/MiMA
# or the original, parent repo:
# https://github.com/mjucker/MiMA
#
# Mazama:
# from prefix: /opt/ohpc/pub/moduledeps/intel
#
# get proper compiler (intel), mpi environment:
#module purge
#module unuse /usr/local/modulefiles
##
##module load intel
##module load mvapich2
#module load intel/19.1.0.166
##
##module load mvapich2/2.3.2
#module load openmpi3
##
#module load netcdf/4.7.1
#module load netcdf-fortran/4.5.2
##module load netcdf-cxx/4.3.1
##module load pnetcdf/1.12.0
#

#
# note also netcdf-cxx configurations:
# ncxx4-config --cflags;
# mvapich2:
#set MIMA_CONFIG_FFLAGS = "`nc-config --fflags; nf-config --fflags; pnetcdf-config --fcflags; pkg-config --cflags mvapich2` -I${HDF5_INC} -I${HDF5_LIB}"
#set MIMA_CONFIG_CFLAGS = "`nc-config --cflags; nf-config --cflags; pnetcdf-config --cflags; pkg-config --cflags mvapich2` -I${HDF5_INC}"
#set MIMA_CONFIG_LDFLAGS = " -shared-intel -L/usr/local/lib `nc-config --libs; nc-config --flibs; nf-config --flibs; ncxx4-config --libs; pnetcdf-config --ldflags; pnetcdf-config --libs; pkg-config --libs mvapich2` -L${HDF5_LIB} -L`pnetcdf-config --libdir`"

# openmpi3:
# pnetcdf-config --fflags; pnetcdf-config --fcflags;
set MIMA_CONFIG_FFLAGS = "`nf-config --cflags pkg-config --cflags ompi-fort` -I${HDF5_INC} -I${HDF5_LIB} -I${NETCDF_FORTRAN_LIB} -I${NETCDF_LIB}"
#$ pnetcdf-config --cflags ompi; pnetcdf-config --cflags; pkg-config --cflags ompi-fort;
set MIMA_CONFIG_CFLAGS = "`nf-config --cflags; pkg-config --cflags ompi`  "
# x; ncxx4-config --libs; pnetcdf-config --ldflags; pnetcdf-config --libs
#set MIMA_CONFIG_LDFLAGS = " -shared-intel -L/usr/local/lib `nc-config --libs; nc-config --flibs; nf-config --flibs; pkg-config --libs ompi; pkg-config --libs ompi-fort` -L${HDF5_LIB}"
set MIMA_CONFIG_LDFLAGS = " -shared-intel -L/usr/local/lib `nf-config --flibs; pkg-config --libs ompi-fort` "
#set LDFLAGS = " -shared-intel -L/usr/local/lib `nf-config --flibs; pkg-config --libs ompi-fort` "
set LDFLAGS = "-shared-intel -L/usr/local/lib  -L/opt/ohpc/pub/libs/intel/openmpi3/netcdf-fortran/4.5.2/lib -lnetcdff  -L/opt/ohpc/pub/libs/intel/openmpi3/hdf5/1.10.5/lib  -L/opt/ohpc/pub/libs/intel/openmpi3/netcdf/4.7.1/lib -lnetcdf -lnetcdf -lm -Wl,-rpath -Wl,/opt/ohpc/pub/mpi/openmpi3-intel/3.1.4/lib -Wl,--enable-new-dtags  -L/opt/ohpc/pub/mpi/openmpi3-intel/3.1.4/lib -lmpi_usempif08 -lmpi_usempi_ignore_tkr  -lmpi_mpifh -lmpi "
#set LDFLAGS = " -shared-intel -L/usr/local/lib `nf-config --flibs; pkg-config --libs ompi-fort` "

# TODO: in the template, LDPATH should be getting set to MIMA_CONFIG_LDFLAGS, but it appears that it is not. It is, in fact, the very final compile step
#  that appears to be failing because the linking paths are not provided... and after all of this, it looks like it may just need the -L and -l prams for
#  NetCDF and NetCDF-fortran. so we can hack it manually if necessary, but it would be nice to make the script(s) work../.
#
#
set echo
#
echo "** cwd: " $cwd
echo "MIMA_CONFIG_FFLAGS: ${MIMA_CONFIG_FFLAGS}"
echo "MIMA_CONFIG_CFLAGS: ${MIMA_CONFIG_CFLAGS}"
echo "MIMA_CONFIG_LDFLAGS: ${MIMA_CONFIG_LDFLAGS}"
#
if (! $?MIMA_nPROCS) then
	set MIMA_nPROCS = 1
endif
echo "Compile on N=${MIMA_nPROCS} process"
#
#
 
#--------------------------------------------------------------------------------------------------------
# define variables
#set platform  = nyu                                     # A unique identifier for your platform
set platform  = SE3Mazama
set npes      = $MIMA_nPROCS                                       # number of processors
set template  = $cwd/../bin/mkmf.template.$platform   # path to template for your platform
set mkmf      = $cwd/../bin/mkmf                      # path to executable mkmf
set sourcedir = $cwd/../src                           # path to directory containing model source code
set mppnccombine = $cwd/../bin/mppnccombine.$platform # path to executable mppnccombine
#--------------------------------------------------------------------------------------------------------
set execdir   = $cwd/exec.$platform       # where code is compiled and executable is created
set workdir   = $cwd/workdir              # where model is run and model output is produced
set pathnames = $cwd/path_names           # path to file containing list of source paths
set namelist  = $cwd/namelists            # path to namelist file
set diagtable = $cwd/diag_table           # path to diagnositics table
set fieldtable = $cwd/field_table         # path to field table (specifies tracers)
#--------------------------------------------------------------------------------------------------------
#
echo **
echo "*** compile step..."
# compile mppnccombine.c, will be used only if $npes > 1
if ( ! -f $mppnccombine ) then
  #icc -O -o $mppnccombine -I$NETCDF_INC -L$NETCDF_LIB $cwd/../postprocessing/mppnccombine.c -lnetcdf
    icc -O -o $mppnccombine -I$NETCDF_INC -I$NETCDF_FORTRAN_INC  -L$NETCDF_LIB -L$NETCDF_FORTRAN_LIB -lnetcdf -lnetcdff $cwd/../postprocessing/mppnccombine.c
endif
#--------------------------------------------------------------------------------------------------------

echo **
echo "*** set up directory structure..."
# setup directory structure
if ( ! -d $execdir ) mkdir $execdir
if ( -e $workdir ) then
  echo "ERROR: Existing workdir may contaminate run. Move or remove $workdir and try again."
  exit 1
endif
#--------------------------------------------------------------------------------------------------------
echo **
echo "*** compile the model code and create executable"
# compile the model code and create executable
cd $execdir
#set cppDefs = "-Duse_libMPI -Duse_netCDF"
set cppDefs = "-Duse_libMPI -Duse_netCDF -DgFortran"
#

#$mkmf -p mima.x -t $template -c "$cppDefs" -a $sourcedir $pathnames /usr/local/include $NETCDF_INC $sourcedir/shared/mpp/include $sourcedir/shared/include
$mkmf -p mima.x -t $template -c "$cppDefs" -a $sourcedir $pathnames /usr/local/include ${NETCDF_INC} ${NETCDF_FORTRAN_INC} ${HDF5_INC} ${MPI_DIR}/include $sourcedir/shared/mpp/include $sourcedir/shared/include
#
# yoder:
# this might work here -- setting LDFLAGS after all other things have passed, but it's breaking -- probably on some sort of c-shell + mkmf combined
#  syntax elsewhere. I think that the pkg-config, nc-config, nf-config, etc. calls produce some content that mkmf+cshell cannot handle, resulting in an
#  error like, "variable name must start with a letter" (so I think it's trying to define a variable in the middle of the variable definition string).
#set LDFLAGS = "${MIMA_CONFIG_LDFLAGS}"
#
#echo "Make with LDFLAGS: $LDFLAGS"
#
make -f Makefile -j $npes
