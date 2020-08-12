#!/bin/bash
#
#SBATCH -n 4
#SBATCH --output=mima_test_run_%j.out
#SBATCH --error=mima_test_run_%j.err
#SBATCH --partition=twohour
#
#
module purge
module unuse /usr/local/modulefiles
#
#module load MiMA/intel19-openmpi3/1.0.0
module load mima/
#
#module load intel/19.1.0.166
#module load openmpi3
#
#module load netcdf/4.7.1
#module load netcdf-fortran/4.5.2
##module load netcdf-cxx/4.3.1
##module load pnetcdf/1.12.0
#
# suppress an ethernet/infiniban MPI related warning/error:
# note:
# see:
# https://groups.io/g/OpenHPC-users/topic/mpi_errors/12423553?p=,,,20,0,0,0::recentpostdate%2Fsticky,,,20,0,0,12423553
# to further diagnose, use
# $ ompi_info
#
# so... do we need these? do we want to suppress these errors?
# these should be done in the module script
#export OMPI_MCA_mca_base_component_show_load_errors=0
#export PMIX_MCA_mca_base_component_show_load_errors=0
#
MIMA_DIR=`pwd`
#
MIMA_PLATFORM=SE3Mazama
#CCOMB=$MIMA_DIR/bin/mppnccombine.${MIMA_PLATFORM}
CCOMB=$MIMA_CCOMB
#
######################
# Begin runtinme script
#EXECDIR=${MIMA_DIR}/test_run
EXECDIR=${SCRATCH}/tmp/mima_test_run
N_PROCS=2

# ... but this will only create the leaf node, so maybe let's just require user
#  intervention.
#if directory does not exist, then create it...
[ ! -d ${EXECDIR} ] && mkdir ${EXECDIR}
#
cp -r $MIMA_DIR/input/* ${EXECDIR}/
#cp exp/exec.${MIMA_PLATFORM}/mima.x $EXECDIR/
cd $EXECDIR

[ ! -d RESTART ] && mkdir RESTART
mpiexec -n $N_PROCS mima.x

#CCOMB=/PATH/TO/MiMA/REPOSITORY/bin/mppnccombine.$PLATFORM
$CCOMB -r atmos_daily.nc atmos_daily.nc.*
$CCOMB -r atmos_avg.nc atmos_avg.nc.*
