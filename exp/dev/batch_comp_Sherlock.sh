#! /bin/bash
#
#SBATCH -n 1
#SBATCH -o /scratch/myoder96/Downloads/MiMA/exp/compileme.out
#SBATCH -e /scratch/myoder96/Downloads/MiMA/exp/compileme.err
#SBATCH --time=02:00:00
#
#
export MIMA_nPROCS=1
MIMA_ROOT=/scratch/myoder96/Downloads/MiMA
#cd /scratch/myoder96/Downloads/MiMA-1.0.1/exp
cd ${MIMA_ROOT}/exp
#
cd exec.SE3Mazama
make clean
rm ${MIMA_ROOT}/bin/mppnccombine.Sherlock
cd ${MIMA_ROOT}/exp
#
#rm /scratch/myoder96/MiMA/exp/compileme.out
#rm /scratch/myoder96/MiMA/exp/compileme.err
#
echo "curr dir:  `pwd`"
srun ./compilescript_Sherlock.sh

