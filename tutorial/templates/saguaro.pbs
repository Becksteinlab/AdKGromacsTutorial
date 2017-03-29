#!/bin/bash
#PBS -N AdK
#PBS -l nodes=32
#PBS -l walltime=00:20:00
#PBS -A phy598s113
#PBS -j oe
#PBS -o md.$PBS_JOBID.out

# host: saguaro
# queuing system: PBS

# max run time in hours, 1 min = 0.0167
WALL_HOURS=0.333

DEFFNM=md
TPR=$DEFFNM.tpr

LIBDIR=/home/obeckste/Library

cd $PBS_O_WORKDIR

. $LIBDIR/Gromacs/versions/4.5.5/bin/GMXRC
module load openmpi/1.4.5-intel-12.1

MDRUN=$LIBDIR/Gromacs/versions/4.5.5/bin/mdrun_mpi

# -noappend because apparently no file locking possible on Lustre (/scratch)
mpiexec $MDRUN -npme 5 -s $TPR -deffnm $DEFFNM -maxh $WALL_HOURS -cpi -noappend
