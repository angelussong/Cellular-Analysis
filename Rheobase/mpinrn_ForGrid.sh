#!/bin/bash
#
#$ -cwd
#$ -j y
#$ -S /bin/bash
#$ -N mpinrn

##########################################
#                                        #
#   Output some useful job information.  #
#                                        #
##########################################

MPI_DIR=/opt/openmpi/


$MPI_DIR/bin/mpirun nrniv -mpi init.hoc
echo The Job is Done

