#!/bin/bash
#SBATCH -p haku1
#SBATCH -n 1
#SBATCH -c 48
#SBATCH -J check
#SBATCH -o stdout.%J
#SBATCH -e stderr.%J
./main
RETCODE=$?
exit ${RETCODE}
