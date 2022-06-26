#!/bin/bash
now=$(date "+%Y-%m-%d_%H:%M:%S")
cmake --build build -j8
mv main main.$now
cat << EOF > run.$now
#!/bin/bash
#SBATCH -p haku0
#SBATCH -n 1
#SBATCH -c 48
#SBATCH -J L
#SBATCH -o stdout.%J
#SBATCH -e stderr.%J
echo $now
./main.$now
RETCODE=\$?
exit \${RETCODE}
EOF
sbatch run.$now
