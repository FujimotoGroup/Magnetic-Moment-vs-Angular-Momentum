#!/bin/bash
now=$(date "+%Y-%m-%d_%H%M%S")
dir="/home/jfujimoto/Work/c++/bismuth"
cmake --build build -j8
mv main main.$now
cat << EOF > run.$now.sh
#!/bin/bash
#PBS -j oe
#PBS -o $now.log
#PBS -l select=1:ncpus=32
echo $now
cd $dir
./main.$now
RETCODE=\$?
exit \${RETCODE}
EOF
qsub -N bismuth ./run.$now.sh
