export OMP_NUM_THREADS=$1
set OMP_NUM_THREADS $1

echo "Number Threads : ${OMP_NUM_THREADS}"
./beshydro --config rhic-conf/ -o output -h
