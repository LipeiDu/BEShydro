rm -R output
mkdir output
rm cpu-vh
make clean
make
#./cpu-vh --config rhic-conf/ -o output -h
