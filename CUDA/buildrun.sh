cd build
module load cmake gnu cuda
cmake ../
make
# cuda-memcheck ./particleSim | tee ../out/heehee.txt
./particleSim