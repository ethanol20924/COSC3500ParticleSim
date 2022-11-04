cd build
module load cmake gnu cuda
cmake ../
make
./particleSim | tee ../out/profile.txt