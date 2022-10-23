cd build

read -p "Full clean and build? " -n 1 -r
if [[ $REPLY =~ ^[Yy]$ ]]
then
    echo ""
    echo "Cleaning build files"
    rm -r CMakeFiles
    rm cmake_install.cmake CMakeCache.txt Makefile particleSim
fi

cmake ../
make
