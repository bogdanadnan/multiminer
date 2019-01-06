if ! [ -d pocl ]; then
    git clone https://github.com/pocl/pocl || exit 1
fi

cd pocl || exit 1
git pull || exit 1
mkdir -p build || exit 1
cd build || exit 1
if ! [ -f CMakeCache.txt ]; then
    cmake -DCMAKE_BUILD_TYPE=Release .. || exit 1
fi
make -j`nproc --all` || exit 1
