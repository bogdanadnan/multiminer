if ! [ -d pocl ]; then
    git clone https://github.com/pocl/pocl || exit 1
fi

cd pocl/build || exit 1
make install || exit 1
mkdir -p /etc/OpenCL/vendors || exit 1
if [ -f /usr/local/etc/OpenCL/vendors/pocl.icd ]; then
    cp /usr/local/etc/OpenCL/vendors/pocl.icd \
        /etc/OpenCL/vendors/pocl.icd || exit 1
fi
cd ../.. || exit 1

cd build/$COMPILER-$CUDA || exit 1

# increase stack size limit, since pocl is very stack-hungry:
ulimit -s $(( 32 * $(ulimit -s) )) || exit 1

./argon2-gpu-test -m opencl -l || exit 1
./argon2-gpu-test -m opencl || exit 1
#CTEST_OUTPUT_ON_FAILURE=1 make test || exit 1
