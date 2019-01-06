case $COMPILER in
    gcc)
        export CC=gcc CXX=g++
        ;;
    clang)
        export CC=clang CXX=clang++
        ;;
    *)
        echo "ERROR: Invalid compiler: $COMPILER" 1>&2
        exit 1
esac

case $CUDA in
    cuda)
        NO_CUDA=FALSE
        ;;
    nocuda)
        NO_CUDA=TRUE
        ;;
    *)
        echo "ERROR: Invalid CUDA mode: $COMPILER" 1>&2
        exit 1
esac

mkdir -p build/$COMPILER-$CUDA || exit 1
cd build/$COMPILER-$CUDA || exit 1

ln -s ../../data data || exit 1

cmake -DNO_CUDA=$NO_CUDA -DCMAKE_BUILD_TYPE=Release ../.. || exit 1
make -j`nproc --all` || exit 1
