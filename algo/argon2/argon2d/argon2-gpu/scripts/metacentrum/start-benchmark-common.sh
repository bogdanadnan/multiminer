#!/bin/bash

mode="$1"; shift 1
machine="$1"; shift 1
machine_spec="$1"; shift 1
branches="$1"; shift 1
ncpus="$1"; shift 1
max_memory_gb="$1"; shift 1
max_batch_size="$1"; shift 1
samples="$1"; shift 1
duration="$1"; shift 1
queue="$1"; shift 1
run_tests="$1"; shift 1

if [ -z "$mode" ]; then
    echo "ERROR: Mode must be specified!" 1>&2
    exit 1
fi

if [ -z "$machine" ]; then
    echo "ERROR: Machine must be specified!" 1>&2
    exit 1
fi

if [ -z "$machine_spec" ]; then
    echo "ERROR: Machine spec must be specified!" 1>&2
    exit 1
fi

if [ -z "$branches" ]; then
    echo "ERROR: Branches must be specified!" 1>&2
    exit 1
fi

if [ -z "$ncpus" ]; then
    echo "ERROR: Number of CPUs must be specified!" 1>&2
    exit 1
fi

if [ -z "$max_memory_gb" ]; then
    max_memory_gb=4
fi

if [ -z "$max_batch_size" ]; then
    max_batch_size=256
fi

if [ -z "$samples" ]; then
    samples=5
fi

if [ -z "$duration" ]; then
    duration=24:00:00
fi

if [ -z "$run_tests" ]; then
    run_tests='yes'
fi

REPO_URL='https://gitlab.com/omos/argon2-gpu.git'

dest_dir="$(pwd)"

task_file="$(mktemp)"

module add pbspro-client

case "$machine_spec" in
    cluster)
        machine_spec=":cl_$machine=True"
        ;;
    node:*)
        machine_spec=":vnode=${machine_spec#node:}"
        ;;
    none)
        machine_spec=""
        ;;
    *)
        echo "ERROR: Invalid machine spec!" 1>&2
        exit 1
        ;;
esac

case "$mode" in
    cpu)
        spec=''
        b_mode='cpu'
        ;;
    gpu)
        ncpus=1
        spec=':ngpus=1'
        b_mode='cuda'
        ;;
    *)
        echo "ERROR: Invalid mode!" 1>&2
        exit 1
        ;;
esac

cat >$task_file <<EOF
#!/bin/bash
#PBS -N argon2-$mode-$machine-${branches// /:}
#PBS -l select=1:ncpus=$ncpus$spec:mem=${max_memory_gb}gb$machine_spec
#PBS -l walltime=$duration
$(if [ -n "$queue" ]; then echo -n "#PBS -q $queue"; fi)

module add gcc-5.3.0
module add cmake-3.6.1
module add cuda-8.0

mkdir -p "$dest_dir/\$PBS_JOBID" || exit 1

cd "$dest_dir/\$PBS_JOBID" || exit 1

mkdir include || exit 1

git clone 'https://github.com/KhronosGroup/OpenCL-Headers.git' include/CL || exit 1
(cd include/CL && git checkout opencl11) || exit 1

export C_INCLUDE_PATH="\$C_INCLUDE_PATH:\$(readlink -f include)"
export CPLUS_INCLUDE_PATH="\$CPLUS_INCLUDE_PATH:\$(readlink -f include)"

git clone --recursive "$REPO_URL" argon2-gpu || exit 1

cd argon2-gpu || exit 1

(cmake -DCMAKE_BUILD_TYPE=Release . && make -j$ncpus) 1>../build.log 2>&1 || exit 1

if [ "$run_tests" == "yes" ]; then
    ./argon2-gpu-test 1>../tests.out 2>../tests.err
fi

bash scripts/benchmark-commits.sh "$mode-$machine" . .. "${max_memory_gb}" "$max_batch_size" "$samples" "$b_mode" '' '' '' '' $branches 1>../bench.log 2>&1
EOF

qsub "$task_file"

rm -f "$task_file"
