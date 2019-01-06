#!/bin/bash

machine="$1"; shift 1
machine_spec="$1"; shift 1
branches="$1"; shift 1
max_memory_gb="$1"; shift 1
max_batch_size="$1"; shift 1
samples="$1"; shift 1
duration="$1"; shift 1
queue="$1"; shift 1
run_tests="$1"; shift 1

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

if [ -z "$max_memory_gb" ]; then
    echo "ERROR: Maximum memory must be specified!" 1>&2
    exit 1
fi

if [ -z "$max_batch_size" ]; then
    max_batch_size=256
fi

if [ -z "$samples" ]; then
    samples=4
fi

if [ -z "$queue" ]; then
    queue=gpu
fi

if [ -z "$run_tests" ]; then
    run_tests='yes'
fi

dirname="$(dirname "$0")"

"$dirname/start-benchmark-common.sh" 'gpu' "$machine" "$machine_spec" "$branches" 1 "$max_memory_gb" "$max_batch_size" "$samples" "$duration" "$queue" "$run_tests"
