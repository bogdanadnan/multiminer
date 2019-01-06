#!/bin/bash

machine="$1"; shift 1
machine_spec="$1"; shift 1
branches="$1"; shift 1
work_factor="$1"; shift 1
max_memory_gb="$1"; shift 1
min_threads="$1"; shift 1
max_threads="$1"; shift 1
min_t_cost="$1"; shift 1
min_m_cost="$1"; shift 1
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

if [ -z "$queue" ]; then
    queue=gpu
fi

if [ -z "$run_tests" ]; then
    run_tests='yes'
fi

dirname="$(dirname "$0")"

"$dirname/start-benchmark-common.sh" 'gpu' "$machine" "$machine_spec"  \
    "$branches" 1 "$work_factor" "$max_memory_gb" \
    "$min_threads" "$max_threads" "$min_t_cost" "$max_t_cost" "$samples" \
    "$duration" "$queue" "$run_tests"
