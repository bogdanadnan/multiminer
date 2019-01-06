#!/bin/bash

machine="$1"; shift 1
machine_spec="$1"; shift 1
branches="$1"; shift 1
ncpus="$1"; shift 1
max_memory_gb="$1"; shift 1
max_batch_size="$1"; shift 1
samples="$1"; shift 1
duration="$1"; shift 1
queue="$1"; shift 1

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

if [ -z "$ncpus" ]; then
    echo "ERROR: Number of CPU cores must be specified!" 1>&2
    exit 1
fi

if [ -z "$samples" ]; then
    samples=8
fi

dirname="$(dirname "$0")"

"$dirname/start-benchmark-common.sh" 'cpu' "$machine" "$machine_spec" "$branches" "$ncpus" "$max_memory_gb" "$max_batch_size" "$samples" "$duration" "$queue" 'no'
