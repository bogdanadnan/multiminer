#!/bin/bash

src_dir="$1"; shift
commits="$1"; shift
command="$1"; shift

if [ -z "$command" ]; then
    echo "ERROR: Command must be specified!" 1>&2
    exit 1
fi

for commit in $commits; do
    (cd "$src_dir" && git checkout "$commit" 1>&2) || exit 1
    commit_hash=$(cd "$src_dir" && git rev-parse --verify HEAD) || exit 1
    
    make 1>&2 || exit 1
    
    "$command" "$commit" "$commit_hash" "$@"
done
