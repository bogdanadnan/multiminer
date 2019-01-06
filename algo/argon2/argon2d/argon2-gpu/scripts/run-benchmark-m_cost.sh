#!/bin/bash

work_factor="$1"; shift 1
max_memory_gb="$1"; shift 1
min_parallel="$1"; shift 1
max_parallel="$1"; shift 1
min_t_cost="$1"; shift 1
samples="$1"; shift 1
modes="$1"; shift 1
kernels="$1"; shift 1
versions="$1"; shift 1
types="$1"; shift 1
precomputes="$1"; shift 1

if [ -z "$work_factor" ]; then
    work_factor=1024
fi

if [ -z "$max_memory_gb" ]; then
    max_memory_gb=32
fi

if [ -z "$min_parallel" ]; then
    min_parallel=16
fi

if [ -z "$max_parallel" ]; then
    max_parallel=1024
fi

if [ -z "$min_t_cost" ]; then
    min_t_cost=4
fi

if [ -z "$samples" ]; then
    samples=5
fi

if [ -z "$modes" ]; then
    modes='cpu opencl cuda'
fi

if [ -z "$kernels" ]; then
    kernels='by-segment oneshot'
fi

if [ -z "$versions" ]; then
    versions='1.3 1.0'
fi

if [ -z "$types" ]; then
    types='i d id'
fi

if [ -z "$precomputes" ]; then
    precomputes='no yes'
fi

work_factor=$(( $work_factor * $max_parallel ))

echo "mode,kernel,version,type,precompute,t_cost,m_cost,lanes,ns_per_hash,batch_size"
for mode_spec in $modes; do
    case "$mode_spec" in
        *:*)
            mode=${mode_spec%%:*}
            device=${mode_spec#*:}
            device=${device:-0}
            ;;
        *)
            mode="$mode_spec"
            device=0
    esac
    
    if [ $mode != 'cpu' ]; then
        kernels2="$kernels"
    else
        kernels2='by-segment'
    fi
    for kernel in $kernels2; do
        for version in $versions; do
            for type in $types; do
                if [ $mode != 'cpu' ] && [ $type != 'd' ]; then
                    precomputes2="$precomputes"
                else
                    precomputes2='no'
                fi
                for precompute in $precomputes2; do
                    if [ $precompute == 'yes' ]; then
                        precompute_flag='-p'
                    else
                        precompute_flag=''
                    fi
                    
                    for (( lanes = 1; lanes <= 32; lanes *= 2 )); do
                        if [ $max_parallel -ge $lanes ]; then
                            batch_size=$(( $max_parallel / $lanes ))
                        else
                            batch_size=1
                        fi
                        
                        for (( m_cost = $((8 * $lanes)); ; m_cost *= 2 )); do
                            ret=1
                            while [ $(( $batch_size * $lanes )) -ge $min_parallel ] \
                                    && [ $(( $batch_size * $m_cost )) \
                                        -le $(($max_memory_gb * 1024 * 1024)) ]; do
                                t_cost=$(( $work_factor / ($m_cost * $batch_size) ))
                                if [ $t_cost -lt $min_t_cost ]; then
                                    t_cost=$min_t_cost
                                fi
                                
                                ns_per_hash=$(./argon2-gpu-bench \
                                    -t $type -v $version \
                                    $precompute_flag \
                                    -m $mode -d $device -k $kernel \
                                    -b $batch_size -s $samples \
                                    -T $t_cost -M $m_cost -L $lanes \
                                    -o ns-per-hash --output-mode mean)
                                ret=$?
                                if [ $ret -eq 0 ]; then
                                    break
                                fi
                                
                                (( batch_size /= 2 ))
                            done
                            
                            if [ $ret -ne 0 ]; then
                                break
                            fi
                            
                            echo "$mode,$kernel,v$version,Argon2$type,$precompute,$t_cost,$m_cost,$lanes,$ns_per_hash,$batch_size"
                        done
                    done
                done
            done
        done
    done
done
