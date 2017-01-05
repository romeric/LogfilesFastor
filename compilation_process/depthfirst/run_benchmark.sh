#!/bin/bash

# CALL THIS SCRIPT WITH ARGUMENTS LIKE SO
# ./run_benchmark.sh clang
# ./run_benchmark.sh icc


if [ "$1" == "gcc" ] || [ "$1" == "g++" ]; then
    COMPILER=g++-6
    CRES_OP=compilation_results_dp_gcc
    BRES_OP=binary_results_dp_gcc
    CRES_NOP=compilation_results_nodp_gcc
    BRES_NOP=binary_results_nodp_gcc
elif [ "$1" == "clang" ] || [ "$1" == "clang++" ]; then
    COMPILER=clang++-3.9
    CRES_OP=compilation_results_dp_clang
    BRES_OP=binary_results_dp_clang
    CRES_NOP=compilation_results_nodp_clang
    BRES_NOP=binary_results_nodp_clang
elif [ "$1" == "icc" ] || [ "$1" == "icpc" ]; then
    COMPILER=/home/roman/intel_2017/bin/icpc
    CRES_OP=compilation_results_dp_icc
    BRES_OP=binary_results_dp_icc
    CRES_NOP=compilation_results_nodp_icc
    BRES_NOP=binary_results_nodp_icc
else
    echo "Don't understand the compiler"
    exit 1
fi


make CXX=$COMPILER OPFLAG="OP" TIME_RES=$CRES_OP > $BRES_OP
make CXX=$COMPILER OPFLAG="NOP" TIME_RES=$CRES_NOP > $BRES_NOP
make clean