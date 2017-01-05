#!/bin/bash

# CALL THIS SCRIPT WITH ARGUMENTS LIKE SO
# ./run_benchmark.sh clang
# ./run_benchmark.sh icc

DIR=../ 

if [ "$1" == "gcc" ] || [ "$1" == "g++" ]; then
    COMPILER=g++-6
    RRES=results_gcc
elif [ "$1" == "clang" ] || [ "$1" == "clang++" ]; then
    # COMPILER=clang++3.9
    COMPILER=/media/MATLAB/clang+llvm-3.9.0-x86_64-linux-gnu-ubuntu-14.04/bin/clang++
    RRES=results_clang
elif [ "$1" == "icc" ] || [ "$1" == "icpc" ]; then
    # COMPILER=/home/roman/intel_2017/bin/icpc
    COMPILER=/media/MATLAB/intel_2017/bin/icpc
    RRES=results_icc
else
    echo "Compiler not understood"
    exit 1
fi


for i in {1..6}
do
    make CXX=$COMPILER bench=scalar_variants PDEG=-DPOLYDEG=$i TCROSS=-DOPTIMISED_CROSS VECT=-DVECTORISED_CLASSIC_OVERLOADS FOLDER=$DIR >> $RRES
done
echo >> $RRES
for i in {1..6}
do
    make CXX=$COMPILER bench=scalar_variants PDEG=-DPOLYDEG=$i TCROSS=-DOPTIMISED_CROSS FOLDER=$DIR >> $RRES
done
echo >> $RRES
for i in {1..6}
do
    make CXX=$COMPILER bench=scalar_variants PDEG=-DPOLYDEG=$i FOLDER=$DIR >> $RRES
done
echo >> $RRES
for i in {1..6}
do
    make CXX=$COMPILER bench=fastor PDEG=-DPOLYDEG=$i FOLDER=$DIR >> $RRES
done

rm -rf a.out