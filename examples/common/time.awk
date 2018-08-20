#!/bin/bash
set -e

#GCC
if [ -s "output-gcc" ]
then
	grep -E 'orig:' output-gcc | awk ' {print $2}' > orig-gcc;
	grep -E 'unroll:' output-gcc | awk ' {print $2}' > unroll-gcc;
	grep -E 'opt:' output-gcc | awk ' {print $2}' > opt-gcc;
fi
if [ -s "output-acc-gcc" ]
then
	grep -E 'opt:' output-acc-gcc | awk ' {print $2}' > acc-gcc;
fi

echo "-------------------- GCC GFLOPS ---------------------"
if [ -s "orig-gcc" ]
then
	sort -nr orig-gcc | awk 'NR==1 {print "Original GFLOPS = " $1 }'
fi
if [ -s "unroll-gcc" ]
then
	sort -nr unroll-gcc | awk 'NR==1 {print "Unrolled GFLOPS = " $1 }'
fi
if [ -s "acc-gcc" ]
then
	sort -nr acc-gcc | awk 'NR==1 {print "Accumulation GFLOPS = " $1 }'
fi
if [ -s "opt-gcc" ]
then
	sort -nr opt-gcc | awk 'NR==1 {print "Reordered GFLOPS = " $1 }'
fi

#LLVM
if [ -s "output-llvm" ]
then
	grep -E 'orig:' output-llvm | awk ' {print $2}' > orig-llvm;
	grep -E 'unroll:' output-llvm | awk ' {print $2}' > unroll-llvm;
	grep -E 'opt:' output-llvm | awk ' {print $2}' > opt-llvm;
fi
if [ -s "output-acc-llvm" ]
then
	grep -E 'opt:' output-acc-llvm | awk ' {print $2}' > acc-llvm;
fi

echo "-------------------- LLVM GFLOPS ---------------------"
if [ -s "orig-llvm" ]
then
	sort -nr orig-llvm | awk 'NR==1 {print "Original GFLOPS = " $1 }'
fi
if [ -s "unroll-llvm" ]
then
	sort -nr unroll-llvm | awk 'NR==1 {print "Unrolled GFLOPS = " $1 }'
fi
if [ -s "acc-llvm" ]
then
	sort -nr acc-llvm | awk 'NR==1 {print "Accumulation GFLOPS = " $1 }'
fi
if [ -s "opt-llvm" ]
then
	sort -nr opt-llvm | awk 'NR==1 {print "Reordered GFLOPS = " $1 }'
fi
