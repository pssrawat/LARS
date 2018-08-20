export GCC=gcc
export CLANG=clang
export TARGET=x86
export VSIZE=256
#TARGET=knl
#VSIZE=512

cur=`pwd`
rm -f output.txt
touch output.txt

for dir in j2d* j3d* chebyshev *-point poisson helmholtz-* *-2d *-3d calc-3d rhs4th3fort hypterm diffterm  
do
	cd ${dir}
	cp gcc-${TARGET}.sh gcc.sh
	cp llvm-${TARGET}.sh llvm.sh
	cp -r opt-${TARGET}/* .
	cd ${cur}
done
 
for dir in j2d* j3d* chebyshev *-point poisson helmholtz-* *-2d *-3d calc-3d rhs4th3fort hypterm diffterm 
do
    echo $'\n\n\n=========================================================='  >> output.txt
    echo $dir >> output.txt
    echo ==========================================================  >> output.txt
    cd ${dir}
    ./run.sh > log 2>errlog  
    ../common/time.awk >> ${cur}/output.txt
    cd ${cur}
done
