rm -f accumulate* slc-* orig_* test*.c stencil*.c reordered-* output-* unroll-* orig-* acc-* opt-* test *.idsl stencils stencilnames unrollfactors printintrinsics accsize
#create files
touch orig-gcc;
touch acc-gcc;
touch unroll-gcc;
touch opt-gcc;
touch output-gcc;
touch output-acc-gcc;
touch orig-llvm;
touch acc-llvm;
touch unroll-llvm;
touch opt-llvm;
touch output-llvm;
touch output-acc-llvm;

#GCC
#do the unroll versions first
for ver in ideal-gas-unroll*.c 
do
if [ -s ${ver} ]
then
	./gcc.sh ideal-gas.c ${ver} >> output-gcc 
fi
done

#do the ilp version
for ver in ideal-gas-ilp-*
do
if [ -s ${ver} ]
then 
	./reorder.sh ${ver}
	for i in reordered-*
	do
		./gcc.sh ideal-gas.c $i >> output-gcc 
	done
fi
done

#do the acc version
for ver in ideal-gas-acc-*.c
do
if [ -s ${ver} ]
then 
	./reorder.sh ${ver}
	./gcc.sh ideal-gas.c accumulate.c >> output-acc-gcc 
	./gcc.sh ideal-gas.c accumulate.c >> output-gcc 
fi
done
 
#do the opt versions
for ver in ideal-gas-opt-*.c
do
if [ -s ${ver} ]
then
	./reorder.sh ${ver}
	for i in reordered-*
	do
		./gcc.sh ideal-gas.c $i >> output-gcc
	done
fi
done 

#LLVM
#do the unroll versions first
for ver in ideal-gas-unroll*.c 
do
if [ -s ${ver} ]
then
	./llvm.sh ideal-gas.c ${ver} >> output-llvm 
fi
done


#do the ilp version
for ver in ideal-gas-ilp-*.c
do
if [ -s ${ver} ]
then 
	./reorder.sh ${ver}
	for i in reordered-*
	do
		./llvm.sh ideal-gas.c $i >> output-llvm 
	done
fi
done

#do the acc version
for ver in ideal-gas-acc-*.c
do
if [ -s ${ver} ]
then 
	./reorder.sh ${ver}
	./llvm.sh ideal-gas.c accumulate.c >> output-acc-llvm 
	./llvm.sh ideal-gas.c accumulate.c >> output-llvm 
fi
done
 
#do the opt versions
for ver in ideal-gas-opt-*.c
do
if [ -s ${ver} ]
then
	./reorder.sh ${ver}
	for i in reordered-*
	do
		./llvm.sh ideal-gas.c $i >> output-llvm
	done
fi
done 
