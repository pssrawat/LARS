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
for ver in calc-dt-unroll*.c 
do
if [ -s ${ver} ]
then
	./gcc.sh calc-dt.c ${ver} >> output-gcc 
fi
done

#do the ilp version
for ver in calc-dt-ilp-*
do
if [ -s ${ver} ]
then 
	./reorder.sh ${ver}
	for i in reordered-*
	do
		./gcc.sh calc-dt.c $i >> output-gcc 
	done
fi
done

#do the acc version
for ver in calc-dt-acc-*.c
do
if [ -s ${ver} ]
then 
	./reorder.sh ${ver}
	for i in reordered-*
	do
		./gcc.sh calc-dt.c $i >> output-gcc
	done
	./gcc.sh calc-dt.c accumulate.c >> output-gcc 
	./gcc.sh calc-dt.c accumulate.c >> output-acc-gcc 
fi
done
 
#do the opt versions
for ver in calc-dt-opt-*.c
do
if [ -s ${ver} ]
then
	./reorder.sh ${ver}
	for i in reordered-*
	do
		./gcc.sh calc-dt.c $i >> output-gcc
	done
fi
done 

#LLVM
#do the unroll versions first
for ver in calc-dt-unroll*.c 
do
if [ -s ${ver} ]
then
	./llvm.sh calc-dt.c ${ver} >> output-llvm 
fi
done


#do the ilp version
for ver in calc-dt-ilp-*.c
do
if [ -s ${ver} ]
then 
	./reorder.sh ${ver}
	for i in reordered-*
	do
		./llvm.sh calc-dt.c $i >> output-llvm 
	done
fi
done

#do the acc version
for ver in calc-dt-acc-*.c
do
if [ -s ${ver} ]
then 
	./reorder.sh ${ver}
	for i in reordered-*
	do
		./llvm.sh calc-dt.c $i >> output-llvm 
	done
	./llvm.sh calc-dt.c accumulate.c >> output-llvm 
	./llvm.sh calc-dt.c accumulate.c >> output-acc-llvm 
fi
done
 
#do the opt versions
for ver in calc-dt-opt-*.c
do
if [ -s ${ver} ]
then
	./reorder.sh ${ver}
	for i in reordered-*
	do
		./llvm.sh calc-dt.c $i >> output-llvm
	done
fi
done 
