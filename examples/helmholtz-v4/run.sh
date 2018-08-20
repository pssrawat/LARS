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
for ver in helmholtz-unroll*.c 
do
if [ -s ${ver} ]
then
	./gcc.sh helmholtz.c ${ver} >> output-gcc 
fi
done

#do the ilp version
for ver in helmholtz-ilp-*
do
if [ -s ${ver} ]
then 
	./reorder.sh ${ver}
	for i in reordered-*
	do
		./gcc.sh helmholtz.c $i >> output-gcc 
	done
fi
done

#do the acc version
for ver in helmholtz-acc-*.c
do
if [ -s ${ver} ]
then 
	./reorder.sh ${ver}
	./gcc.sh helmholtz.c accumulate.c >> output-acc-gcc 
fi
done
 
#do the opt versions
for ver in helmholtz-opt-*.c
do
if [ -s ${ver} ]
then
	./reorder.sh ${ver}
	for i in reordered-*
	do
		./gcc.sh helmholtz.c $i >> output-gcc
	done
fi
done 

#LLVM
#do the unroll versions first
for ver in helmholtz-unroll*.c 
do
if [ -s ${ver} ]
then
	./llvm.sh helmholtz.c ${ver} >> output-llvm 
fi
done


#do the ilp version
for ver in helmholtz-ilp-*.c
do
if [ -s ${ver} ]
then 
	./reorder.sh ${ver}
	for i in reordered-*
	do
		./llvm.sh helmholtz.c $i >> output-llvm 
	done
fi
done

#do the acc version
for ver in helmholtz-acc-*.c
do
if [ -s ${ver} ]
then 
	./reorder.sh ${ver}
	./llvm.sh helmholtz.c accumulate.c >> output-acc-llvm 
fi
done
 
#do the opt versions
for ver in helmholtz-opt-*.c
do
if [ -s ${ver} ]
then
	./reorder.sh ${ver}
	for i in reordered-*
	do
		./llvm.sh helmholtz.c $i >> output-llvm
	done
fi
done 
