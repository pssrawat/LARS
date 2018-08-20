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
for ver in accelerate-unroll*.c 
do
if [ -s ${ver} ]
then
	./gcc.sh accelerate.c ${ver} >> output-gcc 
fi
done

#do the ilp version
for ver in accelerate-ilp-*
do
if [ -s ${ver} ]
then 
	./reorder.sh ${ver}
	for i in reordered-*
	do
		./gcc.sh accelerate.c $i >> output-gcc 
	done
fi
done

#do the acc version
for ver in accelerate-acc-*.c
do
if [ -s ${ver} ]
then 
	./reorder.sh ${ver}
	./gcc.sh accelerate.c accumulate.c >> output-acc-gcc 
fi
done
 
#do the opt versions
for ver in accelerate-opt-*.c
do
if [ -s ${ver} ]
then
	./reorder.sh ${ver}
	for i in reordered-*
	do
		./gcc.sh accelerate.c $i >> output-gcc
	done
fi
done 

#LLVM
#do the unroll versions first
for ver in accelerate-unroll*.c 
do
if [ -s ${ver} ]
then
	./llvm.sh accelerate.c ${ver} >> output-llvm 
fi
done


#do the ilp version
for ver in accelerate-ilp-*.c
do
if [ -s ${ver} ]
then 
	./reorder.sh ${ver}
	for i in reordered-*
	do
		./llvm.sh accelerate.c $i >> output-llvm 
	done
fi
done

#do the acc version
for ver in accelerate-acc-*.c
do
if [ -s ${ver} ]
then 
	./reorder.sh ${ver}
	./llvm.sh accelerate.c accumulate.c >> output-acc-llvm 
fi
done
 
#do the opt versions
for ver in accelerate-opt-*.c
do
if [ -s ${ver} ]
then
	./reorder.sh ${ver}
	for i in reordered-*
	do
		./llvm.sh accelerate.c $i >> output-llvm
	done
fi
done 
