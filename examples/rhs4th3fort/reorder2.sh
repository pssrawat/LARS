#!/bin/bash

awk '/#pragma begin/,/#pragma end/' $1 > stencils
sed '/#pragma begin/,/#pragma end/{//!d}' $1 > accumulate.c 
sed '/#pragma begin/,/#pragma end/{//!d}' $1 > reordered-a.c
sed '/#pragma begin/,/#pragma end/{//!d}' $1 > reordered-b.c
sed '/#pragma begin/,/#pragma end/{//!d}' $1 > reordered-c.c
sed '/#pragma begin/,/#pragma end/{//!d}' $1 > reordered-d.c
sed '/#pragma begin/,/#pragma end/{//!d}' $1 > reordered-e.c
sed '/#pragma begin/,/#pragma end/{//!d}' $1 > reordered-f.c
sed '/#pragma begin/,/#pragma end/{//!d}' $1 > reordered-g.c
sed '/#pragma begin/,/#pragma end/{//!d}' $1 > reordered-h.c
sed '/#pragma begin/,/#pragma end/{//!d}' $1 > reordered-i.c
sed '/#pragma begin/,/#pragma end/{//!d}' $1 > reordered-j.c
sed '/#pragma begin/,/#pragma end/{//!d}' $1 > reordered-k.c
sed '/#pragma begin/,/#pragma end/{//!d}' $1 > reordered-l.c
sed '/#pragma begin/,/#pragma end/{//!d}' $1 > reordered-m.c
sed '/#pragma begin/,/#pragma end/{//!d}' $1 > reordered-n.c
sed '/#pragma begin/,/#pragma end/{//!d}' $1 > reordered-o.c

awk '/#pragma begin/{print $3}' stencils > stencilnames
awk '/unroll/{print $5}' stencils > unrollfactors 
awk '/print-intrinsics/{print $7}' stencils > printintrinsics
awk '/acc-size/{print $9}' stencils > accsize 

while read -r name
do
vs=${VSIZE}
uf=`awk 'NR==1' unrollfactors`
pi=`awk 'NR==1' printintrinsics`
ac=`awk 'NR==1' accsize`
sed -i '1d' unrollfactors
sed -i '1d' printintrinsics 
sed -i '1d' accsize 
awk '/#pragma begin '"$name"'/{flag=1;next} /#pragma end '"$name"'/{flag=0} flag' stencils > $name.idsl
../../test $name.idsl --out-file $name.c --unroll $uf --distribute-rhs false --print-intrinsics false --split-size $ac --sort-function v6
sed -i '/#pragma begin '"$name"'/r '"slc-acc-$name"'.c' accumulate.c
../../test $name.idsl --out-file $name.c --unroll $uf --distribute-rhs false --vect-size $vs --print-intrinsics $pi --split-size $ac --sort-function v0
sed -i '/#pragma begin '"$name"'/r '"slc-$name"'.c' reordered-a.c
sed -i '/#pragma begin '"$name"'/r '"$name"'.c' reordered-b.c
../../test $name.idsl --out-file $name.c --unroll $uf --distribute-rhs false --vect-size $vs --print-intrinsics $pi --split-size $ac --sort-function v0-ilp
sed -i '/#pragma begin '"$name"'/r '"$name"'.c' reordered-c.c
../../test $name.idsl --out-file $name.c --unroll $uf --distribute-rhs false --vect-size $vs --print-intrinsics $pi --split-size $ac --sort-function v1
sed -i '/#pragma begin '"$name"'/r '"$name"'.c' reordered-d.c
../../test $name.idsl --out-file $name.c --unroll $uf --distribute-rhs false --vect-size $vs --print-intrinsics $pi --split-size $ac --sort-function v1-ilp
sed -i '/#pragma begin '"$name"'/r '"$name"'.c' reordered-e.c
../../test $name.idsl --out-file $name.c --unroll $uf --distribute-rhs false --vect-size $vs --print-intrinsics $pi --split-size $ac --sort-function v2
sed -i '/#pragma begin '"$name"'/r '"$name"'.c' reordered-f.c
../../test $name.idsl --out-file $name.c --unroll $uf --distribute-rhs false --vect-size $vs --print-intrinsics $pi --split-size $ac --sort-function v2-ilp
sed -i '/#pragma begin '"$name"'/r '"$name"'.c' reordered-g.c
../../test $name.idsl --out-file $name.c --unroll $uf --distribute-rhs false --vect-size $vs --print-intrinsics $pi --split-size $ac --sort-function v3
sed -i '/#pragma begin '"$name"'/r '"$name"'.c' reordered-h.c
../../test $name.idsl --out-file $name.c --unroll $uf --distribute-rhs false --vect-size $vs --print-intrinsics $pi --split-size $ac --sort-function v3-ilp
sed -i '/#pragma begin '"$name"'/r '"$name"'.c' reordered-i.c
../../test $name.idsl --out-file $name.c --unroll $uf --distribute-rhs false --vect-size $vs --print-intrinsics $pi --split-size $ac --sort-function v4
sed -i '/#pragma begin '"$name"'/r '"$name"'.c' reordered-j.c
../../test $name.idsl --out-file $name.c --unroll $uf --distribute-rhs false --vect-size $vs --print-intrinsics $pi --split-size $ac --sort-function v4-ilp
sed -i '/#pragma begin '"$name"'/r '"$name"'.c' reordered-k.c
../../test $name.idsl --out-file $name.c --unroll $uf --distribute-rhs false --vect-size $vs --print-intrinsics $pi --split-size $ac --sort-function v5
sed -i '/#pragma begin '"$name"'/r '"$name"'.c' reordered-l.c
../../test $name.idsl --out-file $name.c --unroll $uf --distribute-rhs false --vect-size $vs --print-intrinsics $pi --split-size $ac --sort-function v5-ilp
sed -i '/#pragma begin '"$name"'/r '"$name"'.c' reordered-m.c
../../test $name.idsl --out-file $name.c --unroll $uf --distribute-rhs false --vect-size $vs --print-intrinsics $pi --split-size $ac --sort-function v6
sed -i '/#pragma begin '"$name"'/r '"$name"'.c' reordered-n.c
../../test $name.idsl --out-file $name.c --unroll $uf --distribute-rhs false --vect-size $vs --print-intrinsics $pi --split-size $ac --sort-function v6-ilp
sed -i '/#pragma begin '"$name"'/r '"$name"'.c' reordered-o.c
done < stencilnames

for i in reordered-*.c accumulate.c 
do
	sed -i '/#pragma begin stencil/d' $i
	sed -i '/#pragma end stencil/d' $i
done

rm -f *~
