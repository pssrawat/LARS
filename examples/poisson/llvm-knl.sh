${CLANG} -Ofast -fopenmp=libomp -m64 -ffast-math -fstrict-aliasing -march=knl -mfma -ffp-contract=fast -mllvm -pre-RA-sched="source" $1 $2 -o test
./test

${CLANG} -Os -fopenmp=libomp -m64 -ffast-math -fstrict-aliasing -march=knl -mfma -ffp-contract=fast -mllvm -pre-RA-sched="source" $1 $2 -o test
./test

${CLANG} -Ofast -fopenmp=libomp -m64 -ffast-math -fstrict-aliasing -march=knl -mfma -ffp-contract=fast -mllvm -pre-RA-sched="list-ilp" $1 $2 -o test
./test

${CLANG} -Os -fopenmp=libomp -m64 -ffast-math -fstrict-aliasing -march=knl -mfma -ffp-contract=fast -mllvm -pre-RA-sched="list-ilp" $1 $2 -o test
./test

${CLANG} -Ofast -fopenmp=libomp -m64 -ffast-math -fstrict-aliasing -march=knl -mfma -ffp-contract=fast -mllvm -pre-RA-sched="list-hybrid" $1 $2 -o test
./test

${CLANG} -Os -fopenmp=libomp -m64 -ffast-math -fstrict-aliasing -march=knl -mfma -ffp-contract=fast -mllvm -pre-RA-sched="list-hybrid" $1 $2 -o test
./test

${CLANG} -Ofast -fopenmp=libomp -m64 -ffast-math -fstrict-aliasing -march=knl -mfma -ffp-contract=fast -mllvm -pre-RA-sched="list-burr" $1 $2 -o test
./test

${CLANG} -Os -fopenmp=libomp -m64 -ffast-math -fstrict-aliasing -march=knl -mfma -ffp-contract=fast -mllvm -pre-RA-sched="list-burr" $1 $2 -o test
./test

${CLANG} -Ofast -fopenmp=libomp -m64 -lm  -fstrict-aliasing -march=knl -mfma -ffp-contract=fast -mllvm -pre-RA-sched="source" $1 $2 -o test
./test

${CLANG} -Os -fopenmp=libomp -m64 -lm  -fstrict-aliasing -march=knl -mfma -ffp-contract=fast -mllvm -pre-RA-sched="source" $1 $2 -o test
./test

${CLANG} -Ofast -fopenmp=libomp -m64 -lm  -fstrict-aliasing -march=knl -mfma -ffp-contract=fast -mllvm -pre-RA-sched="list-ilp" $1 $2 -o test
./test

${CLANG} -Os -fopenmp=libomp -m64 -lm  -fstrict-aliasing -march=knl -mfma -ffp-contract=fast -mllvm -pre-RA-sched="list-ilp" $1 $2 -o test
./test

${CLANG} -Ofast -fopenmp=libomp -m64 -lm  -fstrict-aliasing -march=knl -mfma -ffp-contract=fast -mllvm -pre-RA-sched="list-hybrid" $1 $2 -o test
./test

${CLANG} -Os -fopenmp=libomp -m64 -lm  -fstrict-aliasing -march=knl -mfma -ffp-contract=fast -mllvm -pre-RA-sched="list-hybrid" $1 $2 -o test
./test

${CLANG} -Ofast -fopenmp=libomp -m64 -lm  -fstrict-aliasing -march=knl -mfma -ffp-contract=fast -mllvm -pre-RA-sched="list-burr" $1 $2 -o test
./test

${CLANG} -Os -fopenmp=libomp -m64 -lm  -fstrict-aliasing -march=knl -mfma -ffp-contract=fast -mllvm -pre-RA-sched="list-burr" $1 $2 -o test
./test
