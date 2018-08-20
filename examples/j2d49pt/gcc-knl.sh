${GCC} -Ofast -fopenmp -m64 -ffast-math -fstrict-aliasing -march=knl -mavx512f -mavx512er -mavx512cd -mavx512pf -fschedule-insns -fschedule-insns2 -fsched-pressure -mfma -ffp-contract=fast $1 $2 -o test
./test

${GCC} -Ofast -fopenmp -m64 -ffast-math -fstrict-aliasing -march=knl -mavx512f -mavx512er -mavx512cd -mavx512pf -fschedule-insns -fschedule-insns2 -mfma -ffp-contract=fast $1 $2 -o test
./test

${GCC} -Ofast -fopenmp -m64 -ffast-math -fstrict-aliasing -march=knl -mavx512f -mavx512er -mavx512cd -mavx512pf -fschedule-insns -fno-schedule-insns2 -mfma -ffp-contract=fast $1 $2 -o test
./test 

${GCC} -Ofast -fopenmp -m64 -ffast-math -fstrict-aliasing -march=knl -mavx512f -mavx512er -mavx512cd -mavx512pf -mfma -fschedule-insns2 -fno-schedule-insns -ffp-contract=fast $1 $2 -o test
./test

${GCC} -Ofast -fopenmp -m64 -ffast-math -fstrict-aliasing -march=knl -mavx512f -mavx512er -mavx512cd -mavx512pf -mfma -fno-schedule-insns -fno-schedule-insns2 -ffp-contract=fast $1 $2 -o test
./test

${GCC} -Os -fopenmp -m64 -ffast-math -fstrict-aliasing -march=knl -mavx512f -mavx512er -mavx512cd -mavx512pf -fschedule-insns -fschedule-insns2 -fsched-pressure -mfma -ffp-contract=fast $1 $2 -o test
./test

${GCC} -Os -fopenmp -m64 -ffast-math -fstrict-aliasing -march=knl -mavx512f -mavx512er -mavx512cd -mavx512pf -fschedule-insns -fschedule-insns2 -mfma -ffp-contract=fast $1 $2 -o test
./test

${GCC} -Os -fopenmp -m64 -ffast-math -fstrict-aliasing -march=knl -mavx512f -mavx512er -mavx512cd -mavx512pf -fschedule-insns -fno-schedule-insns2 -mfma -ffp-contract=fast $1 $2 -o test
./test 

${GCC} -Os -fopenmp -m64 -ffast-math -fstrict-aliasing -march=knl -mavx512f -mavx512er -mavx512cd -mavx512pf -mfma -fschedule-insns2 -fno-schedule-insns -ffp-contract=fast $1 $2 -o test
./test

${GCC} -Os -fopenmp -m64 -ffast-math -fstrict-aliasing -march=knl -mavx512f -mavx512er -mavx512cd -mavx512pf -mfma -fno-schedule-insns -fno-schedule-insns2 -ffp-contract=fast $1 $2 -o test
./test
