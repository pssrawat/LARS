DESCRIPTION

The artifact describes LARS, an automated framework to perform reordering
optimization on straight-line codes; the algorithmic details of the framework
are described in the SC'18 paper "Associative Instruction Reordering to
Alleviate Register Pressure". A copy of the paper is in the main directory
(sc-18.pdf).

The downloaded package comes with
	a. item The source code for LARS 
	b. the benchmarks in the examples/ directory
	c. scripts for code installation and benchmarking




DEPENDENCIES

We tested LARS on ubuntu 16.04 and Red Hat Enterprise Linux Server release 6.7
using GCC 5.3.0, LLVM 5.0.0. The following are hardware requirements for
compiling LARS:
1. flex >= 2.6.0
2. bison >= 3.0.4
3. cmake >= 3.8
4. boost >=1.58
5. GCC >=4.9 with c++11 support for compiling

For benchmarking:
6. GCC >= 7.2.0
7. LLVM >= 5.0




STEPS TO INSTALL AND REPRODUCE RESULTS

1. Download and install LLVM. If you cannot get the latest version of LLVM from apt or any other repo, you can download it from http://releases.llvm.org/download.html. The installation 
   steps are in https://llvm.org/docs/GettingStarted.html. 
   I downloaded LLVM into source directory, and created two separate build and install directories. Then, from
   the build directory, I used the following command to configure LLVM:
	cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=`pwd`/../install ../source/llvm/ -DLLVM_TARGETS_TO_BUILD="X86" -DGCC_INSTALL_PREFIX=<path-to-gcc-install> -DCMAKE_C_COMPILER=<path-to-gcc-install>/bin/gcc -DCMAKE_CXX_COMPILER=<path-to-gcc-install>/bin/g++ 
    You may need to adjust the paths according to your machine configuration.

2. Set the PATH and LD_LIBRARY_PATH variables to point to LLVM installation
	export PATH=<path-to-llvm-install>/include:$PATH
	export LD_LIBRARY_PATH=<path-to-llvm-install>/lib:$LD_LIBRARY_PATH
 
3. Now go to the LARS download, and simply run 'make all' in the main directory. The makefile will create the 'test' executable.

2. Set the paths to benchmarking gcc and clang in the run-benchmarks.sh script. Also, set the target (multi-core CPU or Xeon Phi).

 *Some scripts will not run if the paths are not set*

3. Go to the examples directory, and run the benchmarking script as './run-benchmarks.sh'.
   This will create a file 'output.txt' with all the results. Alternatively, you can set the paths to the benchmarking
   gcc and clang (first 4 lines in run-benchmarks.sh), then go to an independent benchmark
   directory, execute 'run.sh', and see the printed results on standard output. 

   Note that the benchmarks may take >40 hours to run. 




COPYRIGHT

All files in this archive which do not include a prior copyright are by default
included in this tool and copyrighted 2018 Ohio State University.




MORE INFORMATION

For more information on how to add a new benchmark, see the docs/ folder or contact me at <rawat.15@osu.edu>
