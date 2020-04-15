Build instructions

Prefered build method is to use CMake which can build both CPU and GPU
versions.

CMake version 3.13 or newer is required.

Using CMake
 - Make sure the intended compiler is available, supported are currently
   GNU and Intel
 - To build the GPU version the CUDA development kit must also be
   available, especially the nvcc compiler.
 - mkdir build
 - cd build
 - cmake ..
 - make

To build with the Intel compiler use
 - cmake -DCMAKE_C_COMPILER=icc -DCMAKE_CXX_COMPILER=icpc ..
 - make

The CMake file will autodetect if CUDA is available and if so build both
a CPU version, covid19, and a GPU version, covid19-gpu.

The old Makefile and Makefile.icc still works for building the CPU version.

At HPC2N the recommended compiler suite to use is the intelcuda/2019b
module, which contains CUDA, or the intel/2019b moudle for CPU-only
builds. Also load the CMake module after the compiler.

