project(covid19-epi LANGUAGES C)
cmake_minimum_required(VERSION 3.13)
set(CMAKE_C_STANDARD 99)
set(CMAKE_C_STANDARD_REQUIRED ON)

include(CheckLanguage)
check_language(CUDA)

if(CMAKE_CUDA_COMPILER)
    enable_language(CUDA)
endif()

find_package(OpenMP REQUIRED)
find_package(Threads REQUIRED)

set(COMMON_SOURCES COV_rand.c COV_rand.h common.h distance.c distance.h locale.c locale.h prop_distribution.c prop_distribution.h)

add_executable(covid19 covid19.c ${COMMON_SOURCES})
target_compile_options(covid19 PRIVATE "-march=native" "-mtune=native")
target_compile_definitions(covid19 PRIVATE COV_GPU=0 USE_LOCALE_DISTANCE)
target_link_libraries(covid19 PRIVATE OpenMP::OpenMP_C Threads::Threads)

if(UNIX)
    target_link_libraries(covid19 PRIVATE m dl)
endif()

if(CMAKE_CUDA_COMPILER)
    add_executable(covid19-gpu covid19.c ${COMMON_SOURCES} COV_gpu.cu)

    target_compile_options(covid19-gpu PRIVATE
        "$<$<OR:$<COMPILE_LANGUAGE:C>,$<COMPILE_LANGUAGE:CXX>>:-march=native>"
        "$<$<OR:$<COMPILE_LANGUAGE:C>,$<COMPILE_LANGUAGE:CXX>>:-mtune=native>"
        "$<$<COMPILE_LANGUAGE:CUDA>:--compiler-options=-march=native>"
        "$<$<COMPILE_LANGUAGE:CUDA>:--compiler-options=-mtune=native>")

    target_compile_definitions(covid19-gpu PRIVATE COV_GPU=1 USE_LOCALE_DISTANCE)
    target_link_libraries(covid19-gpu PRIVATE OpenMP::OpenMP_C Threads::Threads)

    if(UNIX)
        target_link_libraries(covid19-gpu PRIVATE m dl)
    endif()
endif()
