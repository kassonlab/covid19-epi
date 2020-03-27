project(covid19-epi)
cmake_minimum_required(VERSION 3.8)
set(CMAKE_C_STANDARD 99)
set(CMAKE_C_STANDARD_REQUIRED ON)

find_package(Threads REQUIRED)

set(COMMON_SOURCES COV_rand.c COV_rand.h COV_geo.c COV_geo.h common.h prop_distribution.c prop_distribution.h)

add_executable(covid19 covid19.c ${COMMON_SOURCES})
add_executable(sb spatial_bench.c COV_geo.c COV_geo.h)

target_link_libraries(covid19 PRIVATE Threads::Threads)
target_link_libraries(sb PRIVATE Threads::Threads)

if(UNIX)
    target_link_libraries(covid19 PRIVATE m)
    target_link_libraries(sb PRIVATE m)
endif()