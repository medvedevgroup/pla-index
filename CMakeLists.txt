cmake_minimum_required(VERSION 3.16)

project(pla-index VERSION 1.0)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
# set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_CXX_FLAGS  "-Ofast -mbmi2 -msse4.2 -DNDEBUG -march=native")

add_library(pla STATIC 
  src/BitPacking.cpp
  src/pla_index.cpp
)
target_include_directories(pla PUBLIC ~/include/ )
target_link_directories(pla PUBLIC ~/lib/)
target_link_libraries(pla PUBLIC sdsl divsufsort divsufsort64)

add_executable(build_index 
    src/build_index.cpp    
)
target_link_libraries(build_index PUBLIC pla)

add_executable(benchmark_index
    src/benchmark_index.cpp
)
target_link_libraries(benchmark_index PUBLIC pla)

add_executable(filterNConcat
    src/FilterNConcatFasta.cpp
)

add_executable(create_queries
    src/CreateQueries.cpp
)

