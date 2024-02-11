cmake_minimum_required(VERSION 3.16)

project(pla-index VERSION 1.0)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
# set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_CXX_FLAGS  "-Ofast -mbmi2 -msse4.2 -DNDEBUG -march=native")

add_library(pla STATIC 
  src/BitPacking.cpp
  src/pla_index.cpp
  src/cmdline.cpp
)
target_include_directories(pla PUBLIC  /research/mqa6002/Softwares/Installed/include/ ) 
target_link_directories(pla PUBLIC /research/mqa6002/Softwares/Installed/lib/) 
target_link_libraries(pla PUBLIC sdsl divsufsort divsufsort64)

add_executable(build_pla_index 
    src/build_index.cpp    
)
target_link_libraries(build_pla_index PUBLIC pla)

add_executable(query_pla_index
    src/query_index.cpp
)
target_link_libraries(query_pla_index PUBLIC pla)

add_executable(process_fasta
    src/process_fasta.cpp
)

add_executable(create_queries
    src/CreateQueries.cpp
)

