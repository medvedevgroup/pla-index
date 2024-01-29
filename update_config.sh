#!/bin/bash

incld_path=$1
lib_path=$2

sed -i "72s|.*|target_include_directories(salib PUBLIC src/ ext/ $incld_path \${PROJECT_BINARY_DIR}) |" strobealign/CMakeLists.txt
# sed -i "73s/.*/target_include_directories(salib PUBLIC src\/ ext\/ $lib_path \${PROJECT_BINARY_DIR}) /" strobealign/CMakeLists.txt
sed -i "73s|.*|target_link_directories(salib PUBLIC $lib_path) |" strobealign/CMakeLists.txt