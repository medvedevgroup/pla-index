#!/bin/bash

incld_path=$1
lib_path=$2

sed -i "14s|.*|target_include_directories(pla PUBLIC  $incld_path ) |" CMakeLists.txt
sed -i "15s|.*|target_link_directories(pla PUBLIC $lib_path) |" CMakeLists.txt