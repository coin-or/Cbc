#!/bin/sh
# script to format all sources
# in subdirectories according 
# to .clang-format - 2019 - Haroldo

echo formatting all source using .clang-format
find ./ -iname '*.[ch]pp' -exec clang-format -i {} +
