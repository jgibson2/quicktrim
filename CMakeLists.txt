cmake_minimum_required(VERSION 3.7)
project(quicktrim)

set(CMAKE_C_STANDARD 11)

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wall -Wextra -pedantic -fopenmp -O3")

set(SOURCE_FILES main.c fastqrec.c fastqrec.h trim.h trim.c io.h io.c mem.h mem.c align.c align.h pipelines.h pipelines.c)
add_executable(quicktrim ${SOURCE_FILES})