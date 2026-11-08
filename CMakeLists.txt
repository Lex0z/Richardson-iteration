cmake_minimum_required(VERSION 3.23)
project(Richardson_iteration C)

set(CMAKE_C_STANDARD 11)

add_executable(Richardson_iteration main.c
        Richardson.c
        Richardson.h
        Matrixs.c
        Matrixs.h)
