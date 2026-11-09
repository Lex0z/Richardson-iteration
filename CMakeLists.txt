cmake_minimum_required(VERSION 3.23)
project(Richardson_iteration C)

set(CMAKE_C_STANDARD 11)
set(CMAKE_C_STANDARD_REQUIRED ON)

add_executable(Richardson_iteration
        main.c
        Richardson.c
        Richardson.h
        Matrixs.c
        Matrixs.h)

# Передаем путь к корню проекта в код как макрос PROJECT_ROOT_DIR
target_compile_definitions(${PROJECT_NAME} PRIVATE
        PROJECT_ROOT_DIR="${CMAKE_CURRENT_SOURCE_DIR}"
)

if (NOT MSVC)
    target_link_libraries(${PROJECT_NAME} PRIVATE m)
endif()
