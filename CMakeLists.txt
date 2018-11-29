cmake_minimum_required(VERSION 3.12)
project(CS4340P4)

set(CMAKE_CXX_STANDARD 17)

find_package(PythonLibs 3)

add_executable(main main.cpp)
target_include_directories(main PRIVATE ${PYTHON_INCLUDE_DIRS})
target_link_libraries(main PRIVATE ${PYTHON_LIBRARIES})
