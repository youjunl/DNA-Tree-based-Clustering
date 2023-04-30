# CMakeLists files in this project can
# refer to the root source directory of the project as ${HELLO_SOURCE_DIR} and
# to the root binary directory of the project as ${HELLO_BINARY_DIR}.
cmake_minimum_required(VERSION 3.24)
project(tree)
include_directories(pybind11/include)
include_directories(E:/Python310/include)
add_subdirectory(pybind11)

add_compile_options(-static-libgcc -static-libstdc++)
# Build tree module
add_library(trie trie.cpp)
target_link_libraries(trie -static-libgcc -static-libstdc++)
link_libraries(trie)

# Build python binding
pybind11_add_module(tree module.cpp)

add_executable(test_main test_main.cpp)
