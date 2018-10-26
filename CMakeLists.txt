cmake_minimum_required(VERSION 3.5)


project(proyecto3 VERSION 0.0.0 LANGUAGES CXX)


add_subdirectory (src)

enable_testing ()
add_subdirectory (test)
add_subdirectory(benchmarks)
