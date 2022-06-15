project(sr)

cmake_minimum_required(VERSION 3.8)

find_package(fmt)

add_executable(main main.cpp)
target_link_libraries(main fmt::fmt)

set_property(TARGET main PROPERTY CXX_STANDARD 17)
set_property(TARGET main PROPERTY CXX_STANDARD_REQUIRED YES)
set_property(TARGET main PROPERTY CXX_EXTENSIONS NO)