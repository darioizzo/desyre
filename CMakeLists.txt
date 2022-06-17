project(sr)

cmake_minimum_required(VERSION 3.8)

find_package(fmt)

add_executable(main main.cpp)
target_link_libraries(main fmt::fmt)

set_property(TARGET main PROPERTY CXX_STANDARD 17)
set_property(TARGET main PROPERTY CXX_STANDARD_REQUIRED YES)
set_property(TARGET main PROPERTY CXX_EXTENSIONS NO)

add_executable(koza koza_quintic_simple.cpp)
target_link_libraries(koza fmt::fmt)

set_property(TARGET koza PROPERTY CXX_STANDARD 17)
set_property(TARGET koza PROPERTY CXX_STANDARD_REQUIRED YES)
set_property(TARGET koza PROPERTY CXX_EXTENSIONS NO)

add_executable(autodiff auto_diff_2nd_order.cpp)
target_link_libraries(autodiff fmt::fmt)

set_property(TARGET autodiff PROPERTY CXX_STANDARD 17)
set_property(TARGET autodiff PROPERTY CXX_STANDARD_REQUIRED YES)
set_property(TARGET autodiff PROPERTY CXX_EXTENSIONS NO)