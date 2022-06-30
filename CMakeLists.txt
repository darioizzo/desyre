# NOTE: C++17 supported since CMake 3.8.0:
# https://cmake.org/cmake/help/v3.8/prop_tgt/CXX_STANDARD.html
cmake_minimum_required(VERSION 3.8.0)

# Set default build type to "Release".
# NOTE: this should be done before the project command since the latter can set
# CMAKE_BUILD_TYPE itself (it does so for nmake).
if(NOT CMAKE_BUILD_TYPE)
	set(CMAKE_BUILD_TYPE Release CACHE STRING
		"Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel."
	FORCE)
endif()

project(desyre VERSION 0.0.1 LANGUAGES CXX)
# We allow to find yacma 
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake" "${CMAKE_CURRENT_SOURCE_DIR}/cmake/yacma")

# Print some initial information on the build
message(STATUS "System name: ${CMAKE_SYSTEM_NAME}")
message(STATUS "System processor: ${CMAKE_SYSTEM_PROCESSOR}")
message(STATUS "desyre version: ${desyre_VERSION}")

# Run the YACMA compiler setup.
include(YACMACompilerLinkerSettings)

# Options
option(desyre_BUILD_TESTS "Build unit tests." OFF)
option(desyre_BUILD_STATIC_LIBRARY "Build desyre as a static library, instead of dynamic." OFF)
option(desyre_SETUP_DOCS "Setup the files for building the docs." ON)


# List of source files.
set(desyre_SRC_FILES
    "${CMAKE_CURRENT_SOURCE_DIR}/src/expression.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/kernels.cpp"
)
# desyre library ###################################################################################################################3
if(desyre_BUILD_STATIC_LIBRARY)
    # Setup of the desyre static library.
    message(STATUS "desyre will be built as a static library.")
    add_library(desyre STATIC "${desyre_SRC_FILES}")
else()
    # Setup of the desyre shared library.
    add_library(desyre SHARED "${desyre_SRC_FILES}")
    set_property(TARGET desyre PROPERTY VERSION "0.0")
    set_property(TARGET desyre PROPERTY SOVERSION 0)
    set_target_properties(desyre PROPERTIES CXX_VISIBILITY_PRESET hidden)
    set_target_properties(desyre PROPERTIES VISIBILITY_INLINES_HIDDEN TRUE)
endif()

# Setup common to both the shared and static variants.
target_compile_options(desyre PRIVATE
    "$<$<CONFIG:Debug>:${desyre_CXX_FLAGS_DEBUG}>"
    "$<$<CONFIG:Release>:${desyre_CXX_FLAGS_RELEASE}>"
    "$<$<CONFIG:RelWithDebInfo>:${desyre_CXX_FLAGS_RELEASE}>"
    "$<$<CONFIG:MinSizeRel>:${desyre_CXX_FLAGS_RELEASE}>"
)

# Ensure that C++17 is employed when both compiling and consuming desyre.
target_compile_features(desyre PUBLIC cxx_std_17)
# Enforce vanilla C++17 when compiling desyre.
set_property(TARGET desyre PROPERTY CXX_EXTENSIONS NO)

target_include_directories(desyre PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/include>
    $<INSTALL_INTERFACE:include>)

# Mandatory dependency on fmt.
find_package(fmt)
target_link_libraries(desyre PRIVATE fmt::fmt)

# Mandatory dependency on Boost.
find_package(Boost 1.60 REQUIRED serialization)
target_link_libraries(desyre PUBLIC Boost::boost Boost::serialization)
# NOTE: quench warnings from Boost when building the library.
target_compile_definitions(desyre PRIVATE BOOST_ALLOW_DEPRECATED_HEADERS)


# main.cpp ###################################################################################################################3
add_executable(main main.cpp)
target_link_libraries(main fmt::fmt)

set_property(TARGET main PROPERTY CXX_STANDARD 17)
set_property(TARGET main PROPERTY CXX_STANDARD_REQUIRED YES)
set_property(TARGET main PROPERTY CXX_EXTENSIONS NO)

# koza_quintic_simple.cpp ###################################################################################################################3
add_executable(koza koza_quintic_simple.cpp)
target_link_libraries(koza fmt::fmt)

set_property(TARGET koza PROPERTY CXX_STANDARD 17)
set_property(TARGET koza PROPERTY CXX_STANDARD_REQUIRED YES)
set_property(TARGET koza PROPERTY CXX_EXTENSIONS NO)

# auto_diff_2nd_order.cpp ###################################################################################################################3
add_executable(autodiff auto_diff_2nd_order.cpp)
target_link_libraries(autodiff fmt::fmt)

set_property(TARGET autodiff PROPERTY CXX_STANDARD 17)
set_property(TARGET autodiff PROPERTY CXX_STANDARD_REQUIRED YES)
set_property(TARGET autodiff PROPERTY CXX_EXTENSIONS NO)

# Configuration files #############################################################################################################3
# Configure config.hpp.
configure_file("${CMAKE_CURRENT_SOURCE_DIR}/config.hpp.in" "${CMAKE_CURRENT_BINARY_DIR}/include/desyre/config.hpp" @ONLY)

if(desyre_SETUP_DOCS)
    # Configure the doc files.
    configure_file("${CMAKE_CURRENT_SOURCE_DIR}/doc/conf.py.in" "${CMAKE_CURRENT_SOURCE_DIR}/doc/conf.py" @ONLY)
endif()

# Installation of the header files.
install(DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/include/desyre" DESTINATION include)
install(FILES "${CMAKE_CURRENT_BINARY_DIR}/include/desyre/config.hpp" DESTINATION include/desyre)

# Installation of the library.
install(TARGETS desyre
    EXPORT desyre_export
    LIBRARY DESTINATION "${desyre_INSTALL_LIBDIR}"
    ARCHIVE DESTINATION "${desyre}"
    RUNTIME DESTINATION bin
)

# Setup of the CMake config file.
configure_file("${CMAKE_CURRENT_SOURCE_DIR}/desyre-config.cmake.in" "${CMAKE_CURRENT_BINARY_DIR}/desyre-config.cmake" @ONLY)
install(FILES "${CMAKE_CURRENT_BINARY_DIR}/desyre-config.cmake"
    DESTINATION "${desyre_INSTALL_LIBDIR}/cmake/desyre")
install(EXPORT desyre_export NAMESPACE desyre:: DESTINATION "${desyre_INSTALL_LIBDIR}/cmake/desyre")
# Take care of versioning.
include(CMakePackageConfigHelpers)
# NOTE: SameMinorVersion available only
# since CMake 3.11.
if(${CMAKE_VERSION} VERSION_LESS "3.11.0")
    write_basic_package_version_file("${CMAKE_CURRENT_BINARY_DIR}/desyre-config-version.cmake" COMPATIBILITY SameMajorVersion)
else()
    write_basic_package_version_file("${CMAKE_CURRENT_BINARY_DIR}/desyre-config-version.cmake" COMPATIBILITY SameMinorVersion)
endif()
install(FILES "${CMAKE_CURRENT_BINARY_DIR}/desyre-config-version.cmake" DESTINATION "${desyre_INSTALL_LIBDIR}/cmake/desyre")