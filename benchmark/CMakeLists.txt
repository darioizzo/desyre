# NOTE: don't find a specific version as we already checked
# outside that the Boost version is appropriate.
find_package(Boost REQUIRED COMPONENTS program_options)
if(NOT TARGET Boost::program_options)
    message(STATUS "The 'Boost::program_options' imported target is missing, creating it.")
    add_library(Boost::program_options UNKNOWN IMPORTED)
    set_target_properties(Boost::program_options PROPERTIES
        INTERFACE_INCLUDE_DIRECTORIES "${Boost_INCLUDE_DIRS}"
    )
    set_target_properties(Boost::program_options PROPERTIES
        IMPORTED_LINK_INTERFACE_LANGUAGES "CXX"
        IMPORTED_LOCATION "${Boost_PROGRAM_OPTIONS_LIBRARY}"
    )
endif()

function(ADD_DSYRE_BENCHMARK arg1)
  add_executable(${arg1} ${arg1}.cpp)
  # NOTE: fmt was already located in the main CMakeLists.
  target_link_libraries(${arg1} PRIVATE dsyre Boost::boost Boost::program_options fmt::fmt Pagmo::pagmo)
  target_compile_definitions(${arg1} PRIVATE BOOST_ALLOW_DEPRECATED_HEADERS)
  target_compile_options(${arg1} PRIVATE
    "$<$<CONFIG:Debug>:${DSYRE_CXX_FLAGS_DEBUG}>"
    "$<$<CONFIG:Release>:${DSYRE_CXX_FLAGS_RELEASE}>"
    "$<$<CONFIG:RelWithDebInfo>:${DSYRE_CXX_FLAGS_RELEASE}>"
    "$<$<CONFIG:MinSizeRel>:${DSYRE_CXX_FLAGS_RELEASE}>"
  )
  # Setup the C++ standard.
  target_compile_features(${arg1} PRIVATE cxx_std_17)
  set_property(TARGET ${arg1} PROPERTY CXX_EXTENSIONS NO)
endfunction()

ADD_DSYRE_BENCHMARK(phenotype)
ADD_DSYRE_BENCHMARK(dphenotype)
ADD_DSYRE_BENCHMARK(ddphenotype)
ADD_DSYRE_BENCHMARK(sr_problem_fitness)