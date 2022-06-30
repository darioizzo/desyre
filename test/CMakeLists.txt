add_library(dsyre_test STATIC catch_main.cpp)
target_compile_options(dsyre_test PRIVATE
  "$<$<CONFIG:Debug>:${DSYRE_CXX_FLAGS_DEBUG}>"
  "$<$<CONFIG:Release>:${DSYRE_CXX_FLAGS_RELEASE}>"
  "$<$<CONFIG:RelWithDebInfo>:${DSYRE_CXX_FLAGS_RELEASE}>"
  "$<$<CONFIG:MinSizeRel>:${DSYRE_CXX_FLAGS_RELEASE}>"
)

# Setup the C++ standard.
target_compile_features(dsyre_test PRIVATE cxx_std_17)
set_property(TARGET dsyre_test PROPERTY CXX_EXTENSIONS NO)

function(ADD_DSYRE_TESTCASE arg1)
  add_executable(${arg1} ${arg1}.cpp)
  target_link_libraries(${arg1} PRIVATE dsyre_test dsyre fmt::fmt)
  target_compile_options(${arg1} PRIVATE
    "$<$<CONFIG:Debug>:${DSYRE_CXX_FLAGS_DEBUG}>"
    "$<$<CONFIG:Release>:${DSYRE_CXX_FLAGS_RELEASE}>"
    "$<$<CONFIG:RelWithDebInfo>:${DSYRE_CXX_FLAGS_RELEASE}>"
    "$<$<CONFIG:MinSizeRel>:${DSYRE_CXX_FLAGS_RELEASE}>"
  )
  # Setup the C++ standard.
  target_compile_features(${arg1} PRIVATE cxx_std_17)
  set_property(TARGET ${arg1} PROPERTY CXX_EXTENSIONS NO)
  add_test(${arg1} ${arg1})
endfunction()

ADD_DSYRE_TESTCASE(expression)
