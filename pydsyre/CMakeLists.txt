# Configure the version file.
configure_file("${CMAKE_CURRENT_SOURCE_DIR}/_version.py.in" "${CMAKE_CURRENT_SOURCE_DIR}/_version.py" @ONLY)   

# The list of pydsyre's Python files.
set(PYDSYRE_PYTHON_FILES __init__.py _version.py)

# Copy the python files in the current binary dir,
# so that we can import pydesyre from the build dir.
# NOTE: importing from the build dir will work
# only on single-configuration generators.
foreach(PYDSYRE_PYTHON_FILE ${PYDSYRE_PYTHON_FILES})
    configure_file("${CMAKE_CURRENT_SOURCE_DIR}/${PYDSYRE_PYTHON_FILE}"
        "${CMAKE_CURRENT_BINARY_DIR}/${PYDSYRE_PYTHON_FILE}" COPYONLY)
endforeach()

# Core module.
Python3_add_library(core MODULE WITH_SOABI
    core.cpp
    docstrings.cpp
)

target_link_libraries(core PRIVATE dsyre::dsyre)
target_link_libraries(core PRIVATE "${pybind11_LIBRARIES}")

target_include_directories(core SYSTEM PRIVATE "${pybind11_INCLUDE_DIR}" "${Python3_INCLUDE_DIRS}")
target_compile_definitions(core PRIVATE "${pybind11_DEFINITIONS}")
target_compile_options(core PRIVATE
    "$<$<CONFIG:Debug>:${DSYRE_CXX_FLAGS_DEBUG}>"
    "$<$<CONFIG:Release>:${DSYRE_CXX_FLAGS_RELEASE}>"
    "$<$<CONFIG:RelWithDebInfo>:${DSYRE_CXX_FLAGS_RELEASE}>"
    "$<$<CONFIG:MinSizeRel>:${DSYRE_CXX_FLAGS_RELEASE}>"
)
set_target_properties(core PROPERTIES CXX_VISIBILITY_PRESET hidden)
set_target_properties(core PROPERTIES VISIBILITY_INLINES_HIDDEN TRUE)
target_compile_features(core PRIVATE cxx_std_17)
set_property(TARGET core PROPERTY CXX_EXTENSIONS NO)
