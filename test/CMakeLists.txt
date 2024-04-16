cmake_policy(SET CMP0078 NEW)
cmake_policy(SET CMP0086 NEW)

set (Python_FIND_VIRTUALENV ONLY)
find_package(Python COMPONENTS Interpreter Development NumPy)

find_package(SWIG REQUIRED)
INCLUDE(${SWIG_USE_FILE})

swig_add_library(pycode LANGUAGE python SOURCES pycode.i)
set_property(TARGET pycode PROPERTY SWIG_INCLUDE_DIRECTORIES ${CMAKE_CURRENT_SOURCE_DIR}/../inc)
set_property(TARGET pycode PROPERTY SWIG_GENERATED_INCLUDE_DIRECTORIES ${Python_INCLUDE_DIRS} ${Python_NumPy_INCLUDE_DIRS})
target_link_libraries(pycode libcode ${Python_LIBRARIES})

