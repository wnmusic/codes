cmake_policy(SET CMP0078 NEW)
cmake_policy(SET CMP0086 NEW)

find_package(SWIG REQUIRED)
find_package(Python3 EXACT 3.11.6 COMPONENTS Interpreter Development NumPy)

INCLUDE(${SWIG_USE_FILE})

swig_add_library(pycode LANGUAGE python SOURCES pycode.i)
set_property(TARGET pycode PROPERTY SWIG_INCLUDE_DIRECTORIES ${CMAKE_CURRENT_SOURCE_DIR}/../inc)
set_property(TARGET pycode PROPERTY SWIG_GENERATED_INCLUDE_DIRECTORIES ${Python3_INCLUDE_DIRS} ${Python3_NumPy_INCLUDE_DIRS})
target_link_libraries(pycode libcode ${Python3_LIBRARIES})

