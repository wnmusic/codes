cmake_minimum_required(VERSION 3.0)   
project(libcode)

#set ( CMAKE_MODULE_PATH
#  ${PROJECT_SOURCE_DIR}/cmake
#)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Debug)
endif()


SET(CMAKE_C_FLAGS  "-fPIC")
set(CMAKE_C_FLAGS_RELEASE "-O3 -mavx")
add_library(libcode src/conv.c src/gf2.c src/rs.c src/gf2m.c)

target_include_directories(libcode PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/inc)
target_link_libraries(libcode m)

if (CMAKE_BUILD_TYPE STREQUAL "Debug")
   add_compile_definitions("_DEBUG=1")    
endif()

add_subdirectory(test)
