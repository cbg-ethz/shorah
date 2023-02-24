cmake_minimum_required(VERSION 3.4...3.18)
project(libshorah)


list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake")

# Find dependencies:

# python installation path sheningans
find_package(Python COMPONENTS Interpreter Development)
# NOTE Python_SITEARCH is absolute path, but CMAKE_PREFIX_PATH only works with relative paths
execute_process(
  COMMAND "${Python_EXECUTABLE}" -c "from distutils import sysconfig as sc;print(sc.get_python_lib(prefix='', plat_specific=True))"
  OUTPUT_VARIABLE PYTHON_SITE
  OUTPUT_STRIP_TRAILING_WHITESPACE)
message(STATUS "Python module destination: ${PYTHON_SITE}")

# pybind11
find_package(pybind11 REQUIRED)

# Find Boost
find_package(Boost REQUIRED)

# Tools used to find packages (e.g.: HTSlib)
find_package(HTSlib REQUIRED htslib)

# autodetect popcnt support
include(popcnt)


# Build settings:
# NOTE this file MUST be kept in sync with the build.py script

file(GLOB libshorah_src ${PROJECT_SOURCE_DIR}/lib/src/*.cpp) # build.py's sources
pybind11_add_module(libshorah ${libshorah_src})
install(TARGETS libshorah LIBRARY DESTINATION ${PYTHON_SITE})

include_directories("${PROJECT_SOURCE_DIR}/lib/include" # build.py's  include_dirs
	${pybind11_INCLUDE_DIR} ${Boost_INCLUDE_DIRS} ${HTSlib_INCLUDE_DIRS} )

target_link_libraries(libshorah PRIVATE ${HTSlib_LIBRARIES}) # build.py's libraries
	# NOTE currently boost is used includes only

# NOTE CMake has own mechanism for language standard
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11" ) # build.py's extra_compile_args
set_property(TARGET libshorah PROPERTY CXX_STANDARD 11)

# NOTE for production, we keep the HAVE_POPCNT
#target_compile_definitions(libshorah PRIVATE -DHAVE_POPCNT) # build.py's undef_macros
target_compile_definitions(libshorah PRIVATE HAVE_POPCNT=${HAVE_POPCNT} )
