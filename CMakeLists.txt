cmake_minimum_required(VERSION 3.17.2...3.29)
project(
  ${SKBUILD_PROJECT_NAME}
  LANGUAGES C Fortran
  VERSION ${SKBUILD_PROJECT_VERSION})

set(ModuleName "_thermal_comfort")
set(CMAKE_Fortran_FLAGS
    "${CMAKE_Fortran_FLAGS} -O3 -Wall -march=native -ftree-vectorize -funroll-loops -cpp"
)

find_package(
  Python
  COMPONENTS Interpreter Development.Module NumPy ${SKBUILD_SABI_COMPONENT}
  REQUIRED)

# F2PY headers
execute_process(
  COMMAND "${PYTHON_EXECUTABLE}" -c
          "import numpy.f2py; print(numpy.f2py.get_include())"
  OUTPUT_VARIABLE F2PY_INCLUDE_DIR
  OUTPUT_STRIP_TRAILING_WHITESPACE)

add_library(fortranobject OBJECT "${F2PY_INCLUDE_DIR}/fortranobject.c")
target_link_libraries(fortranobject PUBLIC Python::NumPy)
target_include_directories(fortranobject PUBLIC "${F2PY_INCLUDE_DIR}")
set_property(TARGET fortranobject PROPERTY POSITION_INDEPENDENT_CODE ON)

add_custom_command(
  OUTPUT "${ModuleName}module.c" "${ModuleName}-f2pywrappers2.f90"
  DEPENDS "src/${ModuleName}.f90"
  VERBATIM
  COMMAND
    "${Python_EXECUTABLE}" -m numpy.f2py
    "${CMAKE_CURRENT_SOURCE_DIR}/src/${ModuleName}.f90" -m "${ModuleName}"
    --lower)

if(NOT "${SKBUILD_SABI_COMPONENT}" STREQUAL "")
  python_add_library(
    ${ModuleName}
    MODULE
    "${CMAKE_CURRENT_BINARY_DIR}/${ModuleName}module.c"
    "${CMAKE_CURRENT_BINARY_DIR}/${ModuleName}-f2pywrappers2.f90"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/${ModuleName}.f90"
    WITH_SOABI
    USE_SABI
    3.7)
else()
  python_add_library(
    ${ModuleName} MODULE "${CMAKE_CURRENT_BINARY_DIR}/${ModuleName}module.c"
    "${CMAKE_CURRENT_BINARY_DIR}/${ModuleName}-f2pywrappers2.f90"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/${ModuleName}.f90" WITH_SOABI)
endif()

target_link_libraries("${ModuleName}" PRIVATE fortranobject)

install(TARGETS "${ModuleName}" DESTINATION "thermal_comfort")
