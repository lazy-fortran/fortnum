cmake_minimum_required(VERSION 3.22)

project(fortnum LANGUAGES Fortran)

# Detect standalone vs subproject build so testing defaults sanely.
set(FORTNUM_IS_STANDALONE OFF)
if(CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME)
    set(FORTNUM_IS_STANDALONE ON)
endif()

option(FORTNUM_BUILD_TESTING "Build tests for fortnum" ${FORTNUM_IS_STANDALONE})

# Enzyme autodiff integration lands in #38; keep the switches inert for now.
option(FORTNUM_ENABLE_ENZYME "Enable Enzyme-based reverse-mode autodiff" OFF)
option(FORTNUM_ENZYME_REQUIRED "Fail configure if Enzyme is unavailable" OFF)

# Shared library inclusion needs PIC-clean objects.
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# Pin a modern Fortran standard for free-form, full-feature builds.
set(CMAKE_Fortran_STANDARD 2018)
set(CMAKE_Fortran_STANDARD_REQUIRED ON)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)
set(CMAKE_Fortran_MODULE_DIRECTORY ${PROJECT_BINARY_DIR}/include)

enable_testing()

add_subdirectory(src)

if(FORTNUM_BUILD_TESTING)
    add_subdirectory(test)
endif()

option(FORTNUM_BUILD_EXAMPLES "Build downstream API examples" ${FORTNUM_IS_STANDALONE})
if(FORTNUM_BUILD_EXAMPLES)
    add_subdirectory(examples)
endif()
