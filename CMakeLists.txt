cmake_minimum_required(VERSION 3.22)

project(fortnum LANGUAGES Fortran C)

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

# Internal procedures passed as actual arguments need closures. gfortran's
# default stack trampolines mark the object's .note.GNU-stack executable, which
# makes the linker emit an executable-stack warning. Heap-allocated trampolines
# (GCC 14+) keep the closures working while leaving the stack non-executable, so
# the warning disappears without an executable stack and without breaking the
# procedure-argument paths in ode/roots/integrate/oracle.
include(CheckFortranCompilerFlag)
check_fortran_compiler_flag("-ftrampoline-impl=heap" FORTNUM_HAVE_HEAP_TRAMPOLINES)
if(FORTNUM_HAVE_HEAP_TRAMPOLINES)
    add_compile_options($<$<COMPILE_LANGUAGE:Fortran>:-ftrampoline-impl=heap>)
endif()

# Pin a modern Fortran standard for free-form, full-feature builds.
set(CMAKE_Fortran_STANDARD 2018)
set(CMAKE_Fortran_STANDARD_REQUIRED ON)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)
set(CMAKE_Fortran_MODULE_DIRECTORY ${PROJECT_BINARY_DIR}/include)

enable_testing()

# The complex analytic-zero finder (fortnum_roots_complex) solves a small
# generalized eigenproblem with LAPACK ZGGEV, so the library links LAPACK (and
# its BLAS dependency). This is the only external numerical dependency.
find_package(LAPACK REQUIRED)

add_subdirectory(src)

if(FORTNUM_BUILD_TESTING)
    add_subdirectory(test)
endif()

option(FORTNUM_BUILD_EXAMPLES "Build downstream API examples" ${FORTNUM_IS_STANDALONE})
if(FORTNUM_BUILD_EXAMPLES)
    add_subdirectory(examples)
endif()
