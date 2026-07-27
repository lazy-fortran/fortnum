cmake_minimum_required(VERSION 3.22)

project(fortnum LANGUAGES Fortran C)

# Detect standalone vs subproject build so testing defaults sanely.
set(FORTNUM_IS_STANDALONE OFF)
if(CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME)
    set(FORTNUM_IS_STANDALONE ON)
endif()

option(FORTNUM_BUILD_TESTING "Build tests for fortnum" ${FORTNUM_IS_STANDALONE})
option(FORTNUM_CHECK_GENERATED
    "Regenerate and byte-compare committed fortsym kernels" OFF)

option(FORTNUM_USE_BENCHMARK_SELECTION
    "Select derivative implementations from committed benchmark records" ON)
set(FORTNUM_DAWSON_BENCHMARK_RECORD
    "${CMAKE_CURRENT_SOURCE_DIR}/benchmark/reference/ryzen9_5950x_dawson_outer.json"
    CACHE FILEPATH "Benchmark record used for the Dawson outer JVP")
set(FORTNUM_DERIVATIVE_BENCHMARK_RECORDS
    "${CMAKE_CURRENT_SOURCE_DIR}/benchmark/reference/ryzen9_5950x_dawson_outer.json"
    "${CMAKE_CURRENT_SOURCE_DIR}/benchmark/reference/ryzen9_5950x_linear_solve_jvp_reuse.json"
    "${CMAKE_CURRENT_SOURCE_DIR}/benchmark/reference/ryzen9_5950x_linear_solve_vjp_reuse.json"
    "${CMAKE_CURRENT_SOURCE_DIR}/benchmark/reference/ryzen9_5950x_multiroot_preconditioner.json")

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
include(FortnumBenchmarkSelection)
if(FORTNUM_USE_BENCHMARK_SELECTION)
    fortnum_read_benchmark_selection(
        "${FORTNUM_DAWSON_BENCHMARK_RECORD}"
        FORTNUM_DAWSON_OUTER_JVP_CANDIDATE)
else()
    set(FORTNUM_DAWSON_OUTER_JVP_CANDIDATE "analytical")
endif()
if(NOT FORTNUM_DAWSON_OUTER_JVP_CANDIDATE STREQUAL "analytical")
    message(FATAL_ERROR
        "Dawson outer JVP candidate '${FORTNUM_DAWSON_OUTER_JVP_CANDIDATE}' "
        "is not available in the production library")
endif()

set(FORTNUM_REGISTRY_OPERATORS)
set(FORTNUM_REGISTRY_PRODUCTS)
set(FORTNUM_REGISTRY_WORKLOADS)
set(FORTNUM_REGISTRY_CANDIDATES)
foreach(record IN LISTS FORTNUM_DERIVATIVE_BENCHMARK_RECORDS)
    fortnum_read_benchmark_entry("${record}" operator product workload candidate)
    list(APPEND FORTNUM_REGISTRY_OPERATORS "${operator}")
    list(APPEND FORTNUM_REGISTRY_PRODUCTS "${product}")
    list(APPEND FORTNUM_REGISTRY_WORKLOADS "${workload}")
    list(APPEND FORTNUM_REGISTRY_CANDIDATES "${candidate}")
endforeach()
list(LENGTH FORTNUM_REGISTRY_OPERATORS FORTNUM_REGISTRY_SIZE)
fortnum_format_fortran_strings(
    FORTNUM_REGISTRY_OPERATORS FORTNUM_REGISTRY_OPERATORS_FORTRAN)
fortnum_format_fortran_strings(
    FORTNUM_REGISTRY_PRODUCTS FORTNUM_REGISTRY_PRODUCTS_FORTRAN)
fortnum_format_fortran_strings(
    FORTNUM_REGISTRY_WORKLOADS FORTNUM_REGISTRY_WORKLOADS_FORTRAN)
fortnum_format_fortran_strings(
    FORTNUM_REGISTRY_CANDIDATES FORTNUM_REGISTRY_CANDIDATES_FORTRAN)

set(FORTNUM_BUILD_SELECTION_SOURCE
    "${CMAKE_CURRENT_BINARY_DIR}/generated/fortnum_build_selection.f90")
configure_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/cmake/fortnum_build_selection.f90.in"
    "${FORTNUM_BUILD_SELECTION_SOURCE}"
    @ONLY)

# Enzyme autodiff integration lands in #38; keep the switches inert for now.
option(FORTNUM_ENABLE_ENZYME "Enable Enzyme-based autodiff" OFF)
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
else()
    # Without heap trampolines (gfortran < 14) the closures still compile via
    # stack trampolines. A consumer that turns trampolines into a hard error
    # (libneo builds fortnum with -Werror=trampolines) would then break the
    # build, so keep the warning non-fatal for fortnum's own sources.
    check_fortran_compiler_flag("-Wno-error=trampolines" FORTNUM_HAVE_NO_ERROR_TRAMPOLINES)
    if(FORTNUM_HAVE_NO_ERROR_TRAMPOLINES)
        add_compile_options($<$<COMPILE_LANGUAGE:Fortran>:-Wno-error=trampolines>)
    endif()
endif()

# Pin a modern Fortran standard for free-form, full-feature builds.
set(CMAKE_Fortran_STANDARD 2018)
set(CMAKE_Fortran_STANDARD_REQUIRED ON)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)
set(CMAKE_Fortran_MODULE_DIRECTORY ${PROJECT_BINARY_DIR}/include)

enable_testing()

if(FORTNUM_CHECK_GENERATED)
    find_program(FORTNUM_FO_EXECUTABLE fo REQUIRED)
    if(NOT EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/../lazy-fortran/fortsym/fpm.toml")
        message(FATAL_ERROR
            "FORTNUM_CHECK_GENERATED requires ../lazy-fortran/fortsym")
    endif()
    add_test(
        NAME generated_sources_current
        COMMAND "${CMAKE_CURRENT_SOURCE_DIR}/tools/codegen/check_generated.sh")
    set_tests_properties(generated_sources_current PROPERTIES
        WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/tools/codegen"
        LABELS codegen)
endif()

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
