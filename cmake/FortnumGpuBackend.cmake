include_guard(GLOBAL)

set(FORTNUM_GPU_BACKEND "NONE" CACHE STRING
    "GPU backend for annotated analytical leaves: NONE, OPENACC, or OPENMP")
set_property(CACHE FORTNUM_GPU_BACKEND PROPERTY STRINGS NONE OPENACC OPENMP)
set(FORTNUM_OPENMP_TARGET_FLAGS "" CACHE STRING
    "Compiler and linker flags naming the OpenMP offload target")
string(TOUPPER "${FORTNUM_GPU_BACKEND}" FORTNUM_GPU_BACKEND)
set(FORTNUM_GPU_BACKEND "${FORTNUM_GPU_BACKEND}" CACHE STRING
    "GPU backend for annotated analytical leaves: NONE, OPENACC, or OPENMP"
    FORCE)

set(_fortnum_gpu_backends NONE OPENACC OPENMP)
if(NOT FORTNUM_GPU_BACKEND IN_LIST _fortnum_gpu_backends)
    message(FATAL_ERROR
        "FORTNUM_GPU_BACKEND must be NONE, OPENACC, or OPENMP; got "
        "'${FORTNUM_GPU_BACKEND}'")
endif()

add_library(fortnum_gpu_backend INTERFACE)

if(FORTNUM_GPU_BACKEND STREQUAL "OPENACC")
    find_package(OpenACC QUIET COMPONENTS Fortran)
    if(NOT TARGET OpenACC::OpenACC_Fortran)
        message(FATAL_ERROR
            "FORTNUM_GPU_BACKEND=OPENACC requested, but this released "
            "Fortran compiler has no CMake-detectable OpenACC support")
    endif()
    target_link_libraries(fortnum_gpu_backend INTERFACE
        OpenACC::OpenACC_Fortran)
elseif(FORTNUM_GPU_BACKEND STREQUAL "OPENMP")
    if(NOT FORTNUM_OPENMP_TARGET_FLAGS)
        message(FATAL_ERROR
            "FORTNUM_GPU_BACKEND=OPENMP requires explicit "
            "FORTNUM_OPENMP_TARGET_FLAGS for a real offload target; "
            "host-only OpenMP is not accepted")
    endif()
    find_package(OpenMP QUIET COMPONENTS Fortran)
    if(NOT TARGET OpenMP::OpenMP_Fortran)
        message(FATAL_ERROR
            "FORTNUM_GPU_BACKEND=OPENMP requested, but this released "
            "Fortran compiler has no CMake-detectable OpenMP support")
    endif()
    separate_arguments(_fortnum_openmp_target_flags NATIVE_COMMAND
        "${FORTNUM_OPENMP_TARGET_FLAGS}")
    target_compile_options(fortnum_gpu_backend INTERFACE
        $<$<COMPILE_LANGUAGE:Fortran>:${_fortnum_openmp_target_flags}>)
    target_link_options(fortnum_gpu_backend INTERFACE
        ${_fortnum_openmp_target_flags})
    target_link_libraries(fortnum_gpu_backend INTERFACE
        OpenMP::OpenMP_Fortran)
endif()
