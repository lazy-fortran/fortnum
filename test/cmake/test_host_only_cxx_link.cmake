# With FORTNUM_GPU_BACKEND=NONE, a consumer that links through the C++ driver
# must not build device code. GNU emits it for every "!$omp declare target"
# leaf whenever an offload compiler is installed, and only the link flag
# -foffload=disable stops mkoffload (gfortran 12 nvptx: unresolved sincos).
if(NOT DEFINED FORTNUM_SOURCE_DIR OR NOT DEFINED TEST_BINARY_ROOT)
    message(FATAL_ERROR "FORTNUM_SOURCE_DIR and TEST_BINARY_ROOT are required")
endif()
set(build_dir "${TEST_BINARY_ROOT}/host-only-cxx-link")
file(REMOVE_RECURSE "${build_dir}")
foreach(step configure build run)
    if(step STREQUAL "configure")
        set(command "${CMAKE_COMMAND}"
            -S "${FORTNUM_SOURCE_DIR}/test/cmake/host_only_cxx_link"
            -B "${build_dir}" "-DFORTNUM_SOURCE_DIR=${FORTNUM_SOURCE_DIR}"
            -DFORTNUM_GPU_BACKEND=NONE)
    elseif(step STREQUAL "build")
        set(command "${CMAKE_COMMAND}" --build "${build_dir}")
    else()
        set(command "${build_dir}/host_only_cxx_link")
    endif()
    execute_process(COMMAND ${command}
        RESULT_VARIABLE result OUTPUT_VARIABLE output ERROR_VARIABLE error)
    if(NOT result EQUAL 0)
        message(FATAL_ERROR "${step} failed:\n${output}\n${error}")
    endif()
endforeach()
