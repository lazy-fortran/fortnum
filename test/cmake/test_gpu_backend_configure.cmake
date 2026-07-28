if(NOT DEFINED FORTNUM_SOURCE_DIR)
    message(FATAL_ERROR "FORTNUM_SOURCE_DIR is required")
endif()
if(NOT DEFINED TEST_BINARY_ROOT)
    message(FATAL_ERROR "TEST_BINARY_ROOT is required")
endif()

function(check_configure label expected_result expected_text)
    set(build_dir "${TEST_BINARY_ROOT}/${label}")
    execute_process(
        COMMAND "${CMAKE_COMMAND}"
            -S "${FORTNUM_SOURCE_DIR}"
            -B "${build_dir}"
            -DFORTNUM_BUILD_TESTING=OFF
            -DFORTNUM_BUILD_EXAMPLES=OFF
            ${ARGN}
        RESULT_VARIABLE result
        OUTPUT_VARIABLE output
        ERROR_VARIABLE error)
    set(combined "${output}\n${error}")

    if(expected_result STREQUAL "success")
        if(NOT result EQUAL 0)
            message(FATAL_ERROR
                "${label}: configure unexpectedly failed:\n${combined}")
        endif()
    else()
        if(result EQUAL 0)
            message(FATAL_ERROR "${label}: configure unexpectedly succeeded")
        endif()
    endif()
    if(NOT combined MATCHES "${expected_text}")
        message(FATAL_ERROR
            "${label}: expected '${expected_text}' in output:\n${combined}")
    endif()
endfunction()

check_configure(
    none_accepted success "Build files have been written"
    -DFORTNUM_GPU_BACKEND=NONE)
check_configure(
    invalid_rejected failure "must be NONE, OPENACC, or OPENMP"
    -DFORTNUM_GPU_BACKEND=CUDA)
check_configure(
    host_openmp_rejected failure "requires explicit.*offload target"
    -DFORTNUM_GPU_BACKEND=OPENMP
    -DFORTNUM_OPENMP_TARGET_FLAGS=)
check_configure(
    unvalidated_openmp_compiler_rejected failure
    "supports only the validated NVHPC nvfortran compiler"
    -DFORTNUM_GPU_BACKEND=OPENMP
    -DFORTNUM_OPENMP_TARGET_FLAGS=-fopenmp)
