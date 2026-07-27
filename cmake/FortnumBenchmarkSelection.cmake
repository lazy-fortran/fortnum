function(fortnum_read_benchmark_selection record output_variable)
    if(NOT EXISTS "${record}")
        message(FATAL_ERROR "fortnum benchmark record not found: ${record}")
    endif()

    file(READ "${record}" record_json)
    string(JSON selected ERROR_VARIABLE json_error
        GET "${record_json}" selection selected)
    if(NOT json_error STREQUAL "NOTFOUND")
        message(FATAL_ERROR
            "invalid fortnum benchmark selection in ${record}: ${json_error}")
    endif()
    if(selected STREQUAL "")
        message(FATAL_ERROR "empty fortnum benchmark selection in ${record}")
    endif()

    set(${output_variable} "${selected}" PARENT_SCOPE)
endfunction()
