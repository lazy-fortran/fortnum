include("${FORTNUM_SOURCE_DIR}/cmake/FortnumBenchmarkSelection.cmake")

set(record "${CMAKE_CURRENT_BINARY_DIR}/benchmark-selection-fixture.json")
file(WRITE "${record}"
    "{\"selection\":{\"selected\":\"hybrid\"}}")
fortnum_read_benchmark_selection("${record}" selected)
if(NOT selected STREQUAL "hybrid")
    message(FATAL_ERROR "expected hybrid selection, got '${selected}'")
endif()
