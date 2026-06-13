# Downstream example programs (issue #42). Each example is a self-contained
# program that links fortnum and exercises the optimizer-facing AD API (#41).
# Registered as ctests so they run under the normal test suite.
#
# Wire from the top-level CMakeLists.txt with:
#   option(FORTNUM_BUILD_EXAMPLES "Build downstream examples" ${FORTNUM_IS_STANDALONE})
#   if(FORTNUM_BUILD_EXAMPLES)
#       add_subdirectory(examples)
#   endif()

add_executable(example_scalar_objective
    downstream/scalar_objective/scalar_objective.f90)
target_link_libraries(example_scalar_objective PRIVATE fortnum)
add_test(NAME example_scalar_objective COMMAND example_scalar_objective)

add_executable(example_vector_residual
    downstream/vector_residual/vector_residual.f90)
target_link_libraries(example_vector_residual PRIVATE fortnum)
add_test(NAME example_vector_residual COMMAND example_vector_residual)
