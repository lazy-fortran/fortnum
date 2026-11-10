cmake_minimum_required(VERSION 3.20)
project(fortnum_host_only_cxx_link LANGUAGES Fortran CXX)

find_package(OpenMP REQUIRED COMPONENTS Fortran)
include("${FORTNUM_SOURCE_DIR}/cmake/FortnumGpuBackend.cmake")

add_executable(host_only_cxx_link kernel.f90)
target_compile_options(host_only_cxx_link PRIVATE ${OpenMP_Fortran_FLAGS})
target_link_libraries(host_only_cxx_link PRIVATE
    fortnum_gpu_backend OpenMP::OpenMP_Fortran)
# Mixed-language consumers such as KiLCA link through the C++ driver.
set_target_properties(host_only_cxx_link PROPERTIES LINKER_LANGUAGE CXX)
