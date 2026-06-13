# fortnum core library. Downstream agents append one module source per line
# under the "module sources" block below, kept alphabetical and one-per-line
# to minimise merge conflicts.
add_library(fortnum STATIC)

target_sources(fortnum PRIVATE
    # module sources (alphabetical, one per line; append here)
    fft/fortnum_fft.f90
    fortnum_kinds.f90
    fortnum_status.f90
    fortnum_version.f90
    quadrature/fortnum_integrate_gk.f90
    quadrature/fortnum_quadrature.f90
    special/fortnum_special.f90
    special/fortnum_special_bessel.f90
    special/fortnum_special_dawson.f90
    special/fortnum_special_gamma.f90
    testing/fortnum_oracle.f90
)

# Consumers find .mod files in the module output directory.
target_include_directories(fortnum PUBLIC ${CMAKE_Fortran_MODULE_DIRECTORY})
