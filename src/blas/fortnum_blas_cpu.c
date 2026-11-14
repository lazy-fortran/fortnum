/* The GEMM register tiles for the running CPU, as function pointers, so that
   no kernel of another architecture is referenced: the tile kernel and the
   same tile adding to C in place. FORTNUM_BLAS_ISA (avx2) forces the AVX2
   tile on a CPU that also has AVX-512. */
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

typedef void tile_t(int64_t k, const double *a, const double *b, double *acc);
typedef void tile_c_t(int64_t k, const double *a, const double *b, double *c, int64_t ldc);

#if defined(__x86_64__) && !defined(_WIN32)
tile_t fortnum_dtile_avx2_8x6, fortnum_dtile_avx512_24x8;
tile_c_t fortnum_dtile_avx2_8x6_c, fortnum_dtile_avx512_24x8_c;
#elif defined(__aarch64__) && !defined(_WIN32)
tile_t fortnum_dtile_neon_8x6;
tile_c_t fortnum_dtile_neon_8x6_c;
#endif

/* Returns NULL when no kernel fits this CPU; *mr and *nr give the tile and
   *update the in-place variant. */
tile_t *fortnum_blas_dtile(int32_t *mr, int32_t *nr, int32_t *isa, tile_c_t **update) {
    const char *force = getenv("FORTNUM_BLAS_ISA");
    *mr = 0;
    *nr = 0;
    *isa = 0;
    *update = NULL;
#if defined(__x86_64__) && !defined(_WIN32)
    __builtin_cpu_init();
    int has512 = __builtin_cpu_supports("avx512f");
    int has2 = __builtin_cpu_supports("avx2") && __builtin_cpu_supports("fma");
    if (force && strcmp(force, "avx2") == 0) has512 = 0;
    if (has512) {
        *mr = 24; *nr = 8; *isa = 2;
        *update = fortnum_dtile_avx512_24x8_c;
        return fortnum_dtile_avx512_24x8;
    }
    if (has2) {
        *mr = 8; *nr = 6; *isa = 1;
        *update = fortnum_dtile_avx2_8x6_c;
        return fortnum_dtile_avx2_8x6;
    }
#elif defined(__aarch64__) && !defined(_WIN32)
    (void)force;
    *mr = 8; *nr = 6; *isa = 3;
    *update = fortnum_dtile_neon_8x6_c;
    return fortnum_dtile_neon_8x6;
#else
    (void)force;
#endif
    return NULL;
}
