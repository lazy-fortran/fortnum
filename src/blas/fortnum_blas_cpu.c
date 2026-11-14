/* The GEMM register tile for the running CPU, as a function pointer, so that
   no kernel of another architecture is referenced. FORTNUM_BLAS_ISA
   (avx2, avx512, neon) forces a choice the CPU supports. */
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

typedef void tile_t(int64_t k, const double *a, const double *b, double *acc);

#if defined(__x86_64__) && !defined(_WIN32)
tile_t fortnum_dtile_avx2_8x6, fortnum_dtile_avx512_24x8;
#elif defined(__aarch64__) && !defined(_WIN32)
tile_t fortnum_dtile_neon_8x6;
#endif

/* Returns NULL when no kernel fits this CPU; *mr and *nr give the tile. */
tile_t *fortnum_blas_dtile(int32_t *mr, int32_t *nr, int32_t *isa) {
    const char *force = getenv("FORTNUM_BLAS_ISA");
    *mr = 0;
    *nr = 0;
    *isa = 0;
#if defined(__x86_64__) && !defined(_WIN32)
    __builtin_cpu_init();
    int has512 = __builtin_cpu_supports("avx512f");
    int has2 = __builtin_cpu_supports("avx2") && __builtin_cpu_supports("fma");
    if (force && strcmp(force, "avx2") == 0) has512 = 0;
    if (has512) {
        *mr = 24; *nr = 8; *isa = 2;
        return fortnum_dtile_avx512_24x8;
    }
    if (has2) {
        *mr = 8; *nr = 6; *isa = 1;
        return fortnum_dtile_avx2_8x6;
    }
#elif defined(__aarch64__) && !defined(_WIN32)
    (void)force;
    *mr = 8; *nr = 6; *isa = 3;
    return fortnum_dtile_neon_8x6;
#else
    (void)force;
#endif
    return NULL;
}
