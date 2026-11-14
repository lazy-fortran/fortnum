/* The GEMM register tiles for the running CPU, as function pointers, so that
   no kernel of another architecture is referenced: the tile kernel and the
   same tile adding to C in place. FORTNUM_BLAS_ISA (avx2) forces the AVX2
   tile on a CPU that also has AVX-512. */
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

typedef void tile_t(int64_t k, const double *a, const double *b, double *acc);
typedef void tile_c_t(int64_t k, const double *a, const double *b, double *c, int64_t ldc);
typedef void tile_s_t(int64_t k, const double *a, int64_t lda, const double *b, double *c, int64_t ldc);

#if defined(__x86_64__) && !defined(_WIN32)
tile_t fortnum_dtile_avx2_8x6, fortnum_dtile_avx512_24x8;
tile_c_t fortnum_dtile_avx2_8x6_c, fortnum_dtile_avx512_24x8_c;
tile_s_t fortnum_dtile_avx2_8x6_sc, fortnum_dtile_avx512_24x8_sc, fortnum_dtile_avx2_8x6_sa, fortnum_dtile_avx512_24x8_sa;
#elif defined(__aarch64__) && !defined(_WIN32)
tile_t fortnum_dtile_neon_8x6;
tile_c_t fortnum_dtile_neon_8x6_c;
tile_s_t fortnum_dtile_neon_8x6_sc, fortnum_dtile_neon_8x6_sa;
#endif

/* Returns NULL when no kernel fits this CPU; *mr and *nr give the tile,
   *update the in-place variant, *strided_update and *strided_acc the
   variants reading A unpacked with a leading dimension. */
tile_t *fortnum_blas_dtile(int32_t *mr, int32_t *nr, int32_t *isa, tile_c_t **update,
                           tile_s_t **strided_update, tile_s_t **strided_acc) {
    const char *force = getenv("FORTNUM_BLAS_ISA");
    *mr = 0;
    *nr = 0;
    *isa = 0;
    *update = NULL;
    *strided_update = NULL;
    *strided_acc = NULL;
#if defined(__x86_64__) && !defined(_WIN32)
    __builtin_cpu_init();
    int has512 = __builtin_cpu_supports("avx512f");
    int has2 = __builtin_cpu_supports("avx2") && __builtin_cpu_supports("fma");
    if (force && strcmp(force, "avx2") == 0) has512 = 0;
    if (has512) {
        *mr = 24; *nr = 8; *isa = 2;
        *update = fortnum_dtile_avx512_24x8_c;
        *strided_update = fortnum_dtile_avx512_24x8_sc;
        *strided_acc = fortnum_dtile_avx512_24x8_sa;
        return fortnum_dtile_avx512_24x8;
    }
    if (has2) {
        *mr = 8; *nr = 6; *isa = 1;
        *update = fortnum_dtile_avx2_8x6_c;
        *strided_update = fortnum_dtile_avx2_8x6_sc;
        *strided_acc = fortnum_dtile_avx2_8x6_sa;
        return fortnum_dtile_avx2_8x6;
    }
#elif defined(__aarch64__) && !defined(_WIN32)
    (void)force;
    *mr = 8; *nr = 6; *isa = 3;
    *update = fortnum_dtile_neon_8x6_c;
    *strided_update = fortnum_dtile_neon_8x6_sc;
    *strided_acc = fortnum_dtile_neon_8x6_sa;
    return fortnum_dtile_neon_8x6;
#else
    (void)force;
#endif
    return NULL;
}
