/* The GEMM register tiles for the running CPU, as function pointers, so that
   no kernel of another architecture is referenced: the tile kernel and the
   same tile adding to C in place. FORTNUM_BLAS_ISA forces the AVX2 tile on a
   CPU that also has AVX-512 (avx2) or the NEON tile on an Apple CPU with the
   AMX matrix coprocessor (neon). */
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

typedef void tile_t(int64_t k, const double *a, const double *b, double *acc);
typedef void tile_c_t(int64_t k, const double *a, const double *b, double *c, int64_t ldc);
typedef void tile_s_t(int64_t k, const double *a, int64_t lda, const double *b, double *c, int64_t ldc);
typedef void tile_b_t(int64_t k, const double *a, const double *b, int64_t ldb, double *c, int64_t ldc);

#if defined(__x86_64__) && !defined(_WIN32)
tile_t fortnum_dtile_avx2_8x6, fortnum_dtile_avx512_24x8;
tile_c_t fortnum_dtile_avx2_8x6_c, fortnum_dtile_avx512_24x8_c;
tile_s_t fortnum_dtile_avx2_8x6_sc, fortnum_dtile_avx512_24x8_sc, fortnum_dtile_avx2_8x6_sa, fortnum_dtile_avx512_24x8_sa;
tile_b_t fortnum_dtile_avx2_8x6_bc, fortnum_dtile_avx512_24x8_bc, fortnum_dtile_avx2_8x6_ba, fortnum_dtile_avx512_24x8_ba;
#elif defined(__aarch64__) && !defined(_WIN32)
tile_t fortnum_dtile_neon_8x6;
tile_c_t fortnum_dtile_neon_8x6_c;
tile_s_t fortnum_dtile_neon_8x6_sc, fortnum_dtile_neon_8x6_sa;
tile_b_t fortnum_dtile_neon_8x6_bc, fortnum_dtile_neon_8x6_ba;
#elif defined(__wasm32__) && defined(__wasm_simd128__)
tile_t fortnum_dtile_wasm_8x4;
tile_c_t fortnum_dtile_wasm_8x4_c;
tile_s_t fortnum_dtile_wasm_8x4_sc, fortnum_dtile_wasm_8x4_sa;
tile_b_t fortnum_dtile_wasm_8x4_bc, fortnum_dtile_wasm_8x4_ba;
#endif
#if defined(__aarch64__) && defined(__APPLE__)
#include <setjmp.h>
#include <signal.h>
tile_t fortnum_dtile_amx_32x16;
tile_c_t fortnum_dtile_amx_32x16_c;
tile_s_t fortnum_dtile_amx_32x16_sc, fortnum_dtile_amx_32x16_sa;
static sigjmp_buf amx_jump;
static void amx_missing(int sig) {
    (void)sig;
    siglongjmp(amx_jump, 1);
}
/* Apple publishes no feature flag for AMX: enable and disable it once and
   catch the SIGILL of a CPU without it (M4 and later have SME instead). */
static int amx_present(void) {
    struct sigaction catch_ill, previous;
    int present = 0;
    catch_ill.sa_handler = amx_missing;
    catch_ill.sa_flags = 0;
    sigemptyset(&catch_ill.sa_mask);
    sigaction(SIGILL, &catch_ill, &previous);
    if (!sigsetjmp(amx_jump, 1)) {
        __asm__ volatile("nop\nnop\nnop\n.word 0x201220\nnop\nnop\nnop\n.word 0x201221" ::: "memory");
        present = 1;
    }
    sigaction(SIGILL, &previous, NULL);
    return present;
}
#endif

/* Returns NULL when no kernel fits this CPU; *mr and *nr give the tile,
   *update the in-place variant, *strided_update and *strided_acc the
   variants reading A unpacked with a leading dimension. */
tile_t *fortnum_blas_dtile(int32_t *mr, int32_t *nr, int32_t *isa, tile_c_t **update,
                           tile_s_t **strided_update, tile_s_t **strided_acc,
                           tile_b_t **bstrided_update, tile_b_t **bstrided_acc) {
    const char *force = getenv("FORTNUM_BLAS_ISA");
    *mr = 0;
    *nr = 0;
    *isa = 0;
    *update = NULL;
    *strided_update = NULL;
    *strided_acc = NULL;
    *bstrided_update = NULL;
    *bstrided_acc = NULL;
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
        *bstrided_update = fortnum_dtile_avx512_24x8_bc;
        *bstrided_acc = fortnum_dtile_avx512_24x8_ba;
        return fortnum_dtile_avx512_24x8;
    }
    if (has2) {
        *mr = 8; *nr = 6; *isa = 1;
        *update = fortnum_dtile_avx2_8x6_c;
        *strided_update = fortnum_dtile_avx2_8x6_sc;
        *strided_acc = fortnum_dtile_avx2_8x6_sa;
        *bstrided_update = fortnum_dtile_avx2_8x6_bc;
        *bstrided_acc = fortnum_dtile_avx2_8x6_ba;
        return fortnum_dtile_avx2_8x6;
    }
#elif defined(__aarch64__) && !defined(_WIN32)
#if defined(__APPLE__)
    if (!(force && strcmp(force, "neon") == 0) && amx_present()) {
        *mr = 32; *nr = 16; *isa = 4;
        *update = fortnum_dtile_amx_32x16_c;
        *strided_update = fortnum_dtile_amx_32x16_sc;
        *strided_acc = fortnum_dtile_amx_32x16_sa;
        return fortnum_dtile_amx_32x16;
    }
#endif
    (void)force;
    *mr = 8; *nr = 6; *isa = 3;
    *update = fortnum_dtile_neon_8x6_c;
    *strided_update = fortnum_dtile_neon_8x6_sc;
    *strided_acc = fortnum_dtile_neon_8x6_sa;
    *bstrided_update = fortnum_dtile_neon_8x6_bc;
    *bstrided_acc = fortnum_dtile_neon_8x6_ba;
    return fortnum_dtile_neon_8x6;
#elif defined(__wasm32__) && defined(__wasm_simd128__)
    (void)force;
    *mr = 8; *nr = 4; *isa = 5;
    *update = fortnum_dtile_wasm_8x4_c;
    *strided_update = fortnum_dtile_wasm_8x4_sc;
    *strided_acc = fortnum_dtile_wasm_8x4_sa;
    *bstrided_update = fortnum_dtile_wasm_8x4_bc;
    *bstrided_acc = fortnum_dtile_wasm_8x4_ba;
    return fortnum_dtile_wasm_8x4;
#else
    (void)force;
#endif
    return NULL;
}
