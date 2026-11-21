/* Counts calls of malloc, calloc and realloc, linked with
   -Wl,--wrap=malloc,--wrap=calloc,--wrap=realloc. Without the wrap (fpm
   links every test C source into every test) the __real_ symbols stay
   unresolved weak references and the count reports -1. */
#include <stddef.h>
#include <stdlib.h>

__attribute__((weak)) void *__real_malloc(size_t size);
__attribute__((weak)) void *__real_calloc(size_t count, size_t size);
__attribute__((weak)) void *__real_realloc(void *pointer, size_t size);

static int calls;

void *__wrap_malloc(size_t size) {
    __atomic_add_fetch(&calls, 1, __ATOMIC_RELAXED);
    return __real_malloc(size);
}

void *__wrap_calloc(size_t count, size_t size) {
    __atomic_add_fetch(&calls, 1, __ATOMIC_RELAXED);
    return __real_calloc(count, size);
}

void *__wrap_realloc(void *pointer, size_t size) {
    __atomic_add_fetch(&calls, 1, __ATOMIC_RELAXED);
    return __real_realloc(pointer, size);
}

int fortnum_test_allocations(void) {
    int before = __atomic_load_n(&calls, __ATOMIC_RELAXED);
    void *probe = malloc(1);
    int after = __atomic_load_n(&calls, __ATOMIC_RELAXED);
    free(probe);
    if (after == before) return -1;
    __atomic_sub_fetch(&calls, 1, __ATOMIC_RELAXED);
    return before;
}
