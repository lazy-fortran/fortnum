/* Counts calls of malloc, calloc and realloc, linked with
   -Wl,--wrap=malloc,--wrap=calloc,--wrap=realloc. */
#include <stddef.h>

void *__real_malloc(size_t size);
void *__real_calloc(size_t count, size_t size);
void *__real_realloc(void *pointer, size_t size);

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

int fortnum_test_allocations(void) { return __atomic_load_n(&calls, __ATOMIC_RELAXED); }
