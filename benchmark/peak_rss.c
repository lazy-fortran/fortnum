#include <stdint.h>
#if defined(__EMSCRIPTEN__)
#include <emscripten/heap.h>
#else
#include <sys/resource.h>
#endif

int64_t fortnum_peak_rss_bytes(void) {
#if defined(__EMSCRIPTEN__)
    /* WebAssembly memory only grows, so its size is the peak. */
    return (int64_t)emscripten_get_heap_size();
#else
    struct rusage usage;
    if (getrusage(RUSAGE_SELF, &usage) != 0) {
        return -1;
    }
#if defined(__APPLE__)
    return (int64_t)usage.ru_maxrss;
#else
    return (int64_t)usage.ru_maxrss * 1024;
#endif
#endif
}
