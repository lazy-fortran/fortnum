#include <stdint.h>
#include <sys/resource.h>

int64_t fortnum_peak_rss_bytes(void) {
    struct rusage usage;
    if (getrusage(RUSAGE_SELF, &usage) != 0) {
        return -1;
    }
#if defined(__APPLE__)
    return (int64_t)usage.ru_maxrss;
#else
    return (int64_t)usage.ru_maxrss * 1024;
#endif
}
