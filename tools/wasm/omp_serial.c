/* The OpenMP runtime functions fortnum's tests call, for WebAssembly, which
   has no OpenMP runtime and runs a program on one thread. */
#include <time.h>

int omp_get_max_threads(void) { return 1; }
int omp_get_num_threads(void) { return 1; }
int omp_get_thread_num(void) { return 0; }
int omp_get_num_procs(void) { return 1; }
int omp_in_parallel(void) { return 0; }
void omp_set_num_threads(int n) { (void)n; }

double omp_get_wtime(void)
{
    struct timespec t;
    clock_gettime(CLOCK_MONOTONIC, &t);
    return (double)t.tv_sec + 1.0e-9 * (double)t.tv_nsec;
}
