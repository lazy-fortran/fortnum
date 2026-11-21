/* Where the linker cannot wrap malloc, allocations are not counted. */
__attribute__((weak)) int fortnum_test_allocations(void) { return -1; }
