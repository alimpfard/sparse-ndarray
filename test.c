#include <stdio.h>
#include <stdlib.h>

#include "spndarray.h"

static void test_getset() {
    const double val[] = {4.454, 324, -1231231};
    spndarray *m = spndarray_alloc_nzmax(3, (size_t[]){1, 100, 10}, 10, SPNDARRAY_NTUPLE);
    for (int i = 0; i < 3; i++)
        spndarray_set(m, val[i], (size_t[]){0, i, 1});
    for (int i = 0; i < 3; i++)
        printf("idx %d value put: %f, value got: %f\n", i, val[i], spndarray_get(m, (size_t[]){0, i, 1}));
}

int main() {
    test_getset();
}
