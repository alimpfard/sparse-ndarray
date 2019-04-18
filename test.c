#include <stdio.h>
#include <stdlib.h>

#include "spndarray.h"

static void test_getset() {
  printf(">> Running %s <<\n\n", __FUNCTION__);
    const double val[] = {4.454, 324, -1231231};
    spndarray *m = spndarray_alloc_nzmax(3, (size_t[]){1, 100, 10}, 10, SPNDARRAY_NTUPLE);
    for (int i = 0; i < 3; i++)
        spndarray_set(m, val[i], (size_t[]){0, i, 1});
    for (int i = 0; i < 3; i++)
        printf("idx %d value put: %f, value got: %f\n", i, val[i], spndarray_get(m, (size_t[]){0, i, 1}));
    spndarray_free(m);
    printf("\n>> %s Finished <<\n\n", __FUNCTION__);
}

static void test_incr() {
  printf(">> Running %s <<\n\n", __FUNCTION__);
    const double val[] = {4.454, 324, -1231231};
    spndarray *m = spndarray_alloc_nzmax(3, (size_t[]){1, 100, 10}, 10, SPNDARRAY_NTUPLE);
    for (int i = 0; i < 3; i++)
        spndarray_set(m, val[i], (size_t[]){0, i, 1});
    for (int i = 0; i < 3; i++)
        spndarray_incr(m, (size_t[]){0, i, 1});

    for (int i = 0; i < 3; i++)
        printf("idx %d value put: %f, value got: %f\n", i, val[i], spndarray_get(m, (size_t[]){0, i, 1}));
    spndarray_free(m);
    printf("\n>> %s Finished <<\n\n", __FUNCTION__);
}


int main() {
    test_getset();
    test_incr();
}
