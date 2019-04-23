#include <stdio.h>
#include <stdlib.h>

#include "spndarray.h"

static void test_getset() {
  printf(">> Running %s <<\n\n", __FUNCTION__);
  const double val[] = {4.454, 324, -1231231};
  spndarray *m =
      spndarray_alloc_nzmax(3, (size_t[]){1, 100, 10}, 10, SPNDARRAY_NTUPLE);
  for (int i = 0; i < 3; i++)
    spndarray_set(m, val[i], (size_t[]){0, i, 1});
  for (int i = 0; i < 3; i++)
    printf("idx %d value put: %f, value got: %f\n", i, val[i],
           spndarray_get(m, (size_t[]){0, i, 1}));
  spndarray_free(m);
  printf("\n>> %s Finished <<\n\n", __FUNCTION__);
}

static void test_incr() {
  printf(">> Running %s <<\n\n", __FUNCTION__);
  const double val[] = {4.454, 324, -1231231};
  spndarray *m =
      spndarray_alloc_nzmax(3, (size_t[]){1, 100, 10}, 10, SPNDARRAY_NTUPLE);
  for (int i = 0; i < 3; i++)
    spndarray_set(m, val[i], (size_t[]){0, i, 1});
  for (int i = 0; i < 3; i++)
    spndarray_incr(m, (size_t[]){0, i, 1});

  for (int i = 0; i < 3; i++)
    printf("idx %d value put: %f, value got: %f\n", i, val[i],
           spndarray_get(m, (size_t[]){0, i, 1}));
  spndarray_free(m);
  printf("\n>> %s Finished <<\n\n", __FUNCTION__);
}

static double fmean(double acc, double x, int c) {
  printf("acc: %f, x: %f, c: %d\n", acc, x, c);
  return acc + x / c;
}

static void test_reduce() {
  printf(">> Running %s <<\n\n", __FUNCTION__);
  const double val[] = {4.454, 324, -1231231};
  spndarray *m =
      spndarray_alloc_nzmax(3, (size_t[]){2, 10, 10}, 10, SPNDARRAY_NTUPLE);
  for (int i = 0; i < 2; i++)
    for (int j = 0; j < 10; j++)
      for (int k = 0; k < 10; k++)
        if (!j || (j && (i + k) % j == 0))
          spndarray_set(m, val[(i + j + k) % 3], (size_t[]){i, j, k});
  spndarray *mm = spndarray_reduce(m, 0, fmean);
  for (int i = 0; i < 2; i++)
    for (int j = 0; j < 10; j++)
      for (int k = 0; k < 10; k++)
        printf("idx (%d,)%d,%d value put: %f, value got: %f\n", i, j, k,
               spndarray_get(m, (size_t[]){i, j, k}),
               spndarray_get(mm, (size_t[]){j, k}));
  spndarray_free(m);
  spndarray_free(mm);
  printf("\n>> %s Finished <<\n\n", __FUNCTION__);
}

static void test_op() {
  printf(">> Running %s <<\n\n", __FUNCTION__);
  const double val[] = {4.454, 32, -12, 2, 3, -1};
  spndarray *m = spndarray_alloc_nzmax(3, (size_t[]){2, 10, 10}, 10,
                                       SPNDARRAY_NTUPLE),
            *n = spndarray_alloc_nzmax(3, (size_t[]){2, 10, 10}, 10,
                                       SPNDARRAY_NTUPLE);
  spndarray_set_fillvalue(n, 1);
  for (int i = 0; i < 2; i++)
    for (int j = 0; j < 10; j++)
      for (int k = 0; k < 10; k++)
        if (!j || (j && (i + k) % j == 0)) {
          spndarray_set(m, val[(i + j + k) % 3], (size_t[]){i, j, k});
          spndarray_set(n, val[(i + j + k) % 3 + 3], (size_t[]){i, j, k});
        }
  spndarray *nm = spndarray_mul(n, m, -1);
  printf("SPND mul has %zd nonzero elements\n", nm->nz);
  for (int i = 0; i < 2; i++)
    for (int j = 0; j < 10; j++)
      for (int k = 0; k < 10; k++)
        printf("idx %d,%d,%d %03.3f * %03.3f, value got: %03.3f\n", i, j, k,
               spndarray_get(m, (size_t[]){i, j, k}),
               spndarray_get(n, (size_t[]){i, j, k}),
               spndarray_get(nm, (size_t[]){i, j, k}));
  spndarray_free(nm);
  nm = spndarray_add(n, m);
  printf("SPND add has %zd nonzero elements\n", nm->nz);
  for (int i = 0; i < 2; i++)
    for (int j = 0; j < 10; j++)
      for (int k = 0; k < 10; k++)
        printf("idx %d,%d,%d %03.3f + %03.3f, value got: %03.3f\n", i, j, k,
               spndarray_get(m, (size_t[]){i, j, k}),
               spndarray_get(n, (size_t[]){i, j, k}),
               spndarray_get(nm, (size_t[]){i, j, k}));
  spndarray *copyN = spndarray_memcpy(m, NULL);
  printf("the fv=%f copy has %zd nonzero elements\n", copyN->fill, copyN->nz);
  copyN->fill = 4.454; // should ignore them?
  spndarray *copyV = spndarray_memcpy(m, copyN);
  spndarray_free(m);
  spndarray_free(n);
  spndarray_free(nm);
  printf("the fv=%f copy has %zd nonzero elements\n", copyV->fill, copyV->nz);
  spndarray_free(copyV);
  printf("\n>> %s Finished <<\n\n", __FUNCTION__);
}

int main() {
  test_getset();
  test_incr();
  test_reduce();
  test_op();
}
