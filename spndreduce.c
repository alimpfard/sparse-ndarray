#include "spndarray.h"
#include <math.h>
#include <stdlib.h>

#include "avl.c"

double reduce_sum(double acc, double x, int count) {
  (void)count;
  return acc + x;
}
double reduce_mean(double acc, double x, int count) {
  return acc + (x / count);
}

/*
 * spndarray_reduce()
 * Reduce a dimension by applying a reduction function over it
 *
 * Inputs
 *   dim       - which dimension to reduce over
 *   reduce_fn - the function with which to reduce
 *
 * Output
 *   the new reduced spndarray.
 */
spndarray *spndarray_reduce(spndarray *m, const size_t dim,
                            const reduction_function reduce_fn) {
  if (!SPNDARRAY_ISNTUPLE(m)) {
    fprintf(stderr, "array must be in the ntuple format");
    return NULL;
  }
  // allocate a new spndarray that is missing the given dimension
  size_t mndim = m->ndim, ndim = mndim - 1; // remove one
  size_t dims[ndim], tidx[ndim];
  size_t counters[mndim];
  memset(counters, 0, sizeof(counters));
  memset(tidx, 0, sizeof(tidx));
  size_t lastdimsize = m->dimsizes[ndim];

  for (size_t i = 0, j = 0; i < mndim; i++, j++) {
    if (i != dim) {
      dims[j] = m->dimsizes[i];
    } else
      j--;
  }

  size_t rdimsize = m->dimsizes[dim];
  spndarray *newm =
      spndarray_alloc_nzmax(ndim, dims, m->nzmax, SPNDARRAY_NTUPLE);
  while (counters[ndim] < lastdimsize) {
    size_t i;
    /*
     *  pidx  <- (i, j, k, ...)
     *  x     <- get(m, pidx)
     *  tidx  <- drop(pidx, dim)
     *  pacc  <- ptr(newm, tidx, 0)
     *  *pacc <- reduce_fn(*pacc, x, rdimsize);
     */
    double x = spndarray_get(m, counters);
    if (x == 0.0)
      goto next;
    double *pacc = spndarray_ptr(newm, tidx);
    if (!pacc) {
      spndarray_set(newm, reduce_fn(0, x, rdimsize), tidx);
      goto next;
    }
    *pacc = reduce_fn(*pacc, x, rdimsize);

  next:;
    for (i = 0; counters[i] == m->dimsizes[i]; i++)
      counters[i] = 0;
    ++counters[i];
    size_t j;

    for (i = 0, j = 0; i < mndim; i++, j++)
      if (i == dim)
        j--;
      else
        tidx[j] = counters[i];
  }
  return newm;
}

/*
 * spndarray_extract_dimension()
 *
 * extract the values of a given dimension
 *
 * Inputs
 *  dim - the dimension to extract
 */
spndarray *spndarray_extract_dimension(spndarray *m, size_t dim) {
  spndarray *ex = spndarray_alloc_nzmax(1, &m->dimsizes[dim], m->dimsizes[dim] * 0.1, SPNDARRAY_NTUPLE);
  for (size_t i = 0; i < m->nz; i++) {
    double val = m->data[i];
    spndarray_set(ex, val, &i);
  }
  return ex;
}
