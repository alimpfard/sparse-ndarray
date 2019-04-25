#include "spndarray.h"
#include <math.h>
#include <stdlib.h>

#include "avl.c"

__attribute__((always_inline)) static inline size_t
array_mul(const size_t len, const size_t *arr, const ssize_t skip) {
  size_t prod = 1;
  for (size_t i = 0; i < len; i++)
    if (i != (size_t)skip)
      prod *= arr[i];
  return prod;
}

/*
 * spndarray_mul()
 *
 * calculates the dimension-wise multiplication of m and n
 * Inputs
 *  m - first sparse array of size N
 *  n - second sparse array of size N[-1]
 *  d - which dimension to multiply over (-1 to skip)
 *
 * Outputs
 *  m * n
 *
 * Notes
 *  Currently does not support fillvalues other than zero
 */
spndarray *spndarray_mul(const spndarray *xm, const spndarray *xn,
                         const size_t d) {
  spndarray *res;
  spndarray *m = (spndarray *)xm;
  spndarray *n = (spndarray *)xn;
  ssize_t f = -1;
  if (n->ndim < m->ndim) {
    void *t = m;
    m = n;
    n = t;
  }
  if (m->ndim == n->ndim) {
    if (d != (size_t)-1) // intended underflow don't kill me
      f = d;
    goto sizecheck;
  }
  if (m->ndim != n->ndim - 1) {
    fprintf(stderr, "mul requires dimensions N and N-1, but got %zd and %zd\n",
            m->ndim, n->ndim);
    return NULL;
  }
sizecheck:;
  size_t ms = array_mul(m->ndim, m->dimsizes, f);
  size_t ns = array_mul(n->ndim, n->dimsizes, d);
  if (ms != ns) {
    fprintf(stderr,
            "mul requires projected data count to be equal, but got %zd and %zd"
            " values\n",
            ms, ns);
    return NULL;
  }
  res = spndarray_alloc_nzmax(n->ndim, n->dimsizes, n->nzmax, SPNDARRAY_NTUPLE);
  size_t midx[m->ndim], nidx[n->ndim];
  // TODO reshape
  for (size_t i = 0; i < n->nz; i++) {
    size_t ni = &n->data[i] - n->data;

    for (size_t j = 0; j < n->ndim; j++) {
      nidx[j] = n->dims[j][ni];
      if (j != d)
        midx[j - (j > d)] = nidx[j];
    }
    spndarray_set(res, spndarray_get(m, midx) * spndarray_get(n, nidx), nidx);
  }
  return res;
}

/*
 * spndarray_mul_vec()
 *
 * calculates the dimension-wise multiplication of array ND m and 1D array n
 * Inputs
 *  m - first sparse array of size N
 *  n - second sparse array of size 1
 *  d - which dimension to multiply over
 *
 * Outputs
 *  m * n
 *
 * Notes
 *  Currently does not support fillvalues other than zero
 */
spndarray *spndarray_mul_vec(const spndarray *xm, const spndarray *xn,
                         const size_t d) {
  spndarray *res;
  spndarray *m = (spndarray *)xm;
  spndarray *n = (spndarray *)xn;
  if (m->ndim < n->ndim) {
    void *t = m;
    m = n;
    n = t;
  }
  if (n->ndim != 1) {
    fprintf(stderr, "mul_vec requires a one-dimension vector, but got a %zd onr\n",
            n->ndim);
    return NULL;
  }
  res = spndarray_alloc_nzmax(n->ndim, n->dimsizes, n->nzmax, SPNDARRAY_NTUPLE);
  size_t midx[m->ndim];
  for (size_t i = 0; i < m->nz; i++) {
    size_t ni = &m->data[i] - m->data;
    for (size_t j = 0; j < m->ndim; j++)
      midx[j] = m->dims[j][ni];
    spndarray_set(res, m->data[i] * spndarray_get(n, &m->dims[d][ni]), midx);
  }
  return res;
}

/*
 * spndarray_add()
 *
 * calculates the dimension-wise addition of m and n
 * Inputs
 *  m - first sparse array of size N
 *  n - second sparse array of size N[-1]
 *
 * Outputs
 *  m + n
 */
spndarray *spndarray_add(const spndarray *xm, const spndarray *xn) {
  spndarray *res;
  spndarray *m = (spndarray *)xm;
  spndarray *n = (spndarray *)xn;
  if (m->ndim != n->ndim) {
    fprintf(stderr,
            "add requires dimensions to be equal, but got %zd and %zd\n",
            m->ndim, n->ndim);
    return NULL;
  }
  size_t ms = array_mul(m->ndim, m->dimsizes, -1);
  size_t ns = array_mul(n->ndim, n->dimsizes, -1);
  if (ms != ns) {
    fprintf(stderr,
            "add requires projected data count to be equal, but got %zd and %zd"
            " values\n",
            ms, ns);
    return NULL;
  }
  res = spndarray_alloc_nzmax(n->ndim, n->dimsizes, n->nzmax, SPNDARRAY_NTUPLE);
  spndarray_set_fillvalue(res, n->fill + m->fill);
  size_t midx[m->ndim], nidx[n->ndim];
  // TODO reshape
  for (size_t i = 0; i < n->nz; i++) {
    size_t ni = &n->data[i] - n->data;

    for (size_t j = 0; j < n->ndim; j++) {
      nidx[j] = n->dims[j][ni];
      midx[j] = nidx[j];
    }
    spndarray_set(res, spndarray_get(m, midx) + spndarray_get(n, nidx), nidx);
  }

  for (size_t i = 0; i < m->nz; i++) {
    size_t mi = &m->data[i] - m->data;

    for (size_t j = 0; j < m->ndim; j++) {
      midx[j] = m->dims[j][mi];
      nidx[j] = midx[j];
    }
    spndarray_set(res, spndarray_get(m, midx) + spndarray_get(n, nidx), nidx);
  }
  return res;
}

/*
 * spndarray_sub()
 *
 * calculates the dimension-wise Subtraction of m and n
 * Inputs
 *  m - first sparse array of size N
 *  n - second sparse array of size N[-1]
 *
 * Outputs
 *  m - n
 */
spndarray *spndarray_sub(const spndarray *xm, const spndarray *xn) {
  spndarray *res;
  spndarray *m = (spndarray *)xm;
  spndarray *n = (spndarray *)xn;
  if (m->ndim != n->ndim) {
    fprintf(stderr,
            "sub requires dimensions to be equal, but got %zd and %zd\n",
            m->ndim, n->ndim);
    return NULL;
  }
  size_t ms = array_mul(m->ndim, m->dimsizes, -1);
  size_t ns = array_mul(n->ndim, n->dimsizes, -1);
  if (ms != ns) {
    fprintf(stderr,
            "sub requires projected data count to be equal, but got %zd and %zd"
            " values\n",
            ms, ns);
    return NULL;
  }
  res = spndarray_alloc_nzmax(n->ndim, n->dimsizes, n->nzmax, SPNDARRAY_NTUPLE);
  spndarray_set_fillvalue(res, n->fill + m->fill);
  size_t midx[m->ndim], nidx[n->ndim];
  // TODO reshape
  for (size_t i = 0; i < n->nz; i++) {
    size_t ni = &n->data[i] - n->data;

    for (size_t j = 0; j < n->ndim; j++) {
      nidx[j] = n->dims[j][ni];
      midx[j] = nidx[j];
    }
    spndarray_set(res, spndarray_get(m, midx) - spndarray_get(n, nidx), nidx);
  }

  for (size_t i = 0; i < m->nz; i++) {
    size_t mi = &m->data[i] - m->data;

    for (size_t j = 0; j < m->ndim; j++) {
      midx[j] = m->dims[j][mi];
      nidx[j] = midx[j];
    }
    spndarray_set(res, spndarray_get(m, midx) - spndarray_get(n, nidx), nidx);
  }
  return res;
}

/*
 * spndarray_memcpy()
 *
 * copies an array into another, or allocates a new copy
 *
 * Inputs
 *  src - the array to copy from
 *  dst - the array to copy into, or NULL to allocate a new one
 *
 * Output
 *  the destination array
 *
 * Notes
 *  if a dst is provided, it must have the same dimensionality as the
 *  source array
 *
 *  the fillvalues of the source array will _not_ be transferred over,
 *  they will be simply tranformed to the new fillvalue
 *  [doing so would cause this array to mostly dense]
 */
spndarray *spndarray_memcpy(const spndarray *src, spndarray *dst) {
  if (!dst) {
    dst = spndarray_alloc_nzmax(src->ndim, src->dimsizes, src->nzmax, SPNDARRAY_NTUPLE);
  } else if (dst->ndim != src->ndim) {
    fprintf(stderr,
            "memcpy requires dimensionality to be equal between the source and "
            "the destination array\n");
    return NULL;
  }
  size_t idx[src->ndim];
  spndarray_set_zero(dst);
  for (size_t i = 0; i < src->nz; i++) {
    for (size_t j = 0; j < src->ndim; j++)
      idx[j] = src->dims[j][i];
    spndarray_set(dst, src->data[i], idx);
  }
  return dst;
}

void spndarray_fmap(spndarray *m, double_mapper f) {
  for (size_t i = 0; i < m->nz; i++)
    m->data[i] = f(m->data[i]);
}

void spndarray_negate(spndarray *m) {
  for (size_t i = 0; i < m->nz; i++)
    m->data[i] = -m->data[i];
}
