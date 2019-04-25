#ifndef __SPNDARRAY_H__
#define __SPNDARRAY_H__

#include <stdlib.h>

#ifdef __cplusplus
#define __BEGIN_DECLS extern "C" {
#define __END_DECLS }
#else
#define __BEGIN_DECLS /* empty */
#define __END_DECLS   /* empty */
#endif

__BEGIN_DECLS

/*
 * Binary tree data structure for storing sparse matrix elements
 * in N-tuple format. This is used for efficient detection of
 * duplicates and element retrieval
 */
typedef struct {
  void *tree;       /* tree structure */
  void *node_array; /* preallocated array of tree nodes */
  size_t n;         /* number of tree nodes in use (<= nzmax) */
} spndarray_tree;

/*
 * N-tuple format:
 *
 * if data[n] = A_{i_0i_1i_2...}, then
 *     i_x = A->dims[x][n]
 *
 * Compressed Column Format (CCS):
 *
 * If data[n] = A_{i_0i_1i_2...}, then
 *     i_{j = 0...m-1} = A->dims[j][n]
 *     A->dims[A->ndim-1][m-1] <= n < A->dims[A->ndim-1][m]
 *
 * such that column i_{ndim-1} is stored in
 *     [ data[A->dims[A->ndim-1][m-1]], data[A->dims[A->ndim-1][m-1] + 1], ...,
 * data[A->dims[A->ndim-1][m] - 1] ]
 */
typedef struct {
  size_t *sizes; /* number of indices */
  size_t ndim;   /* number of dimensions */
  double *data;  /* space for elements */

  /* dimsizes (size ndim) contains
   *
   * List of dimension sizes
   */
  size_t *dimsizes;

  /* dims (size ndim) contains
   *
   * List of dimension indices
   */
  size_t **dims;

  size_t nzmax; /* maximum number of array elements */
  size_t nz;    /* current number of non-fillvalue elements */
  double fill;  /* fill value of the array */
  spndarray_tree *tree_data; /* binary tree for sorting N-Tuple data */

  /*
   * workspace of size MAX{sizes} * MAX{sizeof(double), sizeof(size_t)}
   * used in some routine
   */
  union {
    void *work;
    size_t *work_sze;
    double *work_dbl;
  };

  size_t sptype; /* storage type */
} spndarray;

#define SPNDARRAY_NTUPLE (0)
#define SPNDARRAY_CCS (1)

#define SPNDARRAY_ISNTUPLE(m) ((m)->sptype == SPNDARRAY_NTUPLE)
#define SPNDARRAY_ISCCS(m) ((m)->sptype == SPNDARRAY_CCS)

typedef double (*reduction_function)(double acc, double x, int count);
typedef double (*double_mapper)(double value);

/*
 * Prototypes
 */

spndarray *spndarray_alloc(const size_t ndims, const size_t *dimsizes);
spndarray *spndarray_alloc_nzmax(const size_t ndims, const size_t *dimsizes,
                                 const size_t nzmax, const size_t flags);

void spndarray_set_fillvalue(spndarray *m, const double fill);

void spndarray_free(spndarray *m);
int spndarray_realloc(const size_t nzmax, spndarray *m);
int spndarray_set_zero(spndarray *m);
size_t spndarray_nnz(const spndarray *m);

int spndarray_compare_idx(const size_t ndims, const size_t *adims,
                          const size_t *bdims);
int spndarray_tree_rebuild(spndarray *m);

/* spndcopy.c */
spndarray *spndarray_memcpy(const spndarray *src, spndarray *dst);

/* spndgetset.c */
double spndarray_get(const spndarray *m, const size_t *idxs);
double spndarray_getv(const spndarray *m, ...);

int spndarray_set(spndarray *m, double val, const size_t *idxs);
int spndarray_setv(spndarray *m, double val, ...);

double *spndarray_ptr(const spndarray *m, const size_t *idxs);
double *spndarray_ptrv(const spndarray *m, ...);

void spndarray_incr(spndarray *m, const size_t *idxs);
void spndarray_incrv(spndarray *m, ...);

/* spreduce.c */

double reduce_sum(double acc, double x, int count);
double reduce_mean(double acc, double x, int count);

spndarray *spndarray_reduce(spndarray *m, const size_t dim, const reduction_function reduce_fn);
spndarray *spndarray_extract_dimension(spndarray *m, const size_t dim);
// TODO compress, io, operations, prop, swap

/* spndio.c */
int spndarray_fwrite(const spndarray* m, const char* fmt, const char* filepath, const int sparse);

/* spndop.c */
spndarray *spndarray_mul(const spndarray *m, const spndarray *n, const size_t d);
spndarray *spndarray_mul_vec(const spndarray *m, const spndarray *n, const size_t d);
spndarray *spndarray_add(const spndarray *m, const spndarray *n);
void spndarray_mulinverse(spndarray *m);

__END_DECLS
#endif
