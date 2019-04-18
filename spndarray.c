#include <stdlib.h>
#include <math.h>
#include "spndarray.h"

#include "avl.c"

static int compare_ntuple(const void *pa, const void *pb, void *param);
static void *avl_spmalloc(size_t size, void *param);
static void avl_spfree(void *block, void *param);

static struct libavl_allocator avl_allocator_spndarray = { avl_spmalloc, avl_spfree };

__attribute__((always_inline))
static inline size_t array_mul(const size_t len, const size_t* arr) {
    size_t prod = 1;
    for (size_t i = 0; i < len; i++)
        prod *= arr[i];
    return prod;
}

/*
 * spndarray_alloc()
 *
 * Allocate a sparse nd array in ntuple format
 *
 * Inputs
 *   ndims - number of dimensions
 *   dims  - list of dimension sizes
 *
 * Notes
 *   if the dim sizes are not known at allocation time, they can all
 *   be set to 1, and they will be expanded as elements are added
 */
spndarray *
spndarray_alloc(const size_t ndims, const size_t *dims) {
    const double density = 0.1; // estimate
    size_t nzmax = (size_t) floor(array_mul(ndims, dims) * density);

    if (nzmax == 0)
        nzmax = 10;

    return spndarray_alloc_nzmax(ndims, dims, nzmax, SPNDARRAY_NTUPLE);
}

/*
 * spndarray_alloc_nzmax()
 *
 * Allocate a sparse nd array with given nzmax
 *
 * TODO
 */
spndarray *spndarray_alloc_nzmax(const size_t ndims, const size_t *dimsizes,
        const size_t nzmax, const size_t flags) {
    spndarray *m;
    size_t* dimss = calloc(ndims, sizeof(size_t));
    for (size_t i = 0; i < ndims; i++)
        if (dimsizes[i] == 0) {
            fprintf(stderr, "dimension %zd must be a positive integer", i);
            abort();
        } else
            dimss[i] = dimsizes[i];

    m = calloc(1, sizeof(*m));
    if (!m) {
        fprintf(stderr, "not enough space for array");
        abort();
    }

    m->ndim = ndims;
    m->dimsizes = dimss;
    m->nz = 0;
    m->nzmax = (nzmax < 1) ? 1 : nzmax;
    m->sptype = flags;

    m->dims = calloc(ndims, sizeof(size_t*));
    if (!m->dims) {
        // error
        abort();
    }

    if (flags == SPNDARRAY_NTUPLE) {
        m->tree_data = malloc(sizeof(spndarray_tree));
        if (!m->tree_data) {
            fprintf(stderr, "not enough space for AVL tree");
            abort();
        }
        m->tree_data->n = 0;
        m->tree_data->tree = avl_create(compare_ntuple, (void *)m, &avl_allocator_spndarray);
        if(!m->tree_data->tree) {
            fprintf(stderr, "Not enough space for AVL tree");
            abort();
        }

        m->tree_data->node_array = malloc(m->nzmax * sizeof(struct avl_node));
        if (!m->tree_data->node_array) {
            fprintf(stderr, "Not enough space for AVL tree nodes");
            abort();
        }
        for (size_t i = 0; i < ndims; i++) {
            m->dims[i] = malloc(m->nzmax * sizeof(size_t));
            if (!m->dims[i]) {
                fprintf(stderr, "Not enough space for dimension %zd indices", i);
                abort();
            }
        }
    } else if (flags == SPNDARRAY_CCS) {
        // TODO
        fprintf(stderr, "SPNDARRAY_CCS not implemented");
        abort();
    }
    m->data = malloc(m->nzmax * sizeof(double));
    if (!m->data) {
        fprintf(stderr, "Not enough space for the data");
        abort();
    }

    return m;
} /* spndarray_alloc_nzmax() */

/*
 * spndarray_free()
 * Frees the given array
 */
void spndarray_free(spndarray *m) {
    if (m->dims) {
        for (size_t i = 0; i < m->ndim; i++)
            if (m->dims[i])
                free(m->dims[i]);
        free(m->dims);
    }
    if (m->data)
        free(m->data);
    if (m->dimsizes)
        free(m->dimsizes);
    if (m->work)
        free(m->work);
    if (m->tree_data) {
        if (m->tree_data->tree)
            avl_destroy(m->tree_data->tree, NULL);

        if (m->tree_data->node_array)
            free(m->tree_data->node_array);

        free(m->tree_data);
    }
    free(m);
} /* spndarray_free() */

/*
 * spndarray_realloc()
 * As elements are added to the sparse array, it's possible that they
 * will exceed the previously specified nzmax - reallocate the array
 * with a new nzmax
 */
int spndarray_realloc(const size_t nzmax, spndarray* m) {
    int s = 0;
    void *ptr;

    if (nzmax < m->nz) {
        fprintf(stderr, "new nzmax is smaller than the current nz");
        return 1;
    }

    for (size_t i = 0; i < m->ndim; i++) {
        ptr = realloc(m->dims[i], nzmax * sizeof(size_t));
        if (!ptr) {
            fprintf(stderr, "failed to allocate space for dimension %zd indices", i);
            abort();
        }
        m->dims[i] = ptr;
    }
    ptr = realloc(m->data, nzmax * sizeof(double));
    if (!ptr) {
        fprintf(stderr, "failed to allocate space for data");
        abort();
    }

    /* rebuild binary tree */
    if (SPNDARRAY_ISNTUPLE(m)) {
        size_t n;

        // reset tree to empty state, but don't free it
        avl_empty(m->tree_data->tree, NULL);
        m->tree_data->n = 0;

        ptr = realloc(m->tree_data->node_array, nzmax * sizeof(struct avl_node));
        if (!ptr) {
            fprintf(stderr, "failed to allocate space for AVL tree nodes");
            abort();
        }
        m->tree_data->node_array = ptr;
        for (n = 0; n < m->nz; ++n) {
            ptr = avl_insert(m->tree_data->tree, &m->data[n]);
            if (ptr != NULL) {
                fprintf(stderr, "detected duplicate entry while reallocating array");
                abort();
            }
        }
    }
    // update to new nzmax
    m->nzmax = nzmax;
    return s;
} /* spndarray_realloc() */

int spndarray_set_zero(spndarray *m) {
    m->nz = 0;
    if (SPNDARRAY_ISNTUPLE(m)) {
        avl_empty(m->tree_data->tree, NULL);
        m->tree_data->n = 0;
    }
    return 0;
}

size_t spndarray_nnz(const spndarray *m) {
    return m->nz;
}

/*
 * spndarray_compare_idx()
 * Comparison function for searching the binary tree
 * in ntuple format
 *
 * To detect duplicate elements in the tree, we want to determine
 * if there already exists an entry for (i0,i1,...) in the tree.
 * Since the actual tree node stores only the data elements data[n],
 * we will do pointer magick to get from the given data[n] to the
 * indices in dims[...][n].
 *
 * This compare function will sort the tree first by dim 0,
 * then dim1, and so on.
 *
 * Inputs
 *   ndim - number of dimensions
 *   ips  - element a indices
 *   ipb  - element b indices
 *
 * Return
 *   -1 if ipa < ipb
 *   +1 if ipa > ipb
 *    0 if ipa = ipb
 */
int spndarray_compare_idx(const size_t ndims, const size_t* ipa, const size_t* ipb) {
  for (size_t dim = 0 ; dim < ndims; dim++) {
      if (ipa[dim] < ipb[dim])
          return -1;
      else if (ipa[dim] > ipb[dim])
          return 1;
  }
  return 0; // all equal
}

/*
 * spndarray_tree_rebuild()
 * When copying a ntuple array, it is necessary to rebuild
 * the binary tree for element searches
 *
 * Input : m - ntuple array
 */
int spndarray_tree_rebuild(spndarray *m) {
    if (!SPNDARRAY_ISNTUPLE(m)) {
        fprintf(stderr, "m must be in ntuple format");
        return 1;
    }
    size_t n;

    // reset tree to be empty, but leave the root ptr;
    avl_empty(m->tree_data->tree, NULL);
    m->tree_data->n = 0;

    // insert all tree elements
    for (n = 0; n < m->nz; m++) {
        void *ptr = avl_insert(m->tree_data->tree, &m->data[n]);
        if (ptr != NULL) {
            fprintf(stderr, "duplicate entry detected while rebuilding tree");
            return 1;
        }
    }
    return 0;
}

/*
 * compare_ntuple()
 * Comparison function for searching binary tree in
 * ntuple format representation.
 *
 * To detect duplicate entries in the tree, we want
 * to determine if there already exists an entry for (i0,i1,...)
 * in the tree. Since the actual tree node stores only the data
 * elements data[n], we will do pointer magick to get from the
 * given data[n] to the indices dims[...][n]
 *
 * This compare function will sort the tree first by dim 0,
 * then dim1, and so on.
 */

static int compare_ntuple(const void *pa, const void *pb, void *params) {
    spndarray *m = (spndarray *) params;

    const size_t idxa = (const double *) pa - m->data;
    const size_t idxb = (const double *) pb - m->data;

    size_t ipa[m->ndim], ipb[m->ndim];
    for (size_t i = 0; i < m->ndim; i++) {
        ipa[i] = m->dims[i][idxa];
        ipb[i] = m->dims[i][idxb];
    }
    return spndarray_compare_idx(m->ndim, ipa, ipb);
}

static void *
avl_spmalloc (size_t size, void *param) {
    spndarray *m = (spndarray *) param;

    if (size != sizeof(struct avl_node)) {
        fprintf(stderr, "attempting to allocate incorrect node size");
        return NULL;
    }

    // return the next available avl_node slot; index
    // m->tree_data->n keeps track of the next open slot
    if (m->tree_data->n < m->nzmax) {
        unsigned char *node_ptr = (unsigned char *) m->tree_data->node_array;

        // offset in bytes for the next node slot
        size_t offset = (m->tree_data->n)++ * sizeof(struct avl_node);

        return node_ptr + offset;
    } else {
        /* we should never get here - spndarray_realloc() should
         * be called before exceeding nzmax nodes
         */
        fprintf(stderr, "attempting to allocate tree node past nzmax");
        __builtin_unreachable();
    }
}

static void
avl_spfree (void *block, void *params) {
    (void) block;
    (void) params;
    /*
     * do nothing - instead of alloc/free'ing individual nodes.
     * we malloc and free nzmax nodes at a time
     */
}
