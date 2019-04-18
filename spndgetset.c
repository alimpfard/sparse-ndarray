#include <stdlib.h>
#include <math.h>
#include "spndarray.h"

#include "avl.c"

static void *tree_find(const spndarray *m, const size_t ndim, const size_t* idxs);

double spndarray_get(const spndarray *m, const size_t *idxs) {
    if (m->nz == 0)
        return 0.0;

    // out of order...?
    for (size_t i = 0; i < m->ndim; i++)
        if (idxs[i] >= m->dimsizes[i])
            return 0.0;

    if (SPNDARRAY_ISNTUPLE(m)) {
        void *ptr = tree_find(m, m->ndim, idxs);

        double x = ptr ? *(double *) ptr : 0.0;

        return x;
    } else {
        // TODO
        fprintf(stderr, "Not implemented");
        abort();
    }
    // ... how did we get here?
    return 0.0;
}

int spndarray_set(spndarray *m, const double x, const size_t* idxs) {
    if (!SPNDARRAY_ISNTUPLE(m)) {
        fprintf(stderr, "array not in ntuple format");
        return 1;
    } else if (x == 0.0) {
        void *ptr = tree_find(m, m->ndim, idxs);

        /*
         * just set the data element to 0; it'd be simple to
         * delete the node from the avl tree, but deleting the 
         * data from ->data is not so simple
         */
        if (ptr)
            *(double *) ptr = 0.0;

        return 0;
    } else {
        int s = 0;
        if (m->nz >= m->nzmax) {
            s = spndarray_realloc(2 * m->nzmax, m);
            if (s)
                return s;
        }

        // store the ntuple
        for (size_t i = 0; i < m->ndim; i++)
            m->dims[i][m->nz] = idxs[i];

        m->data[m->nz] = x;

        void *ptr = avl_insert(m->tree_data->tree, &m->data[m->nz]);
        if (ptr != NULL) {
            // found duplicate entry, replace it
            *(double *) ptr = x;
        } else {
            // no duplicate found, update indices as needed
            //
            // increase dimensions if needed
            for (size_t i = 0; i < m->ndim; i++)
                m->dimsizes[i] = (m->dimsizes[i] > idxs[i] + 1) ? m->dimsizes[i] : idxs[i] + 1;

            ++(m->nz);
        }
    return s;
    }

}

double *spndarray_ptr (const spndarray *m, const size_t* idxs) {
    for (size_t i = 0; i < m->ndim; i++)
        if (idxs[i] >= m->dimsizes[i])
            return NULL;

    if (SPNDARRAY_ISNTUPLE(m)) {
        void *ptr = tree_find(m, m->ndim, idxs);
        return (double *) ptr;
    } else {
        // TODO
        fprintf(stderr, "Not implemented");
        abort();
    }
}

static void *tree_find(const spndarray *m, const size_t ndim, const size_t *idxs) {
    const struct avl_table *tree = (struct avl_table *) m->tree_data->tree;
    const struct avl_node *p;
    
    for (p = tree->avl_root; p != NULL; ) {
        size_t n = (double *) p->avl_data - m->data;
        size_t pi[ndim];
        for (size_t i = 0; i < ndim; i++)
            pi[i] = m->dims[i][n];

        int cmp = spndarray_compare_idx(ndim, idxs, pi);
        if (cmp < 0)
            p = p->avl_link[0];
        else if (cmp > 0)
            p = p->avl_link[1];
        else
            return p->avl_data;
    }
    return NULL;
}
