#include "spndarray.h"
#include <stdio.h>

/*
 * spndarray_fwrite()
 *
 * writes the sparse array to the given file,
 * applying the given format to the elements
 *
 * Inputs
 *  fmt - element format
 *  filepath - file to save to
 *  full - whether to save only nonzero elements (1) or all the elements (0)
 */
int spndarray_fwrite(const spndarray* m, const char* fmt, const char* filepath, const int full) {
  FILE* fp = fopen(filepath, "w+");
  if (!fmt) fmt = "%f,\n";
  if (!full) {
    size_t ndim = m->ndim;
    for (size_t s = 0; s < m->nz; s++) {
      fprintf(fp, "[");
      for (size_t j = 0; j < ndim; j++)
        fprintf(fp, (ndim-1 == j ? "%ld] = " : "%ld,"), m->dims[j][s]);
      fprintf(fp, fmt, m->data[s]);
    }
    fclose(fp);
    return 0;
  } else {
    fprintf(stderr, "nonsparse saving not implemented yet\n");
    return 1;
  }
}
