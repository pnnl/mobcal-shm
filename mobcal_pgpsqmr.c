#include "system_includes.h"
#include "mobcal_pgpsqmr.h"
void mobcal_pgpsqmr(int n, double *p, double *q, double *r) {
  /*
    p <- p .- (q .* r)
    This routine should get a fused multipliy add if available.
  */
  int i;
  int padi;
  for (i=0;i<8;i++) {
    p[i] = p[i] - (q[i] * r[i]);
  }
}
