#include "system_includes.h"
#include "mobcal_pgq2.h"
void mobcal_pgq2(int ivec_len, double *p, double *q) {
  /*
    p <- q .* q
  */
  int i;
  int padi;
  for (i=0;i<8;i++) {
    p[i] = q[i] * q[i];
  }
}
