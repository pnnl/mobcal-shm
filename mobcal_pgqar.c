#include "system_includes.h"
#include "mobcal_pgqar.h"
void mobcal_pgqar(int ivec_len, double *p, double *q, double *r) {
  /*
    p <- q .+ r 
  */
  int i;
  int padi;
  for (i=0;i<8;i++) {
    p[i] = q[i] + r[i];
  }
}
