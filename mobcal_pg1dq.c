#include "system_includes.h"
#include "mobcal_pg1dq.h"
void mobcal_pg1dq(int n, double *p, double *q) {
  /*
    p <- 1 ./ q
  */
  int i;
  int padi;
  for (i=0;i<8;i++) {
    p[i] = 1.0/q[i];
  }
}
