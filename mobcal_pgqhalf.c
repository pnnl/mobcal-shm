#include "system_includes.h"
#include "mobcal_pgqhalf.h"
void mobcal_pgqhalf(int n, double *p, double *q) {
  /*
    p <- .sqrt(q)
  */
  int i;
  int padi;
  for (i=0;i<8;i++) {
    p[i] = sqrt(q[i]);
  }
}
