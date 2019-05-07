#include "system_includes.h"
#include "mobcal_ew_add.h"
void mobcal_ew_add(int n_coords,double *a, double *b, double *c) {
  /*
    Perform the elment-wise addition a[i] = b[i] + c[i], i=0;i<n
  */
  int i;
  int padi;
  for (i=0;i<n;i++) {
    a[i] = b[i] + c[i];
  }
}
