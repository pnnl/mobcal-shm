#include "system_includes.h"
#include "mobcal_vec_set.h"
void mobcal_vec_set(int n, double *v, double value) {
  int i;
  int padi;
  for (i=0;i<n;i++) {
    v[i] = value;
  }
}
