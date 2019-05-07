#include "system_includes.h"
#include "mobcal_max.h"
double mobcal_max(int n, double *p) {
  double pmax;
  int i;
  int padi;
  pmax = p[0];
  for (i=1;i<n;i++) {
    if (p[i] > pmax) {
      pmax = p[i];
    }
  }
  return(pmax);
}
