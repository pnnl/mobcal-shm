#include "system_includes.h"
#include "mobcal_state_struct.h"
#include "mobcal_amax.h"
double mobcal_amax(double *p, double prev_max) {
  double pmax;
  int i;
  int padi;
  pmax = prev_max;
  for (i=0;i<_IVEC_LEN_;i++) {
    if (p[i] > pmax) {
      pmax = p[i];
    }
  }
  return(pmax);
}
