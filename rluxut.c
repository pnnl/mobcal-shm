#include "system_includes.h"
#include "mobcal_state_struct.h"
#include "rluxut.h"
int rluxut(struct mobcal_state_struct *state, int *isdext) {
  struct ranlux_state_struct *rstate;
  double *seeds;
  double carry;
  double twop24;
  int success;
  int i24;

  int j24;
  int in24;

  int luxlev;
  int i;

  success = 1;
  twop24  = 16777216;
  rstate = state->ranlux_state;
  for (i=0;i<24;i++) {
    isdest[i] = (int)(seeds[i]*twop24);
  }
  i24 = rstate->i24;
  j24 = rstate->j24;
  in24 = rstate->in24;
  luxlev = rstate->luxlev;
  isdext[24] = i24 + (j24 << 5) + (in24 << 10) + (luxlev << 15);
  if (carry > 0) {
    isdext[24] = -isdext[24];
  }
  return(success);
}
