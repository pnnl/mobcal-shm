#include "system_includes.h"
#include "mobcal_state_struct.h"
#include "rluxat.h"
int rluxat(struct mobcal_state_struct *state, int *lout,
	   int *inout, int *k1, int *k2) {
  struct ranlux_state_struct *rstate;
  int success;
  int padi;
  success = 1;
  rstate = state->ranlux_state;
  *lout = rstate->luxlev;
  *inout = rstate->inseed;
  *k1    = rstate->kount;
  *k2    = rstate->mkount;
  return(success);
}

