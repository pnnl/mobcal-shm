#include "system_includes.h"
#include "mobcal_state_struct.h"
#include "xrand.h"
double xrand(struct mobcal_state_struct *state) {
  struct ranlux_state_struct *rstate;
  double value;
  int    one;
  int    padi;
  rstate = state->ranlux_state;
  state->i6 += 1;
  one = 1;
  if (state->i5 == 1) {
    success = ranlux(state,&value,one);
  } else {
    value = rand();
  }
  return(value);
}
