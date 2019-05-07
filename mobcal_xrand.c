#include "system_includes.h"
#include "mobcal_state_struct.h"
#include "mobcal_ranlux.h"
#include "mobcal_xrand.h"
double mobcal_xrand(struct mobcal_state_struct *state) {
  struct ranlux_state_struct *rstate;
  /*
    Extract a single random number from the random number generator.
    Called by: mobcal_rantate
    Calls:     mobcal_ranlux,rand
  */
  double value;
  int    one;
  int    success;
  rstate = state->ranlux_state;
  state->i6 += 1;
  one = 1;
  if (state->i5 == 1) {
    success = mobcal_ranlux(state,&value,one);
  } else {
    value = rand();
  }
  return(value);
}
