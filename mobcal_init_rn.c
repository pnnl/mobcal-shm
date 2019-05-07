#include "system_includes.h"
#include "mobcal_state_struct.h"
#include "init_genrand64.h"
#include "mobcal_rluxgo.h"
#include "mobcal_init_rn.h"
int mobcal_init_rn(struct mobcal_state_struct *state) {
  /*
    Initialize a randomnumber generator for the mobcal code.
    Called by: mobcal_init_constants_1
    Calls:     srand, init_genrand64,mobcal_rluxgo
    Uses i1,i2,i3,i4,i5, and use_mt fields of state.
  */
  struct mt_state_struct *mt_state;
  unsigned long long ull_i2;
  int i1;
  int i2; 
  int i3; 
  int i4; 
  int i5;
  int use_mt;
  int success;
  int padi;
  success = 1;
  i1 = state->i1;
  i2 = state->i2;
  i3 = state->i3;
  i4 = state->i4;
  i5 = state->i5;
  use_mt = state->use_mt;
  if (i5 != 1) {
    srand(i2);
  } else { 
    if (use_mt) {
      mt_state = state->mt_state;
      ull_i2 = (unsigned long long)i2;
      init_genrand64(mt_state,ull_i2);
    } else {
      success = mobcal_rluxgo(state,i1,i2,i3,i4);
    }
  }
  return (success);
}
