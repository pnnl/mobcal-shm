#include "system_includes.h"
#include "mobcal_state_struct.h"
#include "genrand64_real3.h"
#include "mobcal_ranlux.h"
int mobcal_gen_rn(struct mobcal_state_struct *state,
		  double *rvec, int lenv) {
  /*
    Generate a list of random numbers. Assumes mobcal_init_rn has been
    called to initialize the appropriate random number generator.
    Called by: mobcal_init_concstants_1
    Calls:     genrand64_real3, mobcal_ranlux, or rand.
  */
  struct mt_state_struct *mt_state;
  int use_mt;
  int i5;
  int i;
  int success;
  success = 1;
  i5 = state->i5;
  if (i5 != 1) {
    for (i=0;i<lenv;i++) {
      rvec[i] = rand();
    }
  } else {
    use_mt = state->use_mt;
    if (use_mt) {
      mt_state = state->mt_state;
      for (i=0;i<lenv;i++) {
	rvec[i] = genrand64_real3(mt_state);
      }
    } else {
      success = mobcal_ranlux(state,rvec,lenv);
    }
  }
  return(success);
}

