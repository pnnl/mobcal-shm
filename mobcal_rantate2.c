#include "system_includes.h"
#include "mobcal_state_struct.h"
#include "mobcal_rotate.h"
#include "mobcal_rantate2.h"
void mobcal_rantate2 (struct mobcal_state_struct *state, int iran) {
  /*
    Rotate the cluster/molecule to a random orientation.
    Called by: mobcal_mobil2
    Calls:     mobcal_xrand, mobcal_rotate, asin
  */
  double *ranlist;
  double rnt;
  double rnp;
  double rng;
  double theta;
  double phi;
  double gamma;
  double two_pi;
  double half_pi;
  FILE   *ofp;
  FILE   *efp;
  ranlist = state->ranlist;
  /*
  rnt     = mobcal_xrand(state);
  rnp     = mobcal_xrand(state);
  rng     = mobcal_xrand(state);
  */
  rnt     = ranlist[iran];
  rnp     = ranlist[iran+1];
  rng     = ranlist[iran+2];
#ifdef DBG
  fprintf(stdout,"rnt = %le, rnp = %le, rng = %le\n",rnt,rnp,rng);
  fflush(stdout);
#endif
  two_pi  = state->two_pi;
  half_pi = state->half_pi;
  theta = rnt * two_pi;
  gamma = rng * two_pi;
  phi   = asin((rnp*2.0) - 1.0) + half_pi;
  state->theta = theta;
  state->gamma = gamma;
  state->phi   = phi;
  mobcal_rotate(state,theta,phi,gamma);
}
