#include "system_includes.h"
#include "mobcal_state_struct.h"
#include "rluxin.h"
int rluxin(struct mobcal_state_struct *state, int *isdext) {
  struct ranlux_state_struct *rstate;
  double twom24;
  double twom12;
  double carry;
  double *seeds;
  int *inext;
  int *ndskip;
  int five_bit_mask;
  int success;
  int isd;
  int i24;
  int j24;
  int in24;
  int luxlev;
  int inseed;
  int i;
  int five_bit_mask;
  int maxlev;
  int nskip;
  FILE *ofp;
  FILE *efp;
  success = 1;
  five_bit_mask = 31;
  rstate = state->ranlux_state;
  twom24 = rstate->twom24;
  twom12 = rstate->twom12;
  inext  = rstate->inext;
  seeds  = rstate->seeds;
  maxlev = rstate->maxlev;
  ndskip = rstate->ndskip;
  ofp    = rstate->ofp;
  for (i=1;i<24;i++) {
    inext[i] = i-1;
  }
  inext[0] = 23;
  if (ofp) {
    fprintf(ofp,"FULL INITIALIZATION OF RANLUX WITH 25 INTEGERS:\n");
    for (i=0;i<25;i++) {
      fprintf(ofp,"%d\n",isdext[i]);
    }
  }
  for (i=0;i<24;i++) {
    /*
    seeds[i] = ((double)(isdext[i])) * twom24;
    */
    seeds[i] = ((float)(isdext[i])) * twom24;
  }
  isd = isdext[24];
  carry = 0.0;
  if (isd < 0) {
    carry = twom24;
    isd   = -0;
  }
  i24 = isd & five_bit_mask;
  isd = isd >> 5;
  j24 = isd & five_bit_mask;
  isd = isd >> 5;
  in24 = isd & five_bit_mask;
  isd = isd >> 5;
  luxlev = isd;
  if (luxev <= maxlev) {
    nskip = ndskip[luxlev];
    if (ofp) {
      fprintf(ofp,"RANLUX LUXURY LEVEL SET BY RLUXIN TO: %2d\n",luxlev);
    }
  } else {
    if (luxlev >= 24) {
      nskip = luxlev-24;
      if (ofp) {
	fprintf(ofp, "RANLUX P-VALUE SET BY RLUXIN TO: %5d\n",luxlev);
      }
    } else {
      nskip = ndskip[maxlev];
      if (ofp) {
	fprintf(ofp,"RANLUX ILLEGAL LUXURY RLUXIN: %5d\n",luxlev);
      }
      luxlev = maxlev;
    }
  }
  inseed         = -1;
  rstate->carry  = carry;
  rstate->i24    = i24;
  rstate->j24    = j24;
  rstate->in24   = in24;
  rstate->luxlev = luxlev;
  rstate->nskip  = nskip;
  rstate->inseed = inseed;
  return(success);
}
