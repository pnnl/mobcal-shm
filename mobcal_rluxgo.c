#include "system_includes.h"
#include "mobcal_state_struct.h"
#include "mobcal_rluxgo.h"
int mobcal_rluxgo(struct mobcal_state_struct *state,
		  int lux, int ins, int k1, int k2) {
  /*
    C adaptation of fortran, rluxgo.
    Initialize the rlux_state_struct;
    Called by: mobcal_init_constants_1
    Calls:     fprintf, fflush;
  */
  struct ranlux_state_struct *rstate;
  double *seeds;
  double uni;
  double twom12;
  double twom24;
  double twom48;
  double carry;

  int *ndskip;
  int *inext;

  int luxlev;
  int maxlev;

  int ilx;
  int nskip;

  int p;
  int jseed;

  int is1;
  int is2;

  int is3;
  int i;

  int k;
  int icons;

  int iseedi;
  int mod_2_to_24_mask;

  int i24;
  int j24;

  int kount;
  int mkount;

  int in24;
  int nskip_p_24;

  int izip;
  int izip2;

  int igiga;
  int inner;

  int iouter;
  int isk;

  int success;
  

  FILE *ofp;
  FILE *efp;
  success           = 1;
  in24              = 0;
  ofp               = state->ofp;
  rstate            = state->ranlux_state;
  ndskip            = rstate->ndskip;
  inext             = rstate->inext;
  seeds             = rstate->seeds;
  rstate->uninitialized = 0;
  rstate->carry     = 0.0;
  rstate->maxlev    = 4;
  rstate->lxdflt    = 3;
  rstate->igiga     = 1000000000;
  rstate->jsdflt    = 314159265;
  /*
    As icons is used to make a negative numer positive, it really should 
    be the ...647, but we need to preserve the original random number
    gnerator to duplicate it exactly.
  rstate->icons     = 2147483647; 
  */
  rstate->icons     = 2147483563; 
  rstate->luxlev    = rstate->lxdflt;
  rstate->in24      = 0;
  rstate->kount     = (int64_t)0;
  rstate->mkount    = (int64_t)0;
  rstate->i24       = 23;
  rstate->j24       = 9;
  ndskip[0]         = 0;
  ndskip[1]         = 24;
  ndskip[2]         = 73;
  ndskip[3]         = 199;
  ndskip[4]         = 365;
  maxlev            = rstate->maxlev;
  inext             = rstate->inext;
  icons             = rstate->icons;
  igiga             = (int64_t)rstate->igiga;
  mod_2_to_24_mask  = 16777215;
  if (lux < 0) {
    luxlev = rstate->lxdflt;
  } else {
    if (lux <= maxlev) {
      luxlev = lux;
    } else { 
      if ((lux < 24 ) || ( lux > 2000)) {
	luxlev = rstate->maxlev;
	if (ofp) {
	  fprintf(ofp,"mobcal_ranlux: illegal luxury in mobcal_rluxgo: %d\n",lux);
	  fflush(ofp);
	}
      } else {
	luxlev = lux;
	for (ilx=0;ilx<=maxlev;ilx++) {
	  if (lux == (ndskip[ilx] + 24)) {
	    luxlev = ilx;
	  }
	}
      }
    }
  }
  if (luxlev <= maxlev) {
    nskip = ndskip[luxlev];
    p     = nskip + 24;
    if (ofp) {
      fprintf(ofp,"mobcal_ranlux luxury level set by mobcal_rluxgo :%2d     p=%4d\n\n",
	      luxlev,p);
      fflush(ofp);
    }
  } else {
    nskip = luxlev - 24;
    if (ofp) {
      fprintf(ofp,"mobcal_ranlux p-value set by mobcal_rluxgo to:%5d\n\n",luxlev);
      fflush(ofp);
    }
  }
  rstate->luxlev = luxlev;
  rstate->nskip  = nskip;
  rstate->in24   =  0;
  if (ins > 0) {
    jseed = ins;
  } else {
    if (ofp) {
      fprintf(ofp,"mobcal_rluxgo: Illegal initialization, negative input seed: using default\n");
    }
    fprintf(stderr,"rluxgo: Illegal initialization, negative input seed: using default\n");
    jseed = rstate->jsdflt;
  }
  if (ofp) {
    fprintf(ofp,"mobcal_ranlux initialized by mobcal_rluxgo from seeds %2d %2d %2d\n",
	    jseed,k1,k2);
  }
  twom24 = 1.0/16777216;
  twom12 = 1.0/4096.0;
  twom48 = twom24*twom24;
  rstate->twom24 = twom24;
  rstate->twom12 = twom12;
  rstate->twom48 = twom48;
  is1 = 53668;
  is2 = 40014;
  is3 = 12211;
  for (i=0;i<24;i++) {
    k = jseed/is1;
    jseed = is2 * (jseed - (k*is1)) - (k*is3);
    if (jseed < 0) {
      jseed += icons;
    }
    iseedi = jseed & mod_2_to_24_mask;
    seeds[i] = ((double)iseedi)*twom24;
    inext[i] = i-1;
  }
  inext[0] = 23;
  i24      = 23;
  j24      =  9;
  carry    = 0.0;
  if (seeds[23] == 0.0) {
    carry = twom24;
  }
  kount = k1;
  mkount = k2;
  rstate->kount = (int64_t)kount;
  rstate->mkount = (int64_t)mkount;
  if ((k1 + k2) != 0) {
    for (iouter=0;iouter<=k2;iouter+=1) {
      inner = igiga;
      if (iouter == k2) inner = k1;
      for (isk=0;isk<inner;isk+=1) {
	uni = seeds[j24] - seeds[i24] - carry;
	if (uni < 0) {
	  uni = uni + 1.0;
	  carry = twom24;
	} else {
	  carry = 0.0;
	}
	seeds[i24] = uni;
	i24 = inext[i24];
	j24 = inext[j24];
      } /* end for (isk...) */
    }
    nskip_p_24 = (int64_t)nskip + (int64_t)24;
    /*
    in24 = mod(kount,nskip+24)p
    */
    in24 = kount/nskip_p_24;
    in24 = kount - (in24 * nskip_p_24);
    if (mkount > (int64_t)0 ) {
      izip = igiga/nskip_p_24;
      izip = igiga - (izip * nskip_p_24);
      izip2 = mkount * izip + in24;
      in24  = izip2/nskip_p_24;
      in24  = izip2 - (in24 * nskip_p_24);
    }
    if (in24 > 23) {
      if (ofp) {
	fprintf(ofp,"Error in restarting with rmobcal_rluxgo:\n"
		"  The values %d, %d, cannot occur at luxury level %d\n",
		k1,k2,luxlev);
      }
      in24 = 0;
    }
  }
  rstate->i24   = i24;
  rstate->j24   = j24;
  rstate->carry = carry;
  rstate->in24  = (int)in24;
  if (ofp) {
    fflush(ofp);
  }
  return(success);
}
