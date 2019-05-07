#include "system_includes.h"
#include "blas.h"
#include "mobcal_state_struct.h"
#include "mobcal_init_constants_0.h"
int mobcal_init_constants_0 (struct mobcal_state_struct *state) {
  /*
    Initialize print switches and physics/math concstants,
    used in mobcal_ncoord.
    Called by: mobcal
    Sets: ip,it,iu1,iu2,iu3,iv,im2,im4,igs printswitch fields,
          and pi, cang, xe, xk, xn, xeo, and xmv fields of state.
    
    print switches ip=1  print scattering angles
                   it=1  print trajectory
                   iu1=1 print initial coordinates from mobcal_ncoord
                   iu2=1 print angles and coordinates from ROTATE 
                   iu3=1 print angles from ROTATE
                   iv=1  print all potentials from POTENT
                   im2=1 print brief information from MOBIL2
                   im4=1 print brief information from MOBIL4
                   igs=1 print out angles, v, and b from MOBIL4 into 
                         a temporary file called hold
  */
  struct mobcal_diffeq_state_struct *diffeq_state;
  double pi; 
  double quarter_pi;
  double *a;
  double *b;
  double *c;
  double *three_a;
  double *ampc;
  double *amcc;
  double three;
  double xe;
  double four_pi;
  double xeo;

  int    incx;
  int    ifour;
  int    success;
  int    padi;

  state->ip   = 0;
  state->it   = 0;
  state->iu1  = 0;
  state->iu2  = 0;
  state->iu3  = 0;
  state->iv   = 0;
  state->im2  = 0;
  state->im4  = 0;
  state->igs  = 0;
  /*
    Constants from Handbook of Chemistry and Physics, 70th Edition
    (except pi which should always be computed as below)
  */
  quarter_pi = atan(1.0);
  pi          = 4.0*quarter_pi;
  state->quarter_pi = quarter_pi;
  state->half_pi    = 2.0 * quarter_pi;
  state->pi         = pi;
  state->two_pi     = 2.0 *pi;
  four_pi           = 4.0 *pi;
  state->four_pi    = four_pi;
  state->eight_pi   = 8.0 *pi;
  state->cang 	    = 180.0/pi;
  xe    	    = 1.60217733e-19;
  state->xe         = xe;
  state->xk   	    = 1.380658e-23;
  state->xn   	    = 6.0221367e+23;
  xeo  	            = 8.854187817e-12;
  state->xeo        = xeo;
  state->pc         = -0.4825;
  state->xmv  	    = 0.02241410;
  state->romax_adj  = 5.5275e-11; /* 1.1055e-10/2 */
  state->xe2_o_4pi_xeo = (xe * xe) / (four_pi * xeo);
  /*
    Initialize the constants used in the mobcal_diffeq routine.
  */
  diffeq_state = state->diffeq_state;
  a            = diffeq_state->a;
  b            = diffeq_state->b;
  c            = diffeq_state->c;
  three_a      = diffeq_state->three_a;
  ampc         = diffeq_state->ampc;
  amcc         = diffeq_state->amcc;
  a[0] = 0.5;
  a[1] = 0.292893218814;
  a[2] = 1.70710678118;
  a[3] = 0.1666666666667;
  b[0] = 2.0;
  b[1] = 1.0;
  b[2] = 1.0;
  b[3] = 2.0;
  c[0] = -0.5;
  c[1] = -0.292893218814;
  c[2] = -1.70710678118;
  c[3] = -0.5;
  ampc[0] = -0.111059153612;
  ampc[1] = 0.672667757774;
  ampc[2] = -1.70633621697;
  ampc[3] = 2.33387888707;
  ampc[4] = -1.8524668225;
  amcc[0] = 0.0189208128941;
  amcc[1] = -0.121233356692;
  amcc[2] = 0.337771548703;
  amcc[3] = -0.55921513665;
  diffeq_state->var     = 2.97013888888;
  diffeq_state->cvar    = 0.990972222222;
  diffeq_state->acst    = 0.332866152768;
  ifour   = 4;
  incx    = 1;
  three   = 3.0;
  /*
    three_a <- 3 * a;
  */
  dcopy_(&ifour,a,&incx,three_a,&incx);
  /*
  To match single precison 3.0 in diffeq we replace the dscal with
  the following.
  */
  /*
  for (i=0;i<3;i++) {
    three_a[i] = fthree * a[i];
  }
  */
  dscal_(&ifour,&three,three_a,&incx);
  success = 1;


  return(success);
}
