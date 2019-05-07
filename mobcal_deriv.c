#include "system_includes.h"
#include "mobcal_state_struct.h"
#include "mobcal_dljpot.h"
#include "mobcal_deriv.h"
void mobcal_deriv(struct mobcal_state_struct *state,
		  double *w, double *dw, double *pot_p,
		  double *dmax_p) {
  /*
    Defines Hamilton's equations of motion as the time derivatives 
    of the coordinates and momenta.
    Called by: mobcal_diffeq, mobcal_gsang
    Calls:     mobcal_dljpot
    Uses mu field of state. Sets dw components, and contents of the
    pointers passed as arguments.
  */
  double mu;
  double recip_mu;
  double x;
  double y;
  double z;
  double pot;
  double dpotx;
  double dpoty;
  double dpotz;
  double dmax;

  mu = state->mu;
  recip_mu = 1.0/mu;
  
  dw[0] = w[1]*recip_mu;
  dw[2] = w[3]*recip_mu;
  dw[4] = w[5]*recip_mu;

  x  = w[0];
  y  = w[2];
  z  = w[4];
  mobcal_dljpot(state,x,y,z,&pot,&dpotx,&dpoty,&dpotz,&dmax);
  dw[1]    = -(dpotx);
  dw[3]    = -(dpoty);
  dw[5]    = -(dpotz);
  *pot_p   = pot;
  *dmax_p  = dmax;
}
