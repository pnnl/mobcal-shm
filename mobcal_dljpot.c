#include "system_includes.h"
#include "mobcal_state_struct.h"
#include "mobcal_dljpot_dgt.h"
#include "mobcal_dljpot_iu.h"
#include "mobcal_dljpot.h"
void mobcal_dljpot(struct mobcal_state_struct *state,
		   double x, double y, double z,
		   double *pot_p,
		   double *dpotx_p, 
		   double *dpoty_p,
		   double *dpotz_p,
		   double *dmax_p) {
  /*
    Routine to calucate the L-J + ion-dipole potential.
    Called by: mobcal_deriv
    Calls:     exp, ddot_, mobcal_vec_set, mobcal_dljpot_inner, mobcal_max
    Uses:      inatom, xk, xeo, xe, romax, four_pi, coords,
               properties, (fx, fy, fz,
               eox4, ro6lj, ro12lj, dro6, dro12, and pcharge) fields of state,
	       and mthree_pcharge, and pc_scale_pcharge scalings of pcharge.
    Sets       pot, dpotx, dpoty, dpotz and dmax.	       
  */
  if (state->use_iu_dljpot) {
    mobcal_dljpot_iu(state,x,y,z,pot_p,dpotx_p,dpoty_p,dpotz_p,dmax_p);
  } else {
    mobcal_dljpot_dgt(state,x,y,z,pot_p,dpotx_p,dpoty_p,dpotz_p,dmax_p);
  }
}
