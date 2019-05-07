#include "system_includes.h"
#include "mobcal_state_struct.h"
#include "mobcal_dljpot_only_iu.h"
#include "mobcal_dljpot_only_dgt.h"
#include "mobcal_dljpot_only.h"
void mobcal_dljpot_only(struct mobcal_state_struct *state,
			double x, double y, double z,
			double *pot_p) {
  /*
    Routine to calucate the L-J + ion-dipole potential.
    Called by: mobcal_rmax_emax_r00, mobcal_gsang
    Calls:     mobcal_dljpot_only_iu, mobcal_dljpot_only_dgt
    Uses:      inatom, xk, xeo, xe, four_pi, fx, fy, fz,
               eox4, ro6lj, ro12lj, and pcharge fields of state,
	       and pc_scale_pcharge scaling of pcharge.
    Sets       pot
  */
  if (state->use_iu_dljpot) {
    mobcal_dljpot_only_iu(state,x,y,z,pot_p);
  } else {
    mobcal_dljpot_only_dgt(state,x,y,z,pot_p);
  }
}
