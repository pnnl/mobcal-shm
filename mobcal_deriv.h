#ifndef _MOBCAL_DERIV_H_
#define _MOBCAL_DERIV_H_ 1
extern void mobcal_deriv(struct mobcal_state_struct *state,
			 double *w, double *dw, double *pot_p,
			 double *dmax_p);
#endif
