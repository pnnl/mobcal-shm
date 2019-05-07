#ifndef _MOBCAL_DIFFEQ_H_
#define _MOBCAL_DIFFEQ_H_ 1
extern void mobcal_diffeq(struct mobcal_state_struct *state,
			  int *l_p,
			  double *tim_p,
			  double *dt_p,
			  double *w,
			  double *dw,
			  double *pot_p,
			  double *dmax_p);
#endif
