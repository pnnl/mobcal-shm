#include "system_includes.h"
#include "mobcal_state_struct.h"
#include "mobcal_mobil2_result.h"
void mobcal_mobil2_result(struct mobcal_state_struct *state,
			  double *omega_dist,
			  double *mob_p,
			  double *cs_p,
			  double *sdevpc_p) {
  /* 
    Use omegas to obtain higher order correction factor to mobility.
    Called by: mobcal_mobil2
    Calls:     sqrt
  */
  double m1;
  double m2;
  double mom11st;
  double mom12st;
  double mom13st;
  double mom22st;
  double sdom11st;
  double sterr;
  double mob;
  double cs;
  double sdevpc;
  double pi_ro2;
  double ayst;
  double best;
  double cest;
  double term;
  double u2;
  double w;
  double hold1;
  double delta;
  double f;
  double mob_inverse;
  double tm_cross_sctn;
  double mconst;
  double temp;
  double ten_p20;
  double stdv_tm_cross;
  
  int im2;
  FILE *ofp;

  m1       = state->m1;
  m2       = state->m2;
  pi_ro2   = state->pi_ro2;
  im2      = state->im2;
  ofp      = state->ofp;
  mconst   = state->mconst;
  temp     = state->temp;
  mom11st  = omega_dist[0];
  mom12st  = omega_dist[1];
  mom13st  = omega_dist[2];
  mom22st  = omega_dist[3];
  sdom11st = omega_dist[4];
  sterr    = omega_dist[5];
  ten_p20  = 1.0e20;
  cs = mom11st * pi_ro2;
  sdevpc = 100.0 * sdom11st/mom11st;
  ayst = mom22st/mom11st;
  best = ((5.0*mom12st) - (4.0*mom13st))/mom11st;
  cest = mom12st/mom11st;
  term = ((4.0*ayst)/15.0) + (0.5 * (m2-m1) * (m2 - m1)/(m1*m2));
  u2   = term-(.08333*((2.4*best)+1.0)*(m1/m2));
  w    = m1/m2;
  hold1 = (6.0*cest)-5.0;
  hold1 = hold1 * hold1;
  delta = (hold1 * w)/(60.0*(1+u2));
  f     = 1.0/(1.0-delta);
  if (im2 == 0) {
    fprintf(ofp,"\nf value for second order correction = %le\n",f);
    fprintf(ofp,"(integeratios for second order correction are not\n");
    fprintf(ofp,"accurate, check if correction becomes significant)\n\n");
    fprintf(ofp,"omega*12 = %le\nomega*13 = %le\nomega*22 = %le\n",
	    mom12st,mom13st,mom22st);
    fprintf(ofp,"      u2 = %le\nm1overm2 = %le\ndelta = %le\n\n",
	    u2,w,delta);
  }
  mob = mconst*f/(sqrt(temp)*cs);
  mob_inverse = 1.0/mob;
  tm_cross_sctn = cs*ten_p20;
  stdv_tm_cross = sdom11st * pi_ro2 * ten_p20;
  fprintf(ofp,"\naverage (second order) TM mobility = %le\n",mob);
  fprintf(ofp,"inverse average (second order) TM mobility = %le\n",
	  mob_inverse);
  fprintf(ofp,"average TM cross section = %le\n\n",tm_cross_sctn);
  fprintf(ofp,"standard deviation TM cross section = %le\n\n",stdv_tm_cross);
  *mob_p = mob;
  *cs_p  = cs;
  *sdevpc_p = sdevpc;
}
