#include "system_includes.h"
#include "mobcal_state_struct.h"
#include "mobcal_init_gst_intgrl.h"
void mobcal_init_gst_intgrl(struct mobcal_state_struct *state) {
  /*
    Set up integration over gst filling wgst and pgst vectors.
    Called by: mobcal_mobil2
    Calls:     fprintf, fflush, sqrt, exp
    Uses xk, temp, eo, mu, im2, inp, and ofp fields of state.
    Sets the wgst and pgst vectors of state.
    mobil2 600:2000
  */
  double *wgst;
  double *pgst;
  double xk;
  double temp;
  double eo;
  double mu;
  double tst;
  double tst3;
  double recip_tst;
  double recip_3_tst;
  double recip_12_tst_2;
  double hold1;
  double hold2;
  double hold3;
  double dgst;
  double half_dgst;
  double gst;
  double gstt;
  double gst2;
  double gst3;
  double gst5;
  double sum;
  double sum1;
  double recip_sum1;
  double sum2;
  double rti;
  double rtim1;
  double half_mu;
  double recip_half_mu;
  double recip_xk_temp;
  double sum_o_tst3;
  double pgsti;
  double pgst2;
  double pgst3;
  double pgst5;
  double wgsti;

  int im2;
  int inp;
  int i;
  int thread_id;
  FILE *ofp;
  FILE *efp;
  im2   = state->im2;
  inp   = state->inp;
  xk    = state->xk;
  eo    = state->eo;
  mu    = state->mu;
  temp  = state->temp;
  ofp   = state->ofp;
  wgst  = state->wgst;
  pgst  = state->pgst;
  thread_id = state->thread_id;
  tst      = xk * temp / eo;
  if (im2 == 0) {
    if (thread_id == 0) {
      if (ofp) {
	fprintf(ofp," t*= %le\ntemp = %le\n",tst,temp);
      }
    }
  }
  tst3 = tst * tst * tst;
  recip_tst = 1.0/tst;
  recip_3_tst = recip_tst/3.0;
  recip_12_tst_2 = (recip_tst * recip_tst) / 12.0;
  state->tst = tst;
  state->tst3 = tst3;
  state->recip_3_tst = recip_3_tst;
  state->recip_12_tst_2 = recip_12_tst_2;
  
  dgst  = 0.000003*sqrt(tst);
  half_dgst  = dgst * 0.5;
  gst   = dgst;
  sum1  = 0.0;
  for (i=1;i<=inp;i++) {
    sum1 += sqrt((double)i);
  }
  if (im2 == 0) {
    if (thread_id == 0) {
      if (ofp) {
	fprintf(ofp," setup gst integeration - integration over velocity\n");
	fprintf(ofp,"     pgst        wgst         v         ke/kt       gst^5*     frac of\n");
	fprintf(ofp,"                                              exp(gst*2/tst)     sum\n");
      }
    }
  }
  rti = 0.0;
  sum   = 0.0;
  sum2  = 0.0;
  recip_sum1 = 1.0/sum1;
  recip_half_mu = 2.0/mu;
  half_mu       = 0.5*mu;
  recip_xk_temp = 1.0/(xk * temp);
  for (i=0;i<inp;i++) {
    rtim1   = rti;
    rti     = sqrt((double)(i+1));
    hold1   = rti;
    hold2   = rtim1;
    sum2    += hold2;
    wgsti   = hold1 * recip_sum1;
    wgst[i] = wgsti;
    gstt    =  tst3 * (sum2 + (hold1 * 0.5)) * recip_sum1;
    gst2 = gst * gst;
    gst3 = gst2 * gst;
    gst5 = gst2 * gst3;
    sum += (exp( - gst2 * recip_tst) * gst5 * dgst);
    gst += dgst;
      
    while (sum < gstt) {
      gst2 = gst * gst;
      gst3 = gst2 * gst;
      gst5 = gst2 * gst3;
      sum += (exp(-gst2*recip_tst) * gst5 * dgst);
      gst += dgst;
    }
    pgsti = gst - half_dgst;
    pgst[i] = pgsti;

    hold2 = pgsti * pgsti * eo;
    hold1 = hold2 * recip_half_mu;
    hold2 = hold2 * recip_xk_temp;
    /*
    hold1 = pgsti * pgsti * eo * recip_half_mu;
    hold2 = half_mu * hold1 * recip_xk_temp;
    */
    hold1 = sqrt(hold1);
    pgst2 = pgsti * pgsti;
    pgst3 = pgst2 * pgsti;
    pgst5 = pgst2 * pgst3;
    hold3 = exp(-pgst2 * recip_tst) * pgst5;
    if (im2 == 0) {
      if (thread_id == 0) {
	if (ofp) {
	  sum_o_tst3 = sum/tst3;
	  fprintf(ofp,"%le %le %le %le %le %le\n",pgsti,wgsti,hold1,hold2,hold3,sum_o_tst3);
	}
      }
    }
  } /* end for (i...) */
  if (ofp) {
    fflush(ofp);
  }
}
