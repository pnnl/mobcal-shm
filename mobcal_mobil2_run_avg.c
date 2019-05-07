#include "system_includes.h"
#include "mobcal_state_struct.h"
#include "mobcal_mobil2_run_avg.h"
void mobcal_mobil2_run_avg(struct mobcal_state_struct *state) {
  /*
    Calculate running averages. Only called if im2 == 0
    Called by: mobcal_mobil2
    Calls:     fprintf, fflush
    Uses itn, inp, im2, ofp, pi_ro2, temp, recip_mconst, 
         om11st, pgst, and q1st fields
    of state.
    
  */
  double pi_ro2;
  double pi_ro2_ten_p20;
  double temp;
  double recip_mconst;
  double *om11st;
  double *pgst;
  double *q1st;
  double *wgst;
  double hold1;
  double hold2;
  double cs_o_a2;
  double avg_cs_o_a2;
  double avg_recip_ko;
  double recip_ko;
  double recip_icc;
  double scale1;
  double pgst_ig;
  double pgst_ig2;
  double q1st_avg;
  double om11st_icc;
  double recip_inp;
  double scaled_q1st;
  int itn;
  int icc;
  int ig;
  int inp;
  FILE *ofp;
  FILE *efp;
  pi_ro2 = state->pi_ro2;
  temp   = state->temp;
  recip_mconst = state->recip_mconst;
  om11st       = state->om11st;
  pgst         = state->pgst;
  q1st         = state->q1st;
  itn          = state->itn;
  inp          = state->inp;
  ofp          = state->ofp;
  wgst         = state->wgst;
  pi_ro2_ten_p20 = pi_ro2 * 1.0e20;
  recip_inp    = 1.0/((double)inp);
  fprintf(ofp," summary of mobility calculations\n");
  fprintf(ofp," cycle     cs/A^2      avge cs/A^2        Ko^-1"
	  "       avge Ko^-1\n");
  hold1 = 0.0;
  hold2 = 0.0;
  scale1 = sqrt(temp) * pi_ro2 * recip_mconst;
  for (icc=0;icc<itn;icc+=1) {
    om11st_icc = om11st[icc];
    recip_icc = 1.0/((double)(icc+1));
    recip_ko =  om11st_icc*scale1;
    hold1 += om11st_icc;
    /*
      hold3 = 1.0/(mconst/(sqrt(temp)*om11st[icc]*pi*ro*ro));
      hold2 += hold3;
    */
    hold2 = hold1 * scale1;
    cs_o_a2 = om11st_icc * pi_ro2_ten_p20;
    avg_cs_o_a2 = hold1 * pi_ro2_ten_p20 * recip_icc;
    avg_recip_ko = hold2 * recip_icc;
    
    fprintf(ofp,"%d    %le    %le    %le    %le\n",
	    icc, cs_o_a2, avg_cs_o_a2, recip_ko, avg_recip_ko);
    
  }
  fflush(ofp);
  fprintf(ofp,"average values for q1st\n");
  fprintf(ofp,"     gst2        wgst        q1st_avg\n");
  for (ig=0;ig<inp;ig++) {
    pgst_ig = pgst[ig];
    pgst_ig2 = pgst_ig*pgst_ig;
    /*
    q1st_avg = q1st[ig] / ((double)(ig+1));
    */
    q1st_avg = q1st[ig] * recip_inp;
    fprintf(ofp,"%le %le %le\n",pgst_ig2,wgst[ig],q1st_avg);
  }
  fprintf(ofp,"\n");
  fflush(ofp);
}
