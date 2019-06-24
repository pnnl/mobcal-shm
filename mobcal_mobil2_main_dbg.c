#include "system_includes.h"
#include "mobcal_state_struct.h"
/*
#include "mobcal_xrand.h"
#include "mobcal_rantate.h"
*/
#include "mobcal_rantate2.h"
#include "mobcal_gsang.h"
#include "mobcal_mobil2_main_dbg.h"
void mobcal_mobil2_main_dbg(struct mobcal_state_struct *state) {
  /*
    Run the main mobcal loops, over ic[1:itn]  ig[1:inp], and im[1:imp]
    sets the om11st, om12st, om13st and om22st, q1st, w2st, fields of state.
    Assuming ip = 1, or igs = 1.
    Called by: mobcal_mobil2
    Calls:     fprintf, fflsuh, mobcal_rantate, mobcal_gsang
  */
  double *om11st;
  double *om12st;
  double *om13st;
  double *om22st;
  double *q1st;
  double *q2st;
  double *b2max;
  double *pgst;
  double *wgst;
  double *ranlist;
  double pgst_ig;
  double wgst_ig;
  double pgst_ig2;
  double pgst_ig4;
  double gst2;
  double v;
  double eo;
  double b;
  double ro;
  double erat;
  double ang;
  double d1;
  double temp1;
  double temp2;
  double rnb;
  double bst2;
  double hold1;
  double sin_ang;
  double hold2;
  double recip_imp;
  double recip_3_tst;
  double recip_12_tst_2;
  double recip_half_mu;
  double ten_p10;
  double bau;
  double ang_deg;
  double theta_deg;
  double phi_deg;
  double gamma_deg;
  double cang;
  double temp1d;
  double temp2d;

  int im2;
  int ip;

  int inp;
  int itn;

  int imp;
  int igs;

  int ig;
  int ic;

  int im;
  int tot_num_points;

  int istep;
  int iic;

  int iran;
  int itn_inp_imp_4;
  int inp_imp_4;
  int imp_4;

  int iranp1;
  int padi;


  FILE *ofp;
  FILE *hfp;

  /*
  */
#define DBG 1
  
  im2 = state->im2;
  ip  = state->ip;
  inp = state->inp;
  itn = state->itn;
  imp = state->imp;
  igs = state->igs;
  eo  = state->eo;
  ro  = state->ro;
  ofp = state->ofp;
  q1st = state->q1st;
  q2st = state->q2st;
  om11st = state->om11st;
  om12st = state->om12st;
  om13st = state->om13st;
  om22st = state->om22st;
  b2max  = state->b2max;
  pgst   = state->pgst;
  wgst   = state->wgst;
  cang   = state->cang;
  recip_3_tst    = state->recip_3_tst;
  recip_12_tst_2 = state->recip_12_tst_2;
  recip_half_mu  = state->recip_half_mu;
  itn_inp_imp_4  = state->itn_inp_imp_4;
  inp_imp_4      = state->inp_imp_4;
  imp_4          = state->imp_4;
  ranlist        = state->ranlist;
  recip_imp = 1.0/imp;
  ten_p10 = 1.0e10;
  if (im2 == 0) {
    tot_num_points = itn * inp * imp;
    fprintf(ofp,"number of complete cycles (itn) = %d\n",itn);
    fprintf(ofp,"number of velocity points (inp) = %d\n",inp);
    fprintf(ofp,"number of random points (imp) = %d\n",imp);
    fprintf(ofp,"total number of random points = %d\n",tot_num_points);
    fflush(ofp);
  }
  fprintf(ofp,"start mobility calculation\n");
  for (ig=0;ig<inp;ig++) {
    q1st[ig] = 0;
    q2st[ig] = 0;
  }
  for (ic =0;ic <itn;ic++) {
    om11st[ic] = 0.0;
    om12st[ic] = 0.0;
    om13st[ic] = 0.0;
    om22st[ic] = 0.0;
  }    
  hfp = NULL;
  if (igs == 1) {
    hfp = fopen("hold","w");
    if (hfp == NULL) {
      fprintf(ofp,"mobcal_mobil2: Error unable to open hold file\n");
    }
  }
  /*
    main loop.
  */
  iic = state->iic;
  for (ig = 0; ig < inp; ig++) {
    pgst_ig = pgst[ig];
    wgst_ig = wgst[ig];
    pgst_ig2 = pgst_ig * pgst_ig;
    pgst_ig4 = pgst_ig2 * pgst_ig2;
    gst2 = pgst_ig * pgst_ig;
    v    = pgst_ig * sqrt(eo*recip_half_mu);
    state->v = v;
    for (ic = 0;ic < itn;ic++) {
      if (ip == 1) {
	fprintf(ofp,"cycle_number, ic = %d\n",ic);
	fflush(ofp);
      }
      fprintf(ofp," ic = %d ig = %d gst2 = %le v = %le\n",ic,ig,gst2,v);
      temp1 = 0.0;
      temp2 = 0.0;
      fprintf(ofp,"     b/A        ang      (1-cosX)    e ratio"
		  "    theta       phi       gamma\n");
      fflush(ofp);
      for (im = 0; im < imp; im++) {
	/*
	rnb = mobcal_xrand(state);
	mobcal_rantate(state);
	*/
	/*
	iran = (iic * itn_inp_imp_4) + (ic * inp_imp_4) + (ig *  imp_4) +
	  (im * 4) + 1;
	*/
	iran = (iic * itn_inp_imp_4) + (ic * inp_imp_4) + (ig *  imp_4) +
	  (im * 4) + 4;
	rnb = ranlist[iran];
	iranp1 = iran + 1;
	mobcal_rantate2(state,iranp1);

	bst2 = rnb * b2max[ig];
	b    = ro *sqrt(bst2);
	state->b = b;
	theta_deg = state->theta * cang;
	phi_deg   = state->phi * cang;
	gamma_deg  = state->gamma * cang;
	if (igs == 1) {
	  fseek(hfp,(long)0,SEEK_SET);
	  fprintf(hfp,"%d %d %d %d\n",state->iic,ic,ig,im);
	  fprintf(hfp,"%le %le\n",state->v,b);
	  fprintf(hfp,"%le %le %le\n",theta_deg,phi_deg,gamma_deg);
	  fflush(hfp);
	}
	mobcal_gsang(state,v,b,&erat,&ang,&d1,&istep);
	hold1 = 1.0 - cos(ang);
	sin_ang = sin(ang);
	hold2 = sin_ang * sin_ang;
	bau = b*ten_p10;
	ang_deg = ang * cang;
	fprintf(ofp," %le %le %le %le %le %le %le\n",
		bau,ang_deg,hold1,erat,theta_deg,phi_deg,gamma_deg);
	temp1d = ((hold1 * b2max[ig])*recip_imp);
	temp1 += temp1d;
	temp2d = (((1.5 * hold2) * b2max[ig])*recip_imp);
	temp2 += temp2d;
#ifdef DBG
	if (ofp) {
	  fprintf(ofp,"mobcal_mobil2_main: ic = %d, ig = %d, im = %d, ang = %le,temp1d = %le, temp2d = %le\n",ic,ig,im,ang,temp1d,temp2d);
	  fflush(ofp);
	}
#endif

      } /* end for (im...) */
      /*
	NB the following constants after temp1/temp2 could be precomputed.
	They don't depend on ic at all.
      */

      om11st[ic] += (temp1 * wgst_ig);
      om12st[ic] += (temp1 * pgst_ig2 * wgst_ig * recip_3_tst);
      om13st[ic] += (temp1 * pgst_ig4 * wgst_ig * recip_12_tst_2);
      om22st[ic] += (temp2 * pgst_ig2 * wgst_ig * recip_3_tst);
      q1st[ig]   += temp1;
      q2st[ig]   += temp2;
    } /* end for (ig...) */
    fprintf(ofp,"OMEGA(1,1) = %le\n",om11st[ic]);
  } /* end for (ic...) */
}

