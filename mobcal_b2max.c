#include "system_includes.h"
#include "mobcal_state_struct.h"
#include "mobcal_gsang.h"
#include "mobcal_b2max.h"
int mobcal_b2max(struct mobcal_state_struct *state, double rmaxx) {
  /*
    Determine b2max vector. Also sets the cosx vector,  and v and b and 
    cmin fields.
    Called by: mobcal_mobil2
    Calls:     mobcal_gsang, sqrt
    Uses inp,im2,ip,ofp,pgst, eo,mu,ro, v, b fields of state
  */
  double *pgst;
  double *cosx;
  double *b2max;
  double eo;
  double mu;
  double ro;
  double v;
  double b;
  double gst2;
  double recip_half_mu;
  double b2max_ig;
  
  double dbst2;
  double dbst22;
  double cmin;

  double bst2;
  double erat;
  double ang;
  double d1;
  double ten_p10;

  double boa;

  int inp;
  int im2;

  int istep;
  int ibst;

  int ibst_start;
  int ibstmax;
  
  int ibst_save;
  int success;

  int ip;
  int ig;

  int thread_id;
  int padi;

  FILE *ofp;
  FILE *lfp;

  /*
#define DBG 1
  */
  success = 1;
  pgst    = state->pgst;
  cosx    = state->cosx;
  b2max   = state->b2max;
  eo      = state->eo;
  mu      = state->mu;
  ro      = state->ro;
  inp     = state->inp;
  im2     = state->im2;
  ip      = state->ip;
  ibstmax = state->ibstmax;
  ofp     = state->ofp;
  lfp     = state->lfp;
  thread_id = state->thread_id;
  lfp      =  state->lfp;
  b       = 0.0;
  v       = 0.0;
  ten_p10 = 1.0e10;
  recip_half_mu = 2.0/mu;
  ibst   = ((int)(rmaxx/ro)) - 6;
  if (ibst < 0) ibst = 0;
  dbst2  = 1.0;
  dbst22 = dbst2 * 0.1;
  cmin   = 0.0005;
  state->cmin = cmin;
  if (im2 == 0) {
    if (thread_id == 0) {
      if (ofp) {
	fprintf(ofp,
		" set up b2 integration - integration over impact parameter\n");
	fprintf(ofp," Minium value of (1-cosx) = %le\n",cmin);
	fflush(ofp);
      }
    }
  }
  /*
    NB we need to generate cpu_time.
  */
  /*
  cpu_time(start);
  */
  for (ig=inp-1;((ig >= 0) && success);ig = ig-1) {
    /*
    if (thread_id == 0) {
      if (ofp) {
	fprintf(ofp,"ig = %d\n",ig+1);
      }
    }
    */
    gst2 = pgst[ig];
    gst2 = gst2 * gst2;
    v    = sqrt((gst2 * eo) * recip_half_mu);
    if (ig < inp-1) {
      ibst = ((int)b2max[ig+1]/dbst2) - 6;
      if (ibst < 0) {
	ibst = 0;
      }
    }
    if (ip == 1) {
      if (thread_id == 0) {
	if (ofp) {
	  fprintf(ofp," gst2 = %le v = %le\n",gst2,v);
	  fprintf(ofp,"      b          bst2       X ang       cos(X)      e ratio\n");
	}
      }
    }
    
    ibst_start = ibst;
    ibst_save  = -1;
#ifdef DBG
    if (lfp) {
      fprintf(lfp,"mobcal_b2max: before gsang loop ig = %d.\n",ig);
      fflush(lfp);
    }
#endif
    for (ibst = ibst_start; ((ibst <= ibstmax)&&(ibst_save < 0)); ibst += 1) {
      bst2 = dbst2 * ((double)ibst);
      b    = ro * sqrt(bst2);
      mobcal_gsang(state,v,b,&erat,&ang,&d1,&istep);
      cosx[ibst] = 1.0 - cos(ang);
      if (ip == 1) {
	if (thread_id == 0) {
	  if (ofp) {
	    fprintf(ofp," %le %le %le %le %le \n",
		    b, bst2, ang,cosx[ibst],erat);
	    fflush(ofp);
	  }
	}
      }
      if (ibst >= ibst_start+4) {
	if ( (cosx[ibst] < cmin) &&
	     (cosx[ibst-1] < cmin) &&
	     (cosx[ibst-2] < cmin) &&
	     (cosx[ibst-3] < cmin) &&
	     (cosx[ibst-4] < cmin)) {
	  ibst_save = ibst;
	}
      }
    }
#ifdef DBG
    if (lfp) {
      fprintf(lfp,"mobcal_b2max: after gsang loop, ig = %d\n.",ig);
      fflush(lfp);
    }
#endif
    if (ibst_save < 0) {
      success = 0;
      if (thread_id == 0) {
	if (ofp) {
	  fprintf(ofp,"ibst greater than %d\n",ibstmax);
	}
      }
    } else {
      ibst = ibst_save;
    }
    if (success) {
      b2max_ig  = ((double)(ibst - 5)) * dbst2; 
      b2max_ig  += dbst22;
      b = ro*sqrt(b2max_ig);
      mobcal_gsang(state,v,b,&erat, &ang, &d1, &istep);
      while (1.0 - cos(ang) > cmin) {
	b2max_ig = b2max_ig + dbst22;
	b = ro*sqrt(b2max_ig);
	mobcal_gsang(state,v,b,&erat, &ang, &d1, &istep);
      }
      b2max[ig] = b2max_ig;
    }
  } /* end for (ig...) */
#ifdef DBG
    if (lfp) {
      fprintf(lfp,"mobcal_b2max: after ig loop\n");
      fflush(lfp);
    }
#endif
  if (success) {
    state->b = b;
    state->v = v;
  /*
  }
  cpu_time(finish);
  if (ofp) {
    double elapsed;
    elapsed = finish - start;
    fprintf(ofp,"Time (seconds) = %ld\n",elapsed);
  }
  if (success) {
  */
    if (im2 == 0) {
      if (thread_id == 0) {
	if (ofp) {
	  fprintf(ofp,"     gst           b2max/ro2         b/A\n");
	  for (ig=0;ig<inp;ig++) {
	    b2max_ig = b2max[ig];
	    boa    = ro * sqrt(b2max_ig) * ten_p10;
	    fprintf(ofp,"%le %le %le\n",pgst[ig],b2max_ig,boa);
	  }
	  fflush(ofp);
	}
      }
    }
  }
#ifdef DBG
    if (lfp) {
      fprintf(lfp,"mobcal_b2max: before return\n");
      fflush(lfp);
    }
#endif
  return(success);
}
