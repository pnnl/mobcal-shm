#include "system_includes.h"
#include "mobcal_state_struct.h"
#include "mobcal_omega_dist.h"
void mobcal_omega_dist(struct mobcal_state_struct *state, double *omega_dist) {
  /*
    This routine determines the means of the omega(1,1), omega(1,2), omega(1,3)
    and omega(2,2) vectos and the standard deviation and standard error
    of the omega(1,1) vector. These six numbers are returne in 
    omdega_dist[0:5];
    Called by: mobcal_mobil2
    Calls:     sqrt
  */
  double mom11st;
  double mom12st;
  double mom13st;
  double mom22st;
  double sdom11st;
  double sterr;
  double *om11st;
  double *om12st;
  double *om13st;
  double *om22st;
  double recip_itn;
  double hold;
  int ic;
  int itn;
  int im2;
  int padi;
  FILE *ofp;
  FILE *efp;
  itn    = state->itn;
  im2    = state->im2;
  ofp    = state->ofp;
  om11st = state->om11st;
  om12st = state->om12st;
  om13st = state->om13st;
  om22st = state->om22st;
  
  mom11st = 0.0;
  mom12st = 0.0;
  mom13st = 0.0;
  mom22st = 0.0;
  for (ic=0;ic<itn;ic++) {
    mom11st += om11st[ic];
    mom12st += om12st[ic];
    mom13st += om13st[ic];
    mom22st += om22st[ic];
  }      
  recip_itn = 1.0/((double)itn);
  mom11st = mom11st * recip_itn;
  mom12st = mom12st * recip_itn;
  mom13st = mom13st * recip_itn;
  mom22st = mom22st * recip_itn;
  sdom11st = 0.0;
  for(ic=0;ic<itn;ic++) {
    hold = mom11st - om11st[ic];
    sdom11st += hold * hold;
  }
  sdom11st = sqrt(sdom11st * recip_itn);
  sterr    = sdom11st*sqrt(recip_itn);
  if (im2 == 0) {
    fprintf(ofp,"\nmean OMEGA*(1,1) = %le\n",mom11st);
    fprintf(ofp,"standard deviation = %le\n",sdom11st);
    fprintf(ofp,"standeard error of mean = %le\n\n",sterr);
  }
  omega_dist[0] = mom11st;
  omega_dist[1] = mom12st;
  omega_dist[2] = mom13st;
  omega_dist[3] = mom22st;
  omega_dist[4] = sdom11st;
  omega_dist[5] = sterr;
}
  
  
