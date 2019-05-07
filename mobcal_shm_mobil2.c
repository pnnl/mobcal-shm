#include "system_includes.h"
#include "mobcal_state_struct.h"
#include "mobcal_x_orient.h"
#include "mobcal_rmax_emax_r00.h"
#include "mobcal_init_gst_intgrl.h"
#include "mobcal_b2max.h"
/*
#include "mobcal_mobil2_main.h"
*/
#include "mobcal_mobil2_server.h"
#include "mobcal_mobil2_client_dbg.h"
#include "mobcal_mobil2_client.h"
#include "mobcal_mobil2_run_avg.h"
#include "mobcal_omega_dist.h"
#include "mobcal_mobil2_result.h"
#include "mobcal_shm_mobil2.h"
int mobcal_shm_mobil2(struct mobcal_state_struct *state,
		      double *mob_p, double *cs_p, double *sdevpc_p) {
  /*
    mobil2 calculates the mobility using a Lennard-Jones plus
    ion-induced dipole potential. MOBIL2 uses Monte Carlo integrations
    over orientation and impact parameter and a numerical integration
    over g* (the reduced velocity). MOBIL2 includes second order
    corrections, though these don't appear to be important.

    Called by: mobcal_shm
    Calls:     mobcal_x_orient, 
               mobcal_rmax_emax_r00,
               mobcal_init_gst_intgrl,
	       mobcal_b2max
	       mobcal_mobil2_main or mobcal_mobil2_main_dbg,
	       mobcal_mobil2_run_avg,
	       mobcal_omega_dist,
	       mobcal_mobil2_result,
	       fprintf, fflush
	       
    Uses 
    Sets
  */
  double sw1;
  double sw2;
  double dtsf1;
  double dtsf2;
  double rmax;
  double rmaxx;
  double omega_dist[6];

  int success;
  int iu2;

  int im2;
  int inwr;

  int ifail;
  int it;

  int ip;
  int igs;

  int thread_id;
  int padi;

  FILE *ofp;
  FILE *lfp;

  /*
#define DBG 1
  */
  success = 1;
  im2  	= state->im2;
  ip    = state->ip;
  inwr 	= state->inwr;
  ifail = state->ifail;
  sw1   = state->sw1;
  sw2   = state->sw2;
  dtsf1 = state->dtsf1;
  dtsf2 = state->dtsf2;
  igs   = state->igs;
  ofp   = state->ofp;
  lfp   = state->lfp;
  thread_id = state->thread_id;
#ifdef DBG
  if (lfp) {
    fprintf(lfp,"mobcal_shm_mobil2: top \n");
    fflush(lfp);
  }
#endif
  
  if (thread_id == 0) {
    if (im2 == 0) {
      if (ofp) {
	fprintf(ofp,"\nmobility calcuation by MOBIL2 (trajectory method)\n");
	fprintf(ofp,"global tractiory parameters\n");
	fprintf(ofp," sw1 = %le       sw2 = %le\n",sw1,sw2);
	fprintf(ofp," dtsf1 = %le     dtsf2 = %le\n",dtsf1,dtsf2);
	fprintf(ofp," inwr  = %d              ifail = %d\n",inwr,ifail);
	fflush(ofp);
      }
    }
  }
  it  = 0;
  iu2 = 0;
  state->iu2 = iu2;
  /*
    Determine maximum extent and orient molecule along x axis.
    This uses ox,oy and oz fields and modifies fx, fy and fz fields in state.
  */
#ifdef DBG
  if (lfp) {
    fprintf(lfp,"mobcal_shm_mobil2: before mobcal_x_orient\n");
    fflush(lfp);
  }
#endif
  success = mobcal_x_orient(state, &rmax);

  if (success) {
    if (thread_id == 0) {
      if (ip == 1) {
	fprintf(ofp,"\n");
	fflush(ofp);
      }
    }
#ifdef DBG
  if (lfp) {
    fprintf(lfp,"mobcal_shm_mobil2: before mobcal_rmax_emax_r00\n");
    fflush(lfp);
  }
#endif
    /*
      Determine rmax, emax and r00 along x, y, and z directions.
    */
    mobcal_rmax_emax_r00(state, rmax, &rmaxx);
    /*
      Set up integration over gst filling wgst and pgst vectors.
    */
#ifdef DBG
  if (lfp) {
    fprintf(lfp,"mobcal_shm_mobil2: before mobcal_init_gst_intgrl\n");
    fflush(lfp);
  }
#endif
    mobcal_init_gst_intgrl(state);
    /*
      Determine b2max.
    */
#ifdef DBG
  if (lfp) {
    fprintf(lfp,"mobcal_shm_mobil2: before mobcal_b2max\n");
    fflush(lfp);
  }
#endif
    success = mobcal_b2max(state,rmaxx);
  }
#ifdef DBG
  if (lfp) {
    fprintf(lfp,"mobcal_shm_mobil2: after mobcal_b2max\n");
    fflush(lfp);
  }
#endif
  if (success) {
    /*
      Main loop where we determine the q1st, q2st,
      om11st, om12st, om13st and om22st vectors.
      if printing is turned off we can be more efficient.
    */
    if (state->thread_id == 0) {
#ifdef DBG
  if (lfp) {
    fprintf(lfp,"mobcal_shm_mobil2: before mobcal_mobil2_server\n");
    fflush(lfp);
  }
#endif
      mobcal_mobil2_server(state);
#ifdef DBG
  if (lfp) {
    fprintf(lfp,"mobcal_shm_mobil2: after mobcal_mobil2_server\n");
    fflush(lfp);
  }
#endif
      if (im2 == 0) {
	/*
	  Calculate running averages.
	*/
	mobcal_mobil2_run_avg(state);
      }
#ifdef DBG
  if (lfp) {
    fprintf(lfp,"mobcal_shm_mobil2: after mobcal_mobil2_run_avg\n");
    fflush(lfp);
  }
#endif
      /*
	Calulate the means of the omega(1,1), omega(1,2), omega(1,3),
	and omega(2,2) vectors and the standard deviation and standard
	error for the omega(1,1) vector.
      */
      mobcal_omega_dist(state,omega_dist);
#ifdef DBG
  if (lfp) {
    fprintf(lfp,"mobcal_shm_mobil2: after mobcal_omega_dist\n");
    fflush(lfp);
  }
#endif
      /* 
	 Use omegas to obtain higher order correction factor to mobility.
      */
      mobcal_mobil2_result(state,omega_dist,mob_p,cs_p,sdevpc_p);
#ifdef DBG
  if (lfp) {
    fprintf(lfp,"mobcal_shm_mobil2: after mobcal_mobil2_result\n");
    fflush(lfp);
  }
#endif
    } else {
      if ((ip == 1) || (igs == 1)) {
	mobcal_mobil2_client_dbg(state);
      } else {
#ifdef DBG
  if (lfp) {
    fprintf(lfp,"mobcal_shm_mobil2: before mobcal_mobil2_client\n");
    fflush(lfp);
  }
#endif
	mobcal_mobil2_client(state);
#ifdef DBG
  if (lfp) {
    fprintf(lfp,"mobcal_shm_mobil2: after mobcal_mobil2_client\n");
    fflush(lfp);
  }
#endif
      }
    }
  }
  return(success);
}
