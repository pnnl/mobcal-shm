#include "system_includes.h"
#include "mobcal_state_struct.h"
#include "mobcal_dljpot_only.h"
#include "mobcal_rmax_emax_r00.h"
void mobcal_rmax_emax_r00(struct mobcal_state_struct *state, double rmax,
			  double *rmaxx_p) {
  /*
    determine rmax, emax and r00 along x, y, and z directions.
    only rmaxx is actually used for anything other than
    printing.
    Called by: mobcal_mobil2
    Calls      mobcal_dljpot_only
  */
  double emaxx; 
  double rmaxx;
  double r00x;
  double emaxy;
  double rmaxy;
  double r00y;
  double emaxz;
  double rmaxz;
  double r00z;
	
  double xe;
  double recip_xe;
  double ddd;
  double romax;
  
  double pot;
  double x;
  double y;
  double z;
  double xemax;
  double xrmax;
  double xr00;
  double yemax;
  double yrmax;  
  double yr00;
  double zemax;
  double zrmax;  
  double zr00;
  double ten_p10;

  int irn;
  int ir;
  int im2;
  int thread_id;

  FILE *ofp;
  FILE *efp;

  irn      = state->irn;
  xe       = state->xe;
  im2      = state->im2;
  romax    = state->romax;
  ofp      = state->ofp;
  thread_id = state->thread_id;
  recip_xe = 1.0/xe;
  ten_p10  = 1.0e10;
  ddd      = (rmax + romax)/((double)irn);
  y     = 0.0;
  z     = 0.0;
  pot   = 0.0;
  emaxx = 0.0;
  rmaxx = 0.0;
  rmaxy = 0.0;
  rmaxz = 0.0;
  r00x  = 0.0;
  r00y  = 0.0;
  r00z  = 0.0;

  for (ir=1;((ir<=irn) & (pot <= 0.0));ir++) {
    /*
      x = rmax + romax -(dfloat(ir) * ddd))
      why not (irn - ir) * (rmax + romax)/1000
    */
    x = ((double)(irn-ir)) * ddd;
    mobcal_dljpot_only(state,x,y,z,&pot);
    if (pot <= 0.0) {
      r00x  = x;
      if (pot < emaxx) {
	rmaxx = x;
	emaxx = pot;
      }
    }
  } /* end for (ir ...) */
  if (im2 == 0) {
    xemax = emaxx * recip_xe;
    xrmax = rmaxx * ten_p10;
    xr00  = r00x  * ten_p10;
    if (thread_id == 0) {
      fprintf(ofp,"along x axis emax = %le eV, rmax = %leA, r00  = %leA\n",
	      xemax, xrmax, xr00);
    }
    /* NB the rest of this routine only computes values that are
       printetd out if im2 is 0, hence it is included inside this
       if block.
    */
    x     = 0.0;
    z     = 0.0;
    pot   = 0.0;
    emaxy = 0.0;
    for (ir=1;((ir<=irn) & (pot <= 0.0));ir++) {
      y = ((double)(irn-ir)) * ddd;
      mobcal_dljpot_only(state,x,y,z,&pot);
      if (pot <= 0.0) {
	r00y  = y;
	if (pot < emaxy) {
	  rmaxy = y;
	  emaxy = pot;
	}
      }
    } /* end for (ir ...) */
    yemax = emaxy * recip_xe;
    yrmax = rmaxy * ten_p10;
    yr00  = r00y  * ten_p10;
    if (thread_id == 0) {
      fprintf(ofp,"along y axis emax = %le eV, rmax = %leA, r00 = %leA\n",
	      yemax, yrmax, yr00);
    }
    x     = 0.0;
    y     = 0.0;
    pot   = 0.0;
    emaxz = 0.0;
    for (ir=1;((ir<=irn) & (pot <= 0.0));ir++) {
      z = ((double)(irn-ir)) * ddd;
      mobcal_dljpot_only(state,x,y,z,&pot);
      if (pot <= 0.0) {
	r00z  = z;
	if (pot < emaxz) {
	  rmaxz = z;
	  emaxz = pot;
	}
      }
    } /* end for (ir ...) */
    zemax = emaxz * recip_xe;
    zrmax = rmaxz * ten_p10;
    zr00  = r00z  * ten_p10;
    if (thread_id == 0) {
      fprintf(ofp,"along z axis emax = %le eV, rmax = %leA, r00 = %leA\n",
	      zemax, zrmax, zr00);
    }
  } /* end if (im2 == 0) */
  *rmaxx_p = rmaxx;
}
