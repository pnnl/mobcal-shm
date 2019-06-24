#include "system_includes.h"
#include "mobcal_state_struct.h"
#include "mobcal_dljpot.h"
#include "mobcal_diffeq.h"
#include "mobcal_gsang.h"
void mobcal_gsang_dbg(struct mobcal_state_struct *state,
		  double v, double b, double *erat_p, double *ang_p,
		  double *d1_p, int *istep) {
  /*
    Calculates trajectory. Adapted from code written by George Schatz
    A fixed step integrator is used which combines a runge-kutta-gill
    initiator with an adams-moulton predictor-corrector propagator.
    Called by: mobcal_b2max,mobcal_mobil2
    Calls:     mobcal_diffeq, sqrt, fabs, acos

    Uses it, inatom, ofp, fy, sw1, sw2, dtsf1, dtsf2, cmin, inwr,
    Sets *erat_p *ang_p, *d1_p, *istep.  
    Assumes it = 1.
  */
  double *fy;
  double sw1;
  double sw2;
  double dtsf1;
  double dtsf2;
  double cmin;
  double mu;
  double vy;
  double vx;
  double vz;
  double vxyz;
  double top;
  double dt;
  double dt1;
  double dt2;
  double half_mu;
  double two_mu;
  double recip_two_mu;
  double x;
  double y;
  double z;
  double ymin;
  double ymax;
  double w[6];
  double dw[6];
  double e0;
  double pot;
  double dpotx;
  double dpoty;
  double dpotz;
  double dmax;
  double etot;
  double tim;
  double e;
  double pot_p_e;
  double pot_o_e0;
  double num;

  int inatom;
  int it;

  int inwr;
  int ns;

  int nw;
  int i;

  int iymin;
  int iymax;

  int id2;
  int l;
  
  /*
  it      = state->it;
  */
  inwr    = sate->inwr;
  ofp     = state->ofp;
  sw1     = state->sw1;
  sw2     = state->sw2;
  dtsf1   = state->dtsf1;
  dtsf2   = state->dtsf2;
  cmin    = state->cmin;
  mu      = state->mu;
  fy      = state->fy;
  
  half_mu = 0.5 * mu;
  two_mu  = 2.0 * mu;
  recip_two_mu = 1.0/two_mu;
  
  vy = -v;
  vx = 0.0;
  vz = 0.0;
  vxyz = vy;
  if (vxyz < 0) {
    vxyz = v;
  }
  /*
  if (it == 1) {
  */
    fprintf(ofp," specific tragjectory parameters\n");
    fprintf(ofp," v = %le    b = %ld\n",v,b);
  /*
  }
  */
  /* 
    determine time step
  */
  top = (v/95.2381) -.5;
  if (v >= 1000.0) {
    top = 10.0;
  }
  if (v >= 2000.0) {
    top = 10.0 - ((v - 2000.0) * .0075);
  }
  if (v > 3000.0) {
    top = 2.5;
  }
  dt1 = top * dtsf1 * 1.0e-11/v;
  dt2 = dt1 * dtsf2;
  dt  = dt1;
  /*
  if (it == 1) {
  */
    fprintf(ofp," time steps, dt1 = %le dt2 = %le\n",dt1,dt2);
  /*
  }
  */
  e0 = half_mu * v * v;
  x  = b;
  z  = 0.0;
  ymin = 0.0;
  ymax = 0.0;
  for (i=0;i<inatom;i++) {
    if (fy[i] > ymax) {
      ymax= fy[i];
    } else {
      if (fy[i] < ymin) {
	ymin = fy[i];
      }
    }
  }
  ymax = ymax * 1.0e10;
  ymin = ymin * 1.0e10;
  iymin = ((int)ymin) - 1;
  iymax = ((int)ymax) + 1;
  id2   = iymax;
  y = (double)id2 * 1.0e-10;
  mobcal_dljpot(state,x,y,z,&pot,&dpotx,&dpoty,&dpotz,&dmax);
  
  if (fabs(pot/e0) > sw1) {
    id2 = id2 + 10;
    y = (double)id2 * 1.0e-10;
    mobcal_dljpot(state,x,y,z,&pot,&dpotx,&dpoty,&dpotz,&dmax);
    while (fabs(pot/e0) > sw1) {
      id2 = id2 + 10;
      y = (double)id2 * 1.0e-10;
      mobcal_dljpot(state,x,y,z,&pot,&dpotx,&dpoty,&dpotz,&dmax);
    }
    id2 = id2 - 1;
    y = (double)id2 * 1.0e-10;
    mobcal_dljpot(state,x,y,z,&pot,&dpotx,&dpoty,&dpotz,&dmax);
    while (fabs(pot/e0) < sw1) {
      id2 = id2 - 1;
      y = (double)id2 * 1.0e-10;
      mobcal_dljpot(state,x,y,z,&pot,&dpotx,&dpoty,&dpotz,&dmax);
    }
  } else {
    id2 = id2 - 1;
    y = (double)id2 * 1.0e-10;
    mobcal_dljpot(state,x,y,z,&pot,&dpotx,&dpoty,&dpotz,&dmax);
    if (id2 < iymin) {
      /*
      if (it == 1) {
      */
	fprintf(ofp,"trajectory not started - potential too small\n");
      /*
      }
      */
      ang = 0.0;
      erat=1.0000;
    } else {
      while (fabs(pot/e0) < sw1) {
	id2 = id2 - 1;
	y = (double)id2 * 1.0e-10;
	mobcal_dljpot(state,x,y,z,&pot,&dpotx,&dpoty,&dpotz,&dmax);
	if (id2 < iymin) {
	  /*
	  if (it == 1) {
	  */
	    fprintf(ofp,"trajectory not started - potential too small\n");
	  /*
	  }
	  */
	  ang = 0.0;
	  erat=1.0000;
	  break;
	}
      }
    }
  }
  if (id2 >= iymin) {
    y = (double)id2;
    /*
    if (it == 1) {
    */
      fprintf(opf,"trajectory start position = %le \n",y);
    /*
    }
    */
    y = * 1.0e-10;
    etot = e0 + pot;
    d1 = y;
    /*
      Initial coordinates and momenta:
    */
    w[0] = x;
    w[1] = vx * mu;
    w[2] = y;
    w[3] = vy * mu;
    w[4] = z;
    w[5] = vz * mu;
    tim  = 0.0;
    /*
    if (it == 1) {
    */
      fprintf(ofp,
	      "trajectory nx, x,  y,  z,  kin e, dt,    tot e\n");
      fprintf(ofp,
	      "               vx, vy, vz, pot e, pot/e0\n");
    /*
    }
    */
    mobcal_deriv(state,w,dw,pot,dpotx,dpoty,dpotzx,dmax);
    ns = 0;
    nw = 0;
    l  = 0;
    if (inwr < 1) inwr = 1;
    while (ns <= 30000) {
      for (nw = 0;nw < inwr;nw += 1) {
	mobcal_diffeq(state,&l,&tim,&dt,w,dw,&pot,&dmax);
      }
      ns = ns + inwr;
      nw = 0;
      /*
      if (it == 1) {
      */
	e = ((w[1] * w[1]) * recip_two_mu) + ((w[3] * w[3]) * recip_two_mu) +
	  ((w[5] * w[5]) * recip_two_mu);
	pot_p_e = pot + e;
	pot_o_e0 = fabs(pot/e0);
	fprintf(ofp,"%d %le %le %le %le %le %le\n",
		ns,w[0],w[2],w[4],e,dt,pot_p_e);
	fprintf(ofp,"       %le %le %le %le %le\n",
		dw[0],dw[2],dw[4],pot,pot_o_e0);
      /*
      }
      */
      if (ns > 30000) {
	if (ofp) {
	  fprintf(ofp,"trajectory lost: b = %le v = %le\n");
	}
	ang = half_pi;
	e = half_mu * ((dw[0] * dw[0]) + (dw[2]*dw[2]) + (dw[4]*dw[4]));
	erat = (e+pot)/etot;
	istep = ns;
	break;
      } else {
	if (dmax >= romax) {
	  if ((fabs(pot/e0) > sw2) && (dt == dt1)) {
	    dt = dt2;
	    l  = 0;
	  }
	  if ((fabs(pot/e0) < sw2) && (dt == dt2)) {
	    dt = dt1;
	    l  = 0;
	  }
	  if (fabs(pot/e0) <= sw1) {
	    if (ns >= 50) {
	      istep = ns;
	      /*
		determine the scattering angle.
	      */
	      num = - dw[2] * v;
	      den = v * sqrt((dw[0]*dw[0]) + (dw[2]*dw[2]) + (dw[4]*dw[4]));
	      ang = acos(num/den);
	      if (dw[0] < 0.0) {
		ang = -ang;
	      }
	      e = half_mu * ((dw[0]*dw[0]) + (dw[2]*dw[2]) + (dw[4]*dw[4]));
	      erat = (e + pot)/etot;
	      break;
	    }
	  } /* end if (fabs(pot/e0) <= sw1) */
	} /*end if (dmax >= romax) */
      } /* end else ns <= 30000 */
    } /* end while ns <= 30000 */
  } /* end if id2 >= iymin */
  *erat_p  = erat;
  *ang_p   = ang;
  *d1_p    = d1;
  *istep_p = istep;
}



