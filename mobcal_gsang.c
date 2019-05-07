#include "system_includes.h"
#include "mobcal_state_struct.h"
#include "mobcal_dljpot_only.h"
#include "mobcal_diffeq.h"
#include "mobcal_deriv.h"
#include "mobcal_gsang.h"
void mobcal_gsang(struct mobcal_state_struct *state,
		  double v, double b, double *erat_p, double *ang_p,
		  double *d1_p, int *istep_p) {
  /*
    Calculates trajectory. Adapted from code written by George Schatz
    A fixed step integrator is used which combines a runge-kutta-gill
    initiator with an adams-moulton predictor-corrector propagator.
    Called by: mobcal_b2max,mobcal_mobil2
    Calls:     mobcal_diffeq, mobcal_deriv, mobcal_dljpot_only,
               sqrt, fabs, acos

    Uses it, inatom, ofp, fy, sw1, sw2, dtsf1, dtsf2, cmin, inwr,
    Sets *erat_p *ang_p, *d1_p, *istep.  
  */
  double *fy;
  double xfy;
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
  double romax;
  double x;
  double y;
  double z;
  double ymin;
  double ymax;
  double w[6];
  double dw[6];
  double e0;
  double pot;
  double dmax;
  double etot;
  double tim;
  double e;
  double num;
  double ang;
  double erat;
  double d1;
  double half_pi;
  double den;

  int inatom;
  int it;

  int inwr;
  int ns;

  int nw;
  int iatom;

  int iymin;
  int iymax;

  int id2;
  int l;

  int istep;
  int dljpot_input_stream_spacing;
  
  int i;
  int padi;

  FILE *ofp;
  FILE *lfp;

  it      = state->it;
  inwr    = state->inwr;
  ofp     = state->ofp;
  sw1     = state->sw1;
  sw2     = state->sw2;
  dtsf1   = state->dtsf1;
  dtsf2   = state->dtsf2;
  cmin    = state->cmin;
  mu      = state->mu;
  romax   = state->romax;
  inatom  = state->inatom;
  /*
  fy      = state->fy;
  */
  lfp     = state->lfp;
  half_pi   = state->half_pi;
  half_mu = 0.5 * mu;
  two_mu  = 2.0 * mu;
  recip_two_mu = 1.0/two_mu;
  istep   = 0;
  vy = -v;
  vx = 0.0;
  vz = 0.0;
  vxyz = vy;
  if (vxyz < 0) {
    vxyz = v;
  }
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
  e0 = half_mu * v * v;
  x  = b;
  z  = 0.0;
  ymin = 0.0;
  ymax = 0.0;
  dljpot_input_stream_spacing = 11 * _IVEC_LEN_;
  fy = state->dljpot_input_stream + _IVEC_LEN_;
  for (iatom=0;iatom<inatom;iatom+= _IVEC_LEN_) {
    for (i= 0;i<_IVEC_LEN_;i++) {
      xfy = fy[i];
      if (xfy > ymax) {
	ymax= xfy;
      } else {
	if (xfy < ymin) {
	  ymin = xfy;
	}
      }
    }
    fy += dljpot_input_stream_spacing;
  }
  ymax = ymax * 1.0e10;
  ymin = ymin * 1.0e10;
  iymin = ((int)ymin) - 1;
  iymax = ((int)ymax) + 1;
  id2   = iymax;
  y = (double)id2 * 1.0e-10;
  mobcal_dljpot_only(state,x,y,z,&pot);
  
  if (fabs(pot/e0) > sw1) {
    /* 302 */
    id2 = id2 + 10;
    y = (double)id2 * 1.0e-10;
    mobcal_dljpot_only(state,x,y,z,&pot);
    while (fabs(pot/e0) > sw1) {
      id2 = id2 + 10;
      y = (double)id2 * 1.0e-10;
      mobcal_dljpot_only(state,x,y,z,&pot);
    }
    /* 301 */
    id2 = id2 - 1;
    y = (double)id2 * 1.0e-10;
    mobcal_dljpot_only(state,x,y,z,&pot);
    while (fabs(pot/e0) < sw1) {
      id2 = id2 - 1;
      y = (double)id2 * 1.0e-10;
      mobcal_dljpot_only(state,x,y,z,&pot);
    }
  } else {
    /* 300 */
    id2 = id2 - 1;
    y = (double)id2 * 1.0e-10;
    mobcal_dljpot_only(state,x,y,z,&pot);
    if (id2 < iymin) {
      ang = 0.0;
      erat=1.0000;
    } else {
      while (fabs(pot/e0) < sw1) {
	id2 = id2 - 1;
	y = (double)id2 * 1.0e-10;
	mobcal_dljpot_only(state,x,y,z,&pot);
	if (id2 < iymin) {
	  ang = 0.0;
	  erat=1.0000;
	  break;
	}
      }
    }
  }
  if (id2 >= iymin) {
    /* 304*/
    y = (double)id2;
    y = y * 1.0e-10;
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
    mobcal_deriv(state,w,dw,&pot,&dmax);
    ns = 0;
    nw = 0;
    l  = 0;
    istep = 0;
    if (inwr < 1) inwr = 1;
    while (ns <= 30000) {
      for (nw = 0;nw < inwr;nw += 1) {
	mobcal_diffeq(state,&l,&tim,&dt,w,dw,&pot,&dmax);
      }
      ns = ns + inwr;
      istep = ns;
      nw = 0;
      if (ns > 30000) {
	if (lfp) {
	  fprintf(lfp,"trajectory lost: b = %le v = %le\n",b,v);
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
	      ang =  acos(num/den);
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
    *d1_p    = d1;
  } /* end if id2 >= iymin */
  *erat_p  = erat;
  *ang_p   = ang;
  *istep_p = istep;
}



