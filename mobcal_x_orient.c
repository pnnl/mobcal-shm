#include "system_includes.h"
#include "mobcal_state_struct.h"
#include "mobcal_rotate_phi_gamma.h"
#include "mobcal_x_orient.h"
int mobcal_x_orient(struct mobcal_state_struct *state, double *rmax_p) {
  /*
    Orient the molecule along the x axis.
    Called by: mobcal_mobil2
    Calls:     mobcal_rotate_phi_gamma, sqrt, acos, fabs
    Uses ox, oy, yz, im2, inatom, half_pi, two_pi fields of state.
    Sets fx, fy, and fz fields of state and *rmax_p.
    Returns 1 on success, 0 on failure to rotate.
  */
  double *ox;
  double *oy;
  double *oz;
  double *fx;
  double *fy;
  double *fz;
  
  double half_pi;
  double two_pi;

  double rmax;
  double rzy;
  
  double gamma;
  double phi;

  double rmax2;
  double r2;
  double xx;
  double xy;
  double xz;
  double rzy2;
  double rzy2_hold;
  double ox_hold;
  double oy_hold;
  double oz_hold;
  double ten_m10;
  double ten_m20;
  double p;
  double r;
  double rxy;
  double hold;
  double theta;

  int    im2;
  int    ihold;

  int    iatom;
  int    inatom;

  int    success;
  int    ip;

  int    thread_id;
  int    dljpot_input_stream_spacing;

  int    i;
  int    ivec_len;

  int    interval;
  int    hold_pos;

  FILE   *ofp;
  FILE   *efp;

  success = 1;
  dljpot_input_stream_spacing = 11 * _IVEC_LEN_;
  ivec_len                    = _IVEC_LEN_;
  im2     = state->im2;
  inatom  = state->inatom;
  half_pi = state->half_pi;
  two_pi  = state->two_pi;
  theta     = 0.0;
  ox        = state->ox;
  oy        = state->oy;
  oz        = state->oz;
  /*
  fx        = state->fx;
  fy        = state->fy;
  fz        = state->fz;
  */
  fx          = state->dljpot_input_stream;
  fy          = fx + _IVEC_LEN_;
  fz          = fy + _IVEC_LEN_;
  ip          = state->ip;
  ofp         = state->ofp;
  thread_id   = state->thread_id;
  ihold = 0;
  rmax2 = 0.0;
  rzy2_hold = 0.0;
  ox_hold   = 0.0;
  oy_hold   = 0.0;
  oz_hold   = 0.0;
  ten_m20   = 1.0e-20;
  ten_m10   = 1.0e-10;
  hold_pos  = 0;
  interval  = 0;
  for (iatom=0;iatom<inatom;iatom++) {
    xx = ox[iatom];
    xy = oy[iatom];
    xz = oz[iatom];
    rzy2 = (xy*xy) + (xz*xz);
    r2 =  ((xx * xx) + rzy2);
    if (r2 > rmax2) {
      rmax2 = r2;
      rzy2_hold = rzy2;
      oz_hold   = xz;
      oy_hold   = xy;
      ox_hold   = xx;
      ihold = hold_pos;
    }
    hold_pos += 1;
    interval += 1;
    if(interval == _IVEC_LEN_) {
      interval = 0;
      hold_pos += (dljpot_input_stream_spacing - _IVEC_LEN_);
    }
  }
  rmax = sqrt(rmax2);
  rzy = sqrt(rzy2_hold);
  if (rzy == 0) {
    success = 0;
    fprintf(stderr,"mobcal_x_orient: Error rzy = 0\n");
    fflush(stderr);
  } else {
    phi = acos(oz_hold/rzy);
    phi = phi + half_pi;
    /*
      NB tint the original code the if (oy_hold > 0) ...
      is coded as
      if (oy_hold < 0) phi = 2*pi - phi
      phi = 2*pi - phi
      
      below is slightly clearer and takes fewer opearions 
      1 on average (1 mutiply + 1add)/2 vs  3 on average:
      1 multiply plus 1 add + (1multiply + add)/2
    */
    if (oy_hold >= 0) {
      phi = two_pi - phi;
    }
  }
  if (success) {
    /*
      now rotate just the extremum point by phi to determine gamma.
    */
    p = cos(phi);
    r = sin(phi);
    xx = ox_hold;
    xy = oy_hold * p + oz_hold*r;
    xz = -oy_hold * r + oz_hold*p;
    rxy = sqrt(xx * xx + xy * xy);
    gamma = acos(xx/rxy);
    if (xy >= 0.0) {
      gamma = two_pi - gamma;
    }
    if (im2 == 0) {
      state->iu3 = 1;
    }
    if (ip == 1) {
      state->iu2 = 1;
    }
    mobcal_rotate_phi_gamma(state,phi,gamma);
    hold = fx[ihold]/rmax;
    state->iu3 = 0;
    state->iu2 = 0;
    if (thread_id == 0) {
      if (ofp) {
	fprintf(ofp,"hold  = %le, gamma = %le\n",hold,gamma);
	fprintf(ofp,"theta = %le, phi   = %le\n",theta,phi);
      }
    }
    xx = fabs(hold - 1.0);
    xy = fabs(fy[ihold]);
    xz = fabs(fz[ihold]);
    if ((xx > ten_m10) || (xy > ten_m20) || (xz > ten_m20)) {
      success = 0;
      if (thread_id == 0) {
	if (ofp) {
	  fprintf(ofp,"Problem orienting along x axis\n");
	  
	  for (iatom=0;iatom<inatom;iatom+=_IVEC_LEN_) {
	    for (i=0;i<_IVEC_LEN_;i++) {
	      xx = fx[i];
	      xy = fy[i];
	      xz = fz[i];
	      hold = sqrt((xx*xx) + (xy*xy) + (xz*xz));
	      fprintf(ofp,"%d     %le     %le      %le     %le\n",
		      iatom,xx,xy,xz,hold);
	    }
	    fx += dljpot_input_stream_spacing;
	    fy =  fx + _IVEC_LEN_;
	    fz =  fy + _IVEC_LEN_;
	  }
	  fflush(ofp);
	}
      } /* end if (thread_id == 0) */
    }
  }
  *rmax_p = rmax;
  return(success);
}


