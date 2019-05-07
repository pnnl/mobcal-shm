#include "system_includes.h"
#include "mobcal_state_struct.h"
#include "mobcal_struct_asym.h"
int mobcal_struct_asym(struct mobcal_state_struct *state, double *asym_p) {
  /* 
    determine structural asymmetry parameter.
    This is relatively compute intensive, and we have inlined 
    the rotate call.
    Called by: mobcal_ncoord
    Calls:     cos, sin, acos, atan, sqrt
               
  */
  double *work;
  double *ox;
  double *oy;
  double *oz;
  double xfx;
  double xfy;
  double xfz;
  double cang;
  double recip_cang;
  double asymp;
  double theta;
  double phi;
  double gamma;
  double xyzsq;
  double yzsq;
  double xyzsum;
  double yzsum;
  double quarter_pi;
  double two_pi;
  double p;
  double r;
  double g;
  double h;
  double at0;
  double at1;
  double at2;
  double at3;
  double at4;
  double at5;
  double xx;
  double xy;
  double xz;
  double hold;

  int    iphi;
  int    igamma;

  int    inatom;
  int    iatom;

  int    incx;
  int    success;

  char   no_trans;
  char   cdm0;

  char   cdm1;
  char   cdm2;

  success    = 1;
  incx       = 1;
  cang       = state->cang;
  recip_cang = 1.0/cang;
  asymp      = 0.0;
  theta      = 0.0;
  quarter_pi = atan(1.0);
  two_pi     = state->pi + state->pi;
  no_trans   = 'N';
  inatom     = state->inatom;
  work       = state->asym_work;
  ox         = state->ox;
  oy         = state->oy;
  oz         = state->oz;
  /*
    In the doubly nested loop below, with theta = 0
    the first rotation in mobcal_rotate is the identity rotation
    (that is u = 1, v=0 (see mobcal_ncoord documentation)
    Then our rotation matrix becomes 

     g = cos(gamma) h = sin(gamma)
     p = cos(phi)   r = sin(phi)
     
             g      -hp    -hr
	     h       gp     gr
	     0       -r     p

     
    We observe that the last row of the transformation is invariant 
    with respect to gamma,
    hence we need compute fz only once for each phi.

  */
  for (iphi = 0; iphi <= 180; iphi += 2) {
    phi   = ((double)iphi) * recip_cang;
    p     = cos(phi);
    r     = sin(phi);
    for (iatom=0;iatom<inatom;iatom++) {
      xfz       = -r*oy[iatom] + p * oz[iatom];
      work[iatom] = xfz * xfz;
    }
    /*
             g      -hp    -hr
	     h       gp    +gr
	     0       -r     p

    */

    for (igamma = 0; igamma <= 360; igamma += 2) {
      gamma = ((double)igamma) * recip_cang;
      g = cos(gamma);
      h = sin(gamma);
      at0 = g;
      at1 = -h*p;
      at2 = -h*r;
      at3 = h;
      at4 = g*p;
      at5 = g*r;

      xyzsum = 0.0;
      yzsum  = 0.0;
      for (iatom=0;iatom<inatom;iatom++) {
	xx = ox[iatom];
	xy = oy[iatom];
	xz = oz[iatom];
	xfx       = (at0 * xx) + (at1 * xy) + (at2 * xz);
	xfy       = (at3 * xx) + (at4 * xy) + (at5 * xz);
	xfx       = xfx * xfx;
	xfy       = xfy * xfy;
	yzsq      = xfy + work[iatom];
	xyzsq     = xfx + yzsq;
	xyzsum    += sqrt(xyzsq);
	yzsum     += sqrt(yzsq);
      }
      hold = quarter_pi * xyzsum/yzsum;
      if (hold > asymp) {
	asymp = hold;
      }
    }
  }
  *asym_p = asymp;
  return(success);
}
