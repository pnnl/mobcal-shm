#include "system_includes.h"
#include "mobcal_state_struct.h"
#include "dgemm_.h"
#include "mobcal_rotate_phi_gamma.h"
void mobcal_rotate_phi_gamma(struct mobcal_state_struct *state,
			    double phi,double gamma) {
  /*
    from mobcal_rotate:

   new
      x'"     gu-hpv  -gv-hpu  -hr     x
      y'"  = (hu+gpv  -hv+gpu   gr) *  y
      z'"     -rv       -ru     p      z


    let us label the combined rotation matrix a.

     Here if theta = 0 then u = 1, v = 0 and a becomes
     g = cos(gamma) h = sin(gamma)
     p = cos(phi)   r = sin(phi)
     
             g      -hp    -hr
	     h       gp     gr
	     0      -r      p

      Called by: mobcal_x_orient
      Calls:     sin,cos, dgemm_

  */       
  double *ap;
  double *fx;
  /*
  double *fy;
  double *fz;
  */
  double *ox;
  /*
  double *oy;
  double *oz;
  */
  double g;
  double h;
  double p;
  double r;
  double a[9];
  double alpha;
  double beta;

  int inatom;
  int success;

  int m;
  int n;
  int k;
  int dljpot_input_stream_spacing;

  int i;
  int ivec_len;

  char no_trans;
  char c1;
  char c2;
  char c3;
  char c4;
  char c5;
  char c6;
  char c7;

  success = 1;
  ivec_len = _IVEC_LEN_;
  inatom  = state->inatom;
  fx      = state->fx;
  ox      = state->ox;
  ap = (double*)&a[0];
  g = cos(gamma);
  h = sin(gamma);
  p = cos(phi);
  r = sin(phi);
  no_trans = 'N';
  beta = 0.0;
  alpha = 1.0;

  a[0] = g;
  a[1] = -h*p;
  a[2] = -h*r;
  
  a[3] = h;
  a[4] = g*p;
  a[5] = g*r;

  a[6] = 0.0;
  a[7] = -r;
  a[8] = p;

  m     = state->vec_len;
  n     = 3;
  k     = 3;
  /*
  dgemm_(&no_trans,&no_trans,&m,&n,&k,&alpha,ox,&m,ap,&k,&beta,fx,&m);
  */
  fx     = state->dljpot_input_stream;
  dljpot_input_stream_spacing = 11 * _IVEC_LEN_;
  for (i=0;i<m;i+=_IVEC_LEN_) {
    dgemm_(&no_trans,&no_trans,&ivec_len,&n,&k,&alpha,ox,&m,ap,&k,&beta,fx,&ivec_len);
    ox += ivec_len;
    fx += dljpot_input_stream_spacing;
  }

  /*
  for (iatom=0;iatom<inatom;iatom++) {
    xx = ox[iatom];
    xy = oy[iatom];
    xz = oz[iatom];
    fx[iatom] = (at0 * xx) + (at1 * xy) + (at2 * xz);
    fy[iatom] = (at3 * xx) + (at4 * xy) + (at5 * xz);
    fz[iatom] = (at7 * xy) + (at8 * xz);
  }
  */
}

	
