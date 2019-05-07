#include "system_includes.h"
#include "dgemm_.h"
#include "mobcal_state_struct.h"
#include "mobcal_print_rotated_coords.h"
#include "mobcal_rotate.h"
void mobcal_rotate(struct mobcal_state_struct *state,
		  double theta, double phi,double gamma) {
  /*
    Rotates the cluster/molecule.  
    Called by: mobcal_rantate,
    Calls:     sin, cos, dgemm_, mobcal_print_rotated_coords

    The math in the original rotate is absolutely incredibly bad.

    This routine rotates first in the xy plane by theta,
    then in the yz plane by phi and then finally in the xy plane again
    by gamma? while that seems odd to me consider rotating 
    a point x,y,z first in the xy plane by theta  (z is unchanged)

       x'      cos(theta)  -sin(theta)  0     x
       y'  =  (sin(theta)   cos(theta)  0 ) * y
       z'          0           0        1     z

    a matrix (3x3)  vector (3x1) product. The same rotation
    matrix is used for all points.

       Then a rotation in th yz plane by phi is given by (x' is unchanged)

      x"      1       0         0          x'
      y"  =  (0    cos(phi)   sin(phi) ) * y'
      z"      0    -sin(phi)   cos(phi)    z'

      another three by thre matrix product.,

    And then finally another rotation is done in the xy plane by gamma.

      x'"    cos(gamma)  -sin(gamma)  0    x"
      y'" = (sin(gamma)   cos(gamma)  0) * y"
      z'"        0             0      1    z"

     For convenience let g = cos(gamma), h = sin(gamma), 
                         p = cos(phi),   r = sin(phi), 
		         u = cos(theta), v = cos(theta)

     Then stringing together succesive rotations above we get

      x'"     g  -h  0     1  0  0     u -v 0    x
      y'" =  (h   g  0) * (0  p  r) * (v  u 0) * y
      z'"     0   0  1     0 -r  p     0  0 1    z
			  
      Multiplying out the matrices first the left most 2 gives

      x'"    g  -hp  -hr     u -v 0    x
      y'" = (h   gp   gr) * (v  u 0) * y
      z'"    0   -r   p      0  0 1    z

      Combining these remaining two rotation matrices we get

      x'"     gu-hpv  -gv-hpu  -hr     x
      y'"  = (hu+gpv  -hv+gpu   gr) *  y
      z'"     -rv       -ru     p      z

       we note that pv and pu are shared among to expressions for the 
       elements in the combined rotation matrix-matrix product, hence
       its elements may be computed in 14 muliplies and 4 adds (or 6 adds if we 
       count the negations as adds).

      Then to apply this rotation matrix to n points in three space takes
        
      14 + 9n multiplies and 6n + 6 adds and 3 sin and 3 cos evaluations

      Furthermore all these computations maybe composed into a dgemm
      call if we arrange the data correctly.
        
      Compare this to the

      3n cos, sin and acos evaluations, 3n square roots 12n multiplies
      3n divides, 2n compares, and 4.5 adds for the original method
      
  */
  double *ox;
  double *fx;
  double cang;
  double thetad;
  double phid;
  double gammad;
  double g;
  double h;
  double p;
  double r;
  double u;
  double v;
  double pu;
  double pv;
  double *ap;
  double a[9];
  double alpha;
  double beta;

  int iu2;
  int iu3;

  int success;
  int inatom;

  int m;
  int k;

  int n;
  int dljpot_input_stream_spacer;

  int ivec_len;
  int i;

  char no_trans;
  char c1;
  char c2;
  char c3;
  char c4;
  char c5;
  char c6;
  char c7;

  FILE *ofp;
  FILE *efp;
  success = 1;
  ivec_len    = _IVEC_LEN_;
  iu2         = state->iu2;
  iu3         = state->iu3;
  inatom      = state->inatom;
  ofp         = state->ofp;
  ox          = state->ox;
  fx          = state->fx;
  dljpot_input_stream_spacer = 11 * _IVEC_LEN_;
  ap          = (double*)&a[0];
  no_trans    = 'N';
  if ((iu2 == 1) || (iu3 == 1)) {
    cang   = state->cang;
    thetad = theta*cang;
    phid   = phi*cang;
    gammad = gamma*cang;
    fprintf(ofp,"/n coordinates rotated by mobcal_rotate\n"
	    " theta = %le, phi = %le, gamma = %le\n",
	    thetad,phid,gammad);
  }
  /*
    let us label the combined rotation matrix a.
    as we store the coordinates in separate vectors we will
    want to for the transpose of a above. we will store a^t
    in column major order so as to be amenable to the dgemm call.

      x'"     gu-hpv  -gv-hpu  -hr     x
      y'"  = (hu+gpv  -hv+gpu   gr) *  y
      z'"     -rv       -ru     p      z

  */
  g = cos(gamma);
  h = sin(gamma);
  p = cos(phi);
  r = sin(phi);
  u = cos(theta);
  v = sin(theta);

  pu = p*u;
  pv = p*v;

  a[0] = g*u - h*pv;
  a[1] = -(g*v) - h*pu;
  a[2] = -h*r;

  a[3] = h*u + g*pv;
  a[4] = -h*v + g*pu;
  a[5] = g*r;

  a[6] = -r*v;
  a[7] = -r*u;
  a[8] = p;

  alpha = 1.0;
  beta  = 0.0;
  m     = state->vec_len;
  n     = 3;
  k     = 3;
  /*
  dgemm_(&no_trans,&no_trans,&m,&n,&k,&alpha,ox,&m,ap,&k,&beta,fx,&m);
  */
  fx     = state->dljpot_input_stream;
  for (i=0;i<m;i+=_IVEC_LEN_) {
    dgemm_(&no_trans,&no_trans,&ivec_len,&n,&k,&alpha,ox,&m,ap,&k,&beta,fx,&ivec_len);
    ox += _IVEC_LEN_;
    fx += dljpot_input_stream_spacer;
  }
  /*
  for (iatom=0;iatom<inatom;iatom++) {
    xx = ox[iatom];
    xy = oy[iatom];
    xz = oz[iatom];
    fx[iatom] = (at0 * xx) + (at1 * xy) + (at2 * xz);
    fy[iatom] = (at3 * xx) + (at4 * xy) + (at5 * xz);
    fz[iatom] = (at6 * xx) + (at7 * xy) + (at8 * xz);
  }
  */
  /* 
    The old broken way:
c
      if(iu2.eq.1.or.iu3.eq.1) write(8,610) theta*cang,phi*cang,
     ?gamma*cang
  610 format(//1x,'coordinates rotated by ROTATE',//1x,
     ?'theta=',1pe11.4,1x,'phi=',e11.4,1x,'gamma=',1pe11.4,/)
c

      do 1000 iatom=1,inatom
      rxy=dsqrt(ox(iatom)*ox(iatom)+(oy(iatom)*oy(iatom)))
      if(rxy.eq.0.d0) goto 1010
      otheta=dacos(ox(iatom)/rxy)
      if(oy(iatom).lt.0.d0) otheta=(2.d0*pi)-otheta
      ntheta=otheta+theta
 1010 fx(iatom)=dcos(ntheta)*rxy
 1000 fy(iatom)=dsin(ntheta)*rxy
c
      do 2000 iatom=1,inatom
      rzy=dsqrt(oz(iatom)*oz(iatom)+(fy(iatom)*fy(iatom)))
      if(rzy.eq.0.d0) goto 2010
      ophi=dacos(oz(iatom)/rzy)
      if(fy(iatom).lt.0.d0) ophi=(2.d0*pi)-ophi
      nphi=ophi+phi
 2010 fz(iatom)=dcos(nphi)*rzy
 2000 fy(iatom)=dsin(nphi)*rzy
c
      do 3000 iatom=1,inatom
      rxy=dsqrt(fx(iatom)*fx(iatom)+(fy(iatom)*fy(iatom)))
      if(rxy.eq.0.d0) goto 3010
      ogamma=dacos(fx(iatom)/rxy)
      if(fy(iatom).lt.0.d0) ogamma=(2.d0*pi)-ogamma
      ngamma=ogamma+gamma
 3010 fx(iatom)=dcos(ngamma)*rxy
 3000 fy(iatom)=dsin(ngamma)*rxy
c
      if(iu2.eq.0) goto 4000
      write(8,620)
  620 format(9x,'initial coordinates',24x,'new coordinates',/)
      do 4020 iatom=1,inatom
 4020 write(8,600) ox(iatom),oy(iatom),oz(iatom),fx(iatom),
     ?fy(iatom),fz(iatom)
  600 format(1x,1pe11.4,2(1x,e11.4),5x,3(1x,e11.4))
      close (10)
c
 4000 return
      end
  */
  if (iu2 == 1) {
    success = mobcal_print_rotated_coords(state);
  }
}
