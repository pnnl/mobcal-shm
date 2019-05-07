#include "system_includes.h"
#include "mobcal_state_struct.h"
#include "blas.h"
#include "mobcal_vec_set.h"
#include "mobcal_dljpot_inner_iu.h"
#include "mobcal_dljpot_iu.h"
void mobcal_dljpot_iu(struct mobcal_state_struct *state,
		   double x, double y, double z,
		   double *pot_p,
		   double *dpotx_p, 
		   double *dpoty_p,
		   double *dpotz_p,
		   double *dmax_p) {
  /*
    Routine to calucate the L-J + ion-dipole potential via
    the IU algorithm
    Called by: mobcal_dljpot
    Calls:     exp, ddot_, mobcal_vec_set, mobcal_dljpot_inner_iu
    Uses:      inatom, xk, xeo, xe, romax, four_pi, coords,
               properties, (fx, fy, fz,
               eox4, ro6lj, ro12lj, dro6, dro12, and pcharge) fields of state,
	       and mthree_pcharge scalings of pcharge.
    Sets       pot, dpotx, dpoty, dpotz and dmax.	       
  */

  double *px;
  double *py;
  double *pz;

  double *rxv;
  double *ryv;
  double *rzv;
  double *sum1v;
  double *sum2v;
  double *sum3v;
  double *sum4v;
  double *sum5v;
  double *sum6v;
  double *e00v;
  double *de00xv;
  double *de00yv;
  double *de00zv;

  double *workspace;
	
  double dipol;
  double r0;
  double r1;
  double r2;
  double sum1;
  double sum2;
  double sum3;
  double sum4;
  double sum5;
  double sum6;
  double e00;

  double de00x;
  double de00y;
  double de00z;

  double pot;
  double dpotx;
  double dpoty;
  double dpotz;
  double dmax2;
  double recip_dmax2;

  double zero;
  double ones[_IVEC_LEN_];
  
  int incx;
  int ivec_len;

  int avec_len;
  int i;

  int work_space_len;

  zero    = 0.0;
  incx    = 1;
  avec_len = state->avec_len;
  ivec_len = _IVEC_LEN_;
  for (i=0;i<_IVEC_LEN_;i++) {
    ones[i] = 1.0;
  }
  dipol = state->dipol;
  /*
    local vectors:
  */
  workspace      = state->dljpot_workspace;
  /*
    See mobcal_alloc1 for allocation of this workspace.
  */
  work_space_len = avec_len * 64;
  mobcal_vec_set(work_space_len,workspace,zero);

  px = workspace;
  py = &px[avec_len];
  pz = &py[avec_len];

  rxv     = &pz[avec_len];
  ryv     = &rxv[avec_len];
  rzv     = &ryv[avec_len];
  sum1v   = &rzv[avec_len];
  sum2v   = &sum1v[avec_len];	
  sum3v   = &sum2v[avec_len];	
  sum4v   = &sum3v[avec_len];	
  sum5v   = &sum4v[avec_len];	
  sum6v   = &sum5v[avec_len];	

  e00v      = &sum6v[avec_len];
  de00xv    = &e00v[avec_len];
  de00yv    = &de00xv[avec_len];
  de00zv    = &de00yv[avec_len];

  for (i=0;i<_IVEC_LEN_;i++) {
    px[i] = x;
  }
  for (i=0;i<_IVEC_LEN_;i++) {
    py[i] = y;
  }
  for (i=0;i<_IVEC_LEN_;i++) {
    pz[i] = z;
  }
  /*
    Set the work space fields.
  */
  mobcal_dljpot_inner_iu(state, &recip_dmax2);
  /*
    This could be one large dgemm.
  */

  r0     = ddot_(&ivec_len,rxv,&incx,ones,&incx);
  r1     = ddot_(&ivec_len,ryv,&incx,ones,&incx);
  r2     = ddot_(&ivec_len,rzv,&incx,ones,&incx);
  sum1   = ddot_(&ivec_len,sum1v,&incx,ones,&incx);
  sum2   = ddot_(&ivec_len,sum2v,&incx,ones,&incx);
  sum3   = ddot_(&ivec_len,sum3v,&incx,ones,&incx);
  sum4   = ddot_(&ivec_len,sum4v,&incx,ones,&incx);
  sum5   = ddot_(&ivec_len,sum5v,&incx,ones,&incx);
  sum6   = ddot_(&ivec_len,sum6v,&incx,ones,&incx);

  e00    = ddot_(&ivec_len,e00v,&incx,ones,&incx);
  de00x  = ddot_(&ivec_len,de00xv,&incx,ones,&incx);
  de00y  = ddot_(&ivec_len,de00yv,&incx,ones,&incx);
  de00z  = ddot_(&ivec_len,de00zv,&incx,ones,&incx);

  pot = e00 - (dipol*((r0*r0) + (r1*r1) + (r2*r2)));
  dpotx = de00x - (2.0 * dipol * ((r0*sum1) + (r1*sum2) + (r2*sum3)));
  dpoty = de00y - (2.0 * dipol * ((r0*sum2) + (r1*sum4) + (r2*sum5)));
  dpotz = de00z - (2.0 * dipol * ((r0*sum3) + (r1*sum5) + (r2*sum6)));
		       
  *pot_p   = pot;
  *dpotx_p = dpotx;
  *dpoty_p = dpoty;
  *dpotz_p = dpotz;
  dmax2 = 1.0e0/recip_dmax2;
  *dmax_p  = sqrt(dmax2);
}
