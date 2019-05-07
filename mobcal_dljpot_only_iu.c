#include "system_includes.h"
#include "mobcal_state_struct.h"
#include "blas.h"
#include "mobcal_vec_set.h"
#include "mobcal_dljpot_only_inner_iu.h"
#include "mobcal_dljpot_only_iu.h"
void mobcal_dljpot_only_iu(struct mobcal_state_struct *state,
			double x, double y, double z,
			double *pot_p) {
  /*
    Routine to calucate the L-J + ion-dipole potential.
    Called by: mobcal_rmax_emax_r00, mobcal_gsang
    Calls:     exp, ddot_, mobcal_vec_set, mobcal_dljpot_only_inner_iu;
    Uses:      inatom, xk, xeo, xe, four_pi, fx, fy, fz,
               eox4, ro6lj, ro12lj, and pcharge fields of state,
	       and pc_scale_pcharge scaling of pcharge.
    Sets       pot
  */
  double *px;
  double *py;
  double *pz;

  double *rx;
  double *ry;
  double *rz;
  double *e00v;

  double *workspace;
	
  double r0;
  double r1;
  double r2;
  double e00;

  double dipol;
  double pot;

  double zero;
  double ones[_IVEC_LEN_];
  
  int incx;
  int ivec_len;

  int avec_len;
  int i;

  int work_space_len;
  int padi;

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
  work_space_len = avec_len * 64;
  mobcal_vec_set(work_space_len,workspace,zero);
  /*
    local vectors:
  */
  px = workspace;
  py = &px[avec_len];
  pz = &py[avec_len];
  rx      = &pz[avec_len];
  ry      = &rx[avec_len];
  rz      = &ry[avec_len];
  e00v    = &rz[avec_len];	

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
  mobcal_dljpot_only_inner_iu(state);
  /*
    This could be one large dgemm.
  */
  r0     = ddot_(&ivec_len,rx,&incx,ones,&incx);
  r1     = ddot_(&ivec_len,ry,&incx,ones,&incx);
  r2     = ddot_(&ivec_len,rz,&incx,ones,&incx);

  e00   = ddot_(&ivec_len,e00v,&incx,ones,&incx);

  pot    = e00 - (dipol*((r0*r0) + (r1*r1) + (r2*r2)));
  *pot_p   = pot;
}
