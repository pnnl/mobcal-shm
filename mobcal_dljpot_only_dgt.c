#include "system_includes.h"
#include "mobcal_state_struct.h"
#include "blas.h"
#include "mobcal_vec_set.h"
#include "mobcal_dljpot_only_inner.h"
#include "mobcal_dljpot_only_dgt.h"
void mobcal_dljpot_only_dgt(struct mobcal_state_struct *state,
			double x, double y, double z,
			double *pot_p) {
  /*
    Routine to calucate the L-J + ion-dipole potential.
    Called by: mobcal_dljpot_only
    Calls:     exp, ddot_, mobcal_vec_set, mobcal_dljpot_only_inner;
    Uses:      inatom, xk, xeo, xe, four_pi, fx, fy, fz,
               eox4, ro6lj, ro12lj, and pcharge fields of state,
	       and pc_scale_pcharge scaling of pcharge.
    Sets       pot
  */
  double *px;
  double *py;
  double *pz;
  double *bndc;

  double *qpolvi;
  double *rx;
  double *ry;
  double *rz;
  double *qpolv_s0;
  double *e00v_s0;
  double *qpolv_s1;
  double *e00v_s1;
  double *qpolv_s2;
  double *e00v_s2;

  double *workspace;
	
  double qpoli;
  double r0;
  double r1;
  double r2;

  double qpol_s0;
  double e00_s0;

  double qpol_s1;
  double e00_s1;

  double qpol_s2;
  double e00_s2;

  double dipolxx_r0;
  double dipolxx_r1;
  double dipolxx_r2;
  double dipolzz_r0;
  double dipolzz_r1;
  double dipolzz_r2;

  double bond;
  double bond_c;
  double dipolxx;
  double dipolzz;
  double xe2_o_4pi_xeo;
  double xk;
  double xkt;
  double recip_xkt;
  double ptfn;
  double recip_ptfn;
  double temp_pot;
  double weight;
  double pot;

  double zero;
  double romax;
  double pot_min;
  double ones[_IVEC_LEN_];
  double part_ptfn[3];
  
  double pot_mol[3];

  int incx;
  int ivec_len;

  int avec_len;
  int i;

  int isamp;
  int inatom;

  int veclen;
  int work_space_len;

  zero    = 0.0;
  incx    = 1;
  romax = state->romax;
  avec_len = state->avec_len;
  ivec_len = _IVEC_LEN_;
  for (i=0;i<_IVEC_LEN_;i++) {
    ones[i] = 1.0;
  }
  /*
  xe2_o_4pi_xeo = state->xe2_o_4pi_xeo;
  dipolzz = (1.710e-30 * 0.5) * xe2_o_4pi_xeo;
  */
  dipolzz = state->dipol;
  dipolxx = dipolzz;
  pot_min  = 1.0e8;
  /*
    Input vectors.
  */
  xk      = state->xk;
  veclen  = state->vec_len;
  inatom  = state->inatom;
  /*
    local vectors:
  */
  workspace      = state->dljpot_workspace;
  work_space_len = 64*avec_len;
  mobcal_vec_set(work_space_len,workspace,zero);
  /*
    local vectors:
  */
  px = workspace;
  py = &px[avec_len];
  pz = &py[avec_len];
  bndc = &pz[avec_len];
  qpolvi = &bndc[avec_len];		
  rx      = &qpolvi[avec_len];
  ry      = &rx[avec_len];
  rz      = &ry[avec_len];

  qpolv_s0  = &rz[avec_len];	
  e00v_s0   = &qpolv_s0[avec_len];

  qpolv_s1  = &e00v_s0[avec_len];
  e00v_s1   = &qpolv_s1[avec_len];

  qpolv_s2  = &e00v_s1[avec_len];
  e00v_s2   = &qpolv_s2[avec_len];

  bond    = 1.0976e-10;
  bond_c  = bond * 0.5;

  for (i=0;i<_IVEC_LEN_;i++) {
    px[i] = x;
  }
  for (i=0;i<_IVEC_LEN_;i++) {
    py[i] = y;
  }
  for (i=0;i<_IVEC_LEN_;i++) {
    pz[i] = z;
  }
  for (i=0;i<_IVEC_LEN_;i++) {
    bndc[i] = bond_c;
  }

  /*
    Set the work space fields.
  */
  mobcal_dljpot_only_inner(state);
  /*
    This could be one large dgemm.
  */

  qpoli  = ddot_(&ivec_len,qpolvi,&incx,ones,&incx);
  r0     = ddot_(&ivec_len,rx,&incx,ones,&incx);
  r1     = ddot_(&ivec_len,ry,&incx,ones,&incx);
  r2     = ddot_(&ivec_len,rz,&incx,ones,&incx);

  qpol_s0 = ddot_(&ivec_len,qpolv_s0,&incx,ones,&incx);
  e00_s0   = ddot_(&ivec_len,e00v_s0,&incx,ones,&incx);

  qpol_s1 = ddot_(&ivec_len,qpolv_s1,&incx,ones,&incx);
  e00_s1   = ddot_(&ivec_len,e00v_s1,&incx,ones,&incx);

  qpol_s2 = ddot_(&ivec_len,qpolv_s2,&incx,ones,&incx);
  e00_s2   = ddot_(&ivec_len,e00v_s2,&incx,ones,&incx);

  /*
    isamp = 0.
    dpol[0] = dipolzz; dpol[1] = dipolxx dpol[2] = dipolxx
  */

  dipolzz_r0 = dipolzz * r0;
  dipolxx_r1 = dipolxx * r1;
  dipolxx_r2 = dipolxx * r2;
  pot_mol[0] = e00_s0 + qpoli + qpoli + qpol_s0 - (dipolzz_r0 * r0) -
    (dipolxx_r1 * r1) - (dipolxx_r2 * r2);
  if (pot_min > pot_mol[0]) {
    pot_min = pot_mol[0];
  }
  /*
    isamp = 1.
    dpol[0] = dipolxx; dpol[1] = dipolzz dpol[2] = dipolxx
  */
  dipolxx_r0 = dipolxx * r0;
  dipolzz_r1 = dipolzz * r1;
  pot_mol[1] = e00_s1 + qpoli + qpoli + qpol_s1 - (dipolxx_r0 * r0) -
    (dipolzz_r1 * r1) - (dipolxx_r2 * r2);
  if (pot_min > pot_mol[1]) {
    pot_min = pot_mol[1];
  }
  /*
    isamp = 2.
    dpol[0] = dipolxx; dpol[1] = dipolxx dpol[2] = dipolzz
  */
  dipolzz_r2 = dipolzz * r2;
  pot_mol[2] = e00_s2 + qpoli + qpoli + qpol_s2 - (dipolxx_r0 * r0) -
     (dipolxx_r1 * r1) - (dipolzz_r2 * r2);
  if (pot_min > pot_mol[2]) {
    pot_min = pot_mol[2];
  }
  xkt     = 500.0*xk;
  recip_xkt = 1.0/xkt;
  ptfn = 0.0;
  for (isamp = 0;isamp<3;isamp++) {
    temp_pot = exp((pot_min - pot_mol[isamp]) * recip_xkt);
    part_ptfn[isamp] = temp_pot;
    ptfn += temp_pot;
  }
  /*
  for (isamp=0;isamp<3;isamp++) {
    temp_pot = pot_mol[isamp] - pot_min;
    ptfn += ptfn + exp(-temp_pot*recip_xkt);
  }
  */
  recip_ptfn = 1.0/ptfn;
  pot = 0.0;
  for (isamp=0;isamp<3;isamp++) {
    /*
    temp_pot = pot_mol[isamp] - pot_min;
    weight   = exp(-temp_pot*recip_xkt)*recip_ptfn;
    */
    weight = part_ptfn[isamp]*recip_ptfn;
    pot      += (weight * pot_mol[isamp]);
  }    
  *pot_p   = pot;
}
