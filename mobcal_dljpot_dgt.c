#include "system_includes.h"
#include "mobcal_state_struct.h"
#include "blas.h"
#include "mobcal_vec_set.h"
#include "mobcal_dljpot_inner.h"
#include "mobcal_max.h"
#include "mobcal_dljpot_dgt.h"
void mobcal_dljpot_dgt(struct mobcal_state_struct *state,
		       double x, double y, double z,
		       double *pot_p,
		       double *dpotx_p, 
		       double *dpoty_p,
		       double *dpotz_p,
		       double *dmax_p) {
  /*
    Routine to calucate the L-J + ion-dipole potential.
    Called by: mobcal_dljpot
    Calls:     exp, ddot_, mobcal_vec_set, mobcal_dljpot_inner, mobcal_max
    Uses:      inatom, xk, xeo, xe, romax, four_pi, coords,
               properties, (fx, fy, fz,
               eox4, ro6lj, ro12lj, dro6, dro12, and pcharge) fields of state,
	       and mthree_pcharge, and pc_scale_pcharge scalings of pcharge.
    Sets       pot, dpotx, dpoty, dpotz and dmax.	       
  */
  double *fx;
  double *fy;
  double *fz;
  double *pcharge;
  double *eox4;
  double *ro6lj;
  double *ro12lj;
  double *dro6;
  double *dro12;
  double *mthree_pcharge;
  double *pc_scale_pcharge;

  double *px;
  double *py;
  double *pz;
  double *bndc;

  double *qpolvi;
  double *dqpolxi;
  double *dqpolyi;
  double *dqpolzi;
  double *rx;
  double *ry;
  double *rz;
  double *sum1v;
  double *sum2v;
  double *sum3v;
  double *sum4v;
  double *sum5v;
  double *sum6v;
  double *qpolv_s0;
  double *dqpolx_s0;
  double *dqpoly_s0;
  double *dqpolz_s0;
  double *e00v_s0;
  double *de00x_s0;
  double *de00y_s0;
  double *de00z_s0;
  double *qpolv_s1;
  double *dqpolx_s1;
  double *dqpoly_s1;
  double *dqpolz_s1;
  double *e00v_s1;
  double *de00x_s1;
  double *de00y_s1;
  double *de00z_s1;
  double *qpolv_s2;
  double *dqpolx_s2;
  double *dqpoly_s2;
  double *dqpolz_s2;
  double *e00v_s2;
  double *de00x_s2;
  double *de00y_s2;
  double *de00z_s2;

  double *workspace;
	
  double qpoli;
  double dqpolx;
  double dqpoly;
  double dqpolz;
  double r0;
  double r1;
  double r2;
  double sum1;
  double sum2;
  double sum3;
  double sum4;
  double sum5;
  double sum6;

  double qpol_s0;
  double dqpolx_0;
  double dqpoly_0;
  double dqpolz_0;
  double e00_s0;
  double de00x_0;
  double de00y_0;
  double de00z_0;

  double qpol_s1;
  double dqpolx_1;
  double dqpoly_1;
  double dqpolz_1;
  double e00_s1;
  double de00x_1;
  double de00y_1;
  double de00z_1;

  double qpol_s2;
  double dqpolx_2;
  double dqpoly_2;
  double dqpolz_2;
  double e00_s2;
  double de00x_2;
  double de00y_2;
  double de00z_2;

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
  double xk;
  double xkt;
  double recip_xkt;
  double ptfn;
  double recip_ptfn;
  double temp_pot;
  double weight;
  double pot;
  double dpotx;
  double dpoty;
  double dpotz;
  double recip_dmax2;
  double dmax2;

  double zero;
  double romax;
  double pot_min;
  double ones[_IVEC_LEN_];
  double part_ptfn[3];
  
  double pot_mol[3];
  double dpotx_mol[3];
  double dpoty_mol[3];
  double dpotz_mol[3];

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
    ones[i] = 1.0e0;
  }
  dipolzz = state->dipol;
  dipolxx = dipolzz;
  pot_min  = 1.0e8;
  /*
    Input vectors.
  */
  fx      = state->fx;
  fy      = state->fy;
  fz      = state->fz;
  pcharge = state->pcharge;
  eox4    = state->eox4;
  ro6lj   = state->ro6lj;
  ro12lj  = state->ro12lj;
  dro6    = state->dro6;
  dro12   = state->dro12;
  xk      = state->xk;
  veclen  = state->vec_len;
  inatom  = state->inatom;
  mthree_pcharge = state->mthree_pcharge;
  pc_scale_pcharge = state->pc_scale_pcharge;
  /*
    local vectors:
  */
  workspace      = state->dljpot_workspace;
  work_space_len = 64*avec_len;
  mobcal_vec_set(work_space_len,workspace,zero);

  px = workspace;
  py = &px[avec_len];
  pz = &py[avec_len];
  bndc = &pz[avec_len];
  qpolvi = &bndc[avec_len];		
  dqpolxi = &qpolvi[avec_len];
  dqpolyi = &dqpolxi[avec_len];
  dqpolzi = &dqpolyi[avec_len];
  rx      = &dqpolzi[avec_len];
  ry      = &rx[avec_len];
  rz      = &ry[avec_len];
  sum1v   = &rz[avec_len];
  sum2v   = &sum1v[avec_len];	
  sum3v   = &sum2v[avec_len];	
  sum4v   = &sum3v[avec_len];	
  sum5v   = &sum4v[avec_len];	
  sum6v   = &sum5v[avec_len];	

  qpolv_s0  = &sum6v[avec_len];	
  dqpolx_s0 = &qpolv_s0[avec_len];
  dqpoly_s0 = &dqpolx_s0[avec_len];
  dqpolz_s0 = &dqpoly_s0[avec_len];
  e00v_s0   = &dqpolz_s0[avec_len];
  de00x_s0  = &e00v_s0[avec_len];
  de00y_s0  = &de00x_s0[avec_len];
  de00z_s0  = &de00y_s0[avec_len];

  qpolv_s1  = &de00z_s0[avec_len];	
  dqpolx_s1 = &qpolv_s1[avec_len];
  dqpoly_s1 = &dqpolx_s1[avec_len];
  dqpolz_s1 = &dqpoly_s1[avec_len];
  e00v_s1   = &dqpolz_s1[avec_len];
  de00x_s1  = &e00v_s1[avec_len];
  de00y_s1  = &de00x_s1[avec_len];
  de00z_s1  = &de00y_s1[avec_len];

  qpolv_s2  = &de00z_s1[avec_len];	
  dqpolx_s2 = &qpolv_s2[avec_len];
  dqpoly_s2 = &dqpolx_s2[avec_len];
  dqpolz_s2 = &dqpoly_s2[avec_len];
  e00v_s2   = &dqpolz_s2[avec_len];
  de00x_s2  = &e00v_s2[avec_len];
  de00y_s2  = &de00x_s2[avec_len];
  de00z_s2  = &de00y_s2[avec_len];

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
  mobcal_dljpot_inner(state, &recip_dmax2);
  /*
    This could be one large dgemm.
  */

  qpoli  = ddot_(&ivec_len,qpolvi,&incx,ones,&incx);
  dqpolx = ddot_(&ivec_len,dqpolxi,&incx,ones,&incx);
  dqpoly = ddot_(&ivec_len,dqpolyi,&incx,ones,&incx);
  dqpolz = ddot_(&ivec_len,dqpolzi,&incx,ones,&incx);
  r0     = ddot_(&ivec_len,rx,&incx,ones,&incx);
  r1     = ddot_(&ivec_len,ry,&incx,ones,&incx);
  r2     = ddot_(&ivec_len,rz,&incx,ones,&incx);
  sum1   = ddot_(&ivec_len,sum1v,&incx,ones,&incx);
  sum2   = ddot_(&ivec_len,sum2v,&incx,ones,&incx);
  sum3   = ddot_(&ivec_len,sum3v,&incx,ones,&incx);
  sum4   = ddot_(&ivec_len,sum4v,&incx,ones,&incx);
  sum5   = ddot_(&ivec_len,sum5v,&incx,ones,&incx);
  sum6   = ddot_(&ivec_len,sum6v,&incx,ones,&incx);

  qpol_s0 = ddot_(&ivec_len,qpolv_s0,&incx,ones,&incx);
  dqpolx_0 = ddot_(&ivec_len,dqpolx_s0,&incx,ones,&incx);
  dqpoly_0 = ddot_(&ivec_len,dqpoly_s0,&incx,ones,&incx);
  dqpolz_0 = ddot_(&ivec_len,dqpolz_s0,&incx,ones,&incx);
  e00_s0   = ddot_(&ivec_len,e00v_s0,&incx,ones,&incx);
  de00x_0  = ddot_(&ivec_len,de00x_s0,&incx,ones,&incx);
  de00y_0  = ddot_(&ivec_len,de00y_s0,&incx,ones,&incx);
  de00z_0  = ddot_(&ivec_len,de00z_s0,&incx,ones,&incx);

  qpol_s1 = ddot_(&ivec_len,qpolv_s1,&incx,ones,&incx);
  dqpolx_1 = ddot_(&ivec_len,dqpolx_s1,&incx,ones,&incx);
  dqpoly_1 = ddot_(&ivec_len,dqpoly_s1,&incx,ones,&incx);
  dqpolz_1 = ddot_(&ivec_len,dqpolz_s1,&incx,ones,&incx);
  e00_s1   = ddot_(&ivec_len,e00v_s1,&incx,ones,&incx);
  de00x_1  = ddot_(&ivec_len,de00x_s1,&incx,ones,&incx);
  de00y_1  = ddot_(&ivec_len,de00y_s1,&incx,ones,&incx);
  de00z_1  = ddot_(&ivec_len,de00z_s1,&incx,ones,&incx);

  qpol_s2 = ddot_(&ivec_len,qpolv_s2,&incx,ones,&incx);
  dqpolx_2 = ddot_(&ivec_len,dqpolx_s2,&incx,ones,&incx);
  dqpoly_2 = ddot_(&ivec_len,dqpoly_s2,&incx,ones,&incx);
  dqpolz_2 = ddot_(&ivec_len,dqpolz_s2,&incx,ones,&incx);
  e00_s2   = ddot_(&ivec_len,e00v_s2,&incx,ones,&incx);
  de00x_2  = ddot_(&ivec_len,de00x_s2,&incx,ones,&incx);
  de00y_2  = ddot_(&ivec_len,de00y_s2,&incx,ones,&incx);
  de00z_2  = ddot_(&ivec_len,de00z_s2,&incx,ones,&incx);

  /*
    isamp = 0.
    dpol[0] = dipolzz; dpol[1] = dipolxx; dpol[2] = dipolxx
  */

  dipolzz_r0 = dipolzz * r0;
  dipolxx_r1 = dipolxx * r1;
  dipolxx_r2 = dipolxx * r2;
  pot_mol[0] = e00_s0 + qpoli + qpoli + qpol_s0 - (dipolzz_r0 * r0) -
    (dipolxx_r1 * r1) - (dipolxx_r2 * r2);
  if (pot_min > pot_mol[0]) {
    pot_min = pot_mol[0];
  }

  dpotx_mol[0] = (de00x_0 + dqpolx + dqpolx + dqpolx_0) -
   2.0 * ((dipolzz_r0 * sum1) + (dipolxx_r1 * sum2) + (dipolxx_r2 * sum3));

  dpoty_mol[0] = (de00y_0 + dqpoly + dqpoly + dqpoly_0) -
   2.0 * ((dipolzz_r0 * sum2) + (dipolxx_r1 * sum4) + (dipolxx_r2 *sum5));

  dpotz_mol[0] = (de00z_0 + dqpolz + dqpolz + dqpolz_0) -
   2.0 * ((dipolzz_r0 * sum3) + (dipolxx_r1 * sum5) + (dipolxx_r2 *sum6));

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
  dpotx_mol[1] = (de00x_1 + dqpolx + dqpolx + dqpolx_1) -
   2.0 * ((dipolxx_r0 * sum1) + (dipolzz_r1 * sum2) + (dipolxx_r2 * sum3));

  dpoty_mol[1] = (de00y_1 + dqpoly + dqpoly + dqpoly_1) -
   2.0 * ((dipolxx_r0 * sum2) + (dipolzz_r1 * sum4) + (dipolxx_r2 *sum5));

  dpotz_mol[1] = (de00z_1 + dqpolz + dqpolz + dqpolz_1) -
   2.0 * ((dipolxx_r0 * sum3) + (dipolzz_r1 * sum5) + (dipolxx_r2 *sum6));

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

  dpotx_mol[2] = (de00x_2 + dqpolx + dqpolx + dqpolx_2) -
   2.0 * ((dipolxx_r0 * sum1) + (dipolxx_r1 * sum2) + (dipolzz_r2 * sum3));

  dpoty_mol[2] = (de00y_2 + dqpoly + dqpoly + dqpoly_2) -
   2.0 * ((dipolxx_r0 * sum2) + (dipolxx_r1 * sum4) + (dipolzz_r2 *sum5));

  dpotz_mol[2] = (de00z_2 + dqpolz + dqpolz + dqpolz_2) -
   2.0 * ((dipolxx_r0 * sum3) + (dipolxx_r1 * sum5) + (dipolzz_r2 *sum6));


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
  dpotx = 0.0;
  dpoty = 0.0;
  dpotz = 0.0;
  for (isamp=0;isamp<3;isamp++) {
    /*
    temp_pot = pot_mol[isamp] - pot_min;
    weight   = exp(-temp_pot*recip_xkt)*recip_ptfn;
    */
    weight = part_ptfn[isamp]*recip_ptfn;
    pot      += (weight * pot_mol[isamp]);
    dpotx    += (weight * dpotx_mol[isamp]);
    dpoty    += (weight * dpoty_mol[isamp]);
    dpotz    += (weight * dpotz_mol[isamp]);
  }    
  *pot_p   = pot;
  *dpotx_p = dpotx;
  *dpoty_p = dpoty;
  *dpotz_p = dpotz;
  dmax2    = 1.0e0/recip_dmax2;
  *dmax_p  = sqrt(dmax2);
}
