#include "system_includes.h"
#include "mobcal_state_struct.h"
/*
#include "mobcal_amax.h"
*/
#include "mobcal_dljpot_inner.h"
void mobcal_dljpot_inner(struct mobcal_state_struct *state, 
			 double *recip_dmax2_p) {
  /*
    Workspace needs to be 328 + 6* vec_len  elements long.
    the px, py, pz, and bndc vectors are assumed to be set by 
    the calling routine and are in the state->dljpot_workspace vector.
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


  /*
    Workspace vectors.
  */

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
	
  /*
  double *recip_rxyz2_00;
  double *recip_rxyz2_01;
  double *recip_rxyz2_10;
  double *recip_rxyz2_11;
  double *recip_rxyz2_20;
  double *recip_rxyz2_21;
  */

  double yregister[128];
  double *xx_center;
  double *yy_center;
  double *zz_center;
  double *xx_center2;
  double *yy_center2;
  double *zz_center2;
  double *rxyz_center2;
  double *recip_rxyz_center2;
  double *recip_rxyz_center;
  double *recip_rxyz_center3;
  double *lpc_scale_pcharge;
  double *lqpolvi;
  double *vptemp1;
  double *vptemp2;
  double *lpcharge;
  double *ldqpolx;
  double *ldqpoly;
  double *ldqpolz;
  double *rxyz3i;
  double *lrx;
  double *lry;
  double *lrz;
  double *recip_rxyz_center5;
  double *lmthree_pcharge;
  double *rxyz5i;
  double *lsum1v;
  double *lsum2v;
  double *lsum3v;
  double *lsum4v;
  double *lsum5v;
  double *lsum6v;
  double *pconv;
  double *lro12lj;
  double *lro6lj;
  double *leox4;
  double *ldro6;
  double *ldro12;
  double *lbndc;
  double *xx;
  double *yy;
  double *zz;
  double *rxyz2;
  double *recip_rxyz2;
  double *recip_rxyz;
  double *recip_rxyz3;
  double *recip_rxyz6;
  double *recip_rxyz12;
  double *recip_rxyz8;
  double *recip_rxyz14;
  double *de00g;
  double *rv;
  double *workspace;

  double romax;
  double dmax;
  double recip_dmax;
  double recip_dmax2;

  int64_t vec_len;
  int avec_len;
  int iatom;

  int ivec_len;
  int i;

  int dmax_offset;
  int dljpot_input_stream_spacing;
  dljpot_input_stream_spacing = 11 * _IVEC_LEN_;
  /*
    Input vectors.
  fx      = state->fx;
  fy      = state->fy;
  fz      = state->fz;
  pcharge = state->pcharge;
  eox4    = state->eox4;
  ro6lj   = state->ro6lj;
  ro12lj  = state->ro12lj;
  dro6    = state->dro6;
  dro12   = state->dro12;
  mthree_pcharge = state->mthree_pcharge;
  pc_scale_pcharge = state->pc_scale_pcharge;
  */
  fx = state->dljpot_input_stream;
  /*
    The follwing set of pointer assignments uses address arithmetic.
  */
  fy = fx + _IVEC_LEN_;
  fz = fy + _IVEC_LEN_;
  pc_scale_pcharge = fz + _IVEC_LEN_;
  pcharge          = pc_scale_pcharge + _IVEC_LEN_;
  mthree_pcharge   = pcharge + _IVEC_LEN_;
  eox4             = mthree_pcharge + _IVEC_LEN_;
  ro6lj            = eox4 + _IVEC_LEN_;
  ro12lj           = ro6lj + _IVEC_LEN_;
  dro6             = ro12lj + _IVEC_LEN_;
  dro12            = dro6 + _IVEC_LEN_;

  romax   = state->romax;
  dmax    = 2.0 * romax;
  recip_dmax = 1.0/dmax;
  recip_dmax2 = recip_dmax * recip_dmax;


  dmax_offset = (int)state->dmax_offset;
  vec_len  = state->vec_len;
  workspace = state->dljpot_workspace;
  avec_len = state->avec_len;
  ivec_len = _IVEC_LEN_;
  /*
    local vectors:
  */
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


  /*
  recip_rxyz2_00   = &workspace[dmax_offset];
  recip_rxyz2_01   = &recip_rxyz2_00[vec_len];
  recip_rxyz2_10   = &recip_rxyz2_01[vec_len];
  recip_rxyz2_11   = &recip_rxyz2_10[vec_len];
  recip_rxyz2_20   = &recip_rxyz2_11[vec_len];
  recip_rxyz2_21   = &recip_rxyz2_20[vec_len];
  */
  /*
    Start of local ymm registers.
  */
  /* ymm1 */
  lro12lj   = &yregister[0];

  /* ymm2 */
  lro6lj    = &lro12lj[avec_len];
  /* ymm3 */
  leox4     = &lro6lj[avec_len];
  /* ymm4 */		  
  xx_center = &leox4[avec_len];
  /* ymm5 */
  yy_center = &xx_center[avec_len];
  /* ymm6 */
  zz_center = &yy_center[avec_len];
  /* ymm7 */
  xx_center2 = &zz_center[avec_len];
  ldro6      = xx_center2;
  /* ymm8 */
  yy_center2 = &xx_center2[avec_len];
  ldro12     = yy_center2;
  /* ymm9 */
  zz_center2 = &yy_center2[avec_len];
  lbndc       = zz_center2;
  /* ymm10 */
  rxyz_center2 = &zz_center2[avec_len];
  recip_rxyz_center2  = rxyz_center2;
  vptemp2             = rxyz_center2;
  recip_rxyz12        = vptemp2;
  /* ymm11 */
  recip_rxyz_center = &rxyz_center2[avec_len];
  rxyz3i             = recip_rxyz_center;
  pconv              = recip_rxyz_center;
  xx                 = recip_rxyz_center;
  yy                 = recip_rxyz_center;
  zz                 = recip_rxyz_center;
  /* ymm12 */
  lpc_scale_pcharge = &recip_rxyz_center[avec_len];
  /* ymm13 */
  recip_rxyz_center3 = &lpc_scale_pcharge[avec_len];
  lqpolvi            = recip_rxyz_center3;
  recip_rxyz_center5 = recip_rxyz_center3;	
  rxyz5i             = recip_rxyz_center3;
  rxyz2              = recip_rxyz_center3;
  recip_rxyz2        = recip_rxyz_center3;
  recip_rxyz8        = recip_rxyz_center3;
  de00g              = recip_rxyz_center3;
  /* ymm14 */
  vptemp1            = &recip_rxyz_center3[avec_len];
  lpcharge           = vptemp1;
  recip_rxyz         = vptemp1;
  recip_rxyz3        = vptemp1;
  recip_rxyz6        = vptemp1;
  recip_rxyz14       = vptemp1;
  /* ymm15 */
  rv                 = &vptemp1[avec_len];
  lmthree_pcharge    = rv;
  ldqpolx            = rv;
  ldqpoly            = rv;
  ldqpolz            = rv;
  lrx                = rv;
  lry                = rv;
  lrz                = rv;
  lsum1v	     = rv;
  lsum2v	     = rv;
  lsum3v	     = rv;
  lsum4v	     = rv;
  lsum5v	     = rv;
  lsum6v	     = rv;

  for (iatom=0;iatom < vec_len;iatom += _IVEC_LEN_) {
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_;i++) {
      xx_center[i] = px[i] - fx[i];
    }
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_;i++) {
      yy_center[i] = py[i] - fy[i];
    }
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_;i++) {
      zz_center[i] = pz[i] - fz[i];
    }
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_;i++) {
      xx_center2[i] = xx_center[i] * xx_center[i];
    }
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_;i++) {
      yy_center2[i] = yy_center[i] * yy_center[i];
    }
    /*
      rxyz2 = xx_center2 + yy_center2.
    */
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_;i++) {
      recip_rxyz_center2[i] = xx_center2[i] + yy_center2[i];
    }
    /*
      zz_center2 = zz_center .* zz_center
    */
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_;i++) {
      zz_center2[i] = zz_center[i] * zz_center[i];
    }
    /*
      rxyz2 += zz_center2
    */
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_;i++) {
      recip_rxyz_center2[i] = rxyz_center2[i] + zz_center2[i];
    }
    /*
      recip_rxyz_center2 = 1.0 ./ rxyz_center2;
      Here we take advantage of the face that recip_rxyz_center2 points rxyz_center2.
      mobcal_pg1dq(ivec_len,recip_rxyz_center2,rxyz_center2);
    */
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      recip_rxyz_center2[i] = 1.0/recip_rxyz_center2[i];
    }
    /*
      recip_rxyz_center <- sqrt(recip_rxyz_center2)
    */
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      recip_rxyz_center[i] = sqrt(recip_rxyz_center2[i]);
    }
    /*
      Load qpolvi
    */
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      lqpolvi[i] = qpolvi[i];
    }
    /*
      Load pc_scale_pcharge.
    */
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      lpc_scale_pcharge[i] = pc_scale_pcharge[i];
    }
    /*
      qpolvi -= (recip_rxyz_center .* pc_scale_pcharge)
    */
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      lqpolvi[i] = lqpolvi[i] - (recip_rxyz_center[i] * lpc_scale_pcharge[i]);
    }
    /*
      Store qpolvi.
    */
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      qpolvi[i] = lqpolvi[i];
    }
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      recip_rxyz_center3[i] = recip_rxyz_center[i] * recip_rxyz_center2[i];
    }
    /* 
      vptemp1 <- pc_scale_pcharge .* recip_rxyz_center3
    */
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      vptemp1[i] = recip_rxyz_center3[i] * lpc_scale_pcharge[i];
    }
    /*
      load dqpolxi
    */
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      rv[i] = dqpolxi[i];
    }
    /*
      dqpolxi += vptemp1 .* xx_center
    */
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      rv[i] = rv[i] + (vptemp1[i] * xx_center[i]);
    }
    /*
      store dqpolxi  (dqpolxi <- rv);
    */
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      dqpolxi[i] = rv[i];
    }

    /*
      load dqpolyi
    */
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      rv[i] = dqpolyi[i];
    }
    /*
      dqpolyi += vptemp1 .* yy_center
    */
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      rv[i] = rv[i] + (vptemp1[i] * yy_center[i]);
    }
    /*
      store dqpolyi  (dqpolyi <- rv);
    */
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      dqpolyi[i] = rv[i];
    }

    /*
      load dqpolzi
    */
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      rv[i] = dqpolzi[i];
    }
    /*
      dqpolzi += vptemp1 .* zz_center
    */
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      rv[i] = rv[i] + (vptemp1[i] * zz_center[i]);
    }
    /*
      store dqpolzi  (dqpolzi <- rv);
    */
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      dqpolzi[i] = rv[i];
    }
    /*
      Load pcharge.
    */
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      lpcharge[i] = pcharge[i];
    }
    /*
      rxyz3i <- pcharge .* recip_rxyz_center3.
    */
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      rxyz3i[i] = lpcharge[i] * recip_rxyz_center3[i];
    }

    /*
      Load rx.
    */
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      rv[i] = rx[i];
    }
    /*
      rx += rxyz3i .* xx_center.
    */
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      rv[i] += rxyz3i[i] * xx_center[i];
    }
    /*
      Store rx.
    */
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      rx[i] = rv[i];
    }

    /*
      Load ry.
    */
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      rv[i] = ry[i];
    }
    /*
      ry += rxyz3i .* yy_center.
    */
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      rv[i] += rxyz3i[i] * yy_center[i];
    }
    /*
      Store ry.
    */
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      ry[i] = rv[i];
    }

    /*
      Load rz.
    */
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      rv[i] = rz[i];
    }
    /*
      rz += rxyz3i .* zz_center.
    */
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      rv[i] += rxyz3i[i] * zz_center[i];
    }
    /*
      Store rz.
    */
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      rz[i] = rv[i];
    }
    /*
      recip_rxyz_center5 <- recip_rxyz_center3 * recip_rxyz_center2.
      Here we make use of &recip_rxyz_center5 = &recip_rxyz_center3.
    */
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      recip_rxyz_center5[i] = recip_rxyz_center5[i] * recip_rxyz_center2[i];
    }
    /*
      Load mthree_pcharge
    */
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      lmthree_pcharge[i] = mthree_pcharge[i];
    }
    /*
      rxyz5i = mthree_pccharge .* recip_rxyx_center5.
      Here we take advantage of &rxyz5i = &recip_rxyz_center5.
    */
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      rxyz5i[i] = rxyz5i[i] * lmthree_pcharge[i];
    }
    /*
      Load sum1v
    */
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      rv[i] = sum1v[i];
    }
    /*
      sum1v += rxyz3i;
    */
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      rv[i] += rxyz3i[i];
    }
    /*
      sum1v += rxyz5i .* xx_center2;
    */
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      rv[i] += (rxyz5i[i] * xx_center2[i]);
    }
    /*
      Store sum1v
    */
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      sum1v[i] = rv[i];
    }

    /*
      Load sum4v
    */
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      rv[i] = sum4v[i];
    }
    /*
      sum4v += rxyz3i;
    */
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      rv[i] += rxyz3i[i];
    }
    /*
      sum4v += rxyz5i .* yy_center2;
    */
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      rv[i] += (rxyz5i[i] * yy_center2[i]);
    }
    /*
      Store sum4v
    */
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      sum4v[i] = rv[i];
    }

    /*
      Load sum6v
    */
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      rv[i] = sum6v[i];
    }
    /*
      sum6v += rxyz3i;
    */
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      rv[i] += rxyz3i[i];
    }
    /*
      sum6v += rxyz5i .* zz_center2;
    */
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      rv[i] += (rxyz5i[i] * zz_center2[i]);
    }
    /*
      Store sum6v
    */
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      sum6v[i] = rv[i];
    }

    /*
      pconv <- xx_center .* yy_center
    */
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      pconv[i] = xx_center[i] * yy_center[i];
    }
    /*
      Load sum2v
    */
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      rv[i] = sum2v[i];
    }
    /*
      sum2v += pconv .* rxyz5i
    */
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      rv[i] += (pconv[i] * rxyz5i[i]);
    }
    /*
      Store sum2v
    */
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      sum2v[i] = rv[i];
    }

    /*
      pconv <- xx_center .* zz_center
    */
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      pconv[i] = xx_center[i] * zz_center[i];
    }
    /*
      Load sum3v
    */
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      rv[i] = sum3v[i];
    }
    /*
      sum3v += pconv .* rxyz5i
    */
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      rv[i] += (pconv[i] * rxyz5i[i]);
    }
    /*
      Store sum3v
    */
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      sum3v[i] = rv[i];
    }

    /*
      pconv <- yy_center .* zz_center
    */
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      pconv[i] = yy_center[i] * zz_center[i];
    }
    /*
      Load sum5v
    */
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      rv[i] = sum5v[i];
    }
    /*
      sum5v += pconv .* rxyz5i
    */
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      rv[i] += (pconv[i] * rxyz5i[i]);
    }
    /*
      Store sum5v
    */
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      sum5v[i] = rv[i];
    }

    /*
      Load ro12lj
    */
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      lro12lj[i] = ro12lj[i];
    }
    /*
      Load ro6lj
    */
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      lro6lj[i]  = ro6lj[i];
    }
    /*
      Load eox4
    */
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      leox4[i]   = eox4[i];
    }
    /*
      Load dro6
    */
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      ldro6[i]   = dro6[i];
    }
    /*
      Load dro12
    */
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      ldro12[i]  = dro12[i];
    }
    /*
      Load bndc
    */
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      lbndc[i]   = bndc[i];
    }
    /*
      Now for the unrolling of the nested isamp,ibatom loops.
      The computations above are shared by all six iterations of the
      unrolling.
    */
    /*
      isamp = 0, ibatom = 0. 
      xx = xx_center - bndc, yy = yy_center, zz = zz_center
    */
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      xx[i] = xx_center[i] - lbndc[i];
    }
    /*
      rxyz2 <- xx .* xx
    */
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      recip_rxyz2[i] = xx[i] * xx[i];
    }
    /*
      rxyz2 += yy .* yy
    */
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      recip_rxyz2[i] += (yy_center[i] * yy_center[i]);
    }
    /*
      rxyz2 += zz .* zz
    */
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      recip_rxyz2[i] += (zz_center[i] * zz_center[i]);
    }
    /*
      recip_rxyz2 = 1.0/recip_rxyz2.
    */
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      recip_rxyz2[i] = 1.0/recip_rxyz2[i];
    }
    /*
      track recip_dmax2.
    recip_dmax2 = mobcal_amax(recip_rxyz2,recip_dmax2);
    */
    for (i=0;i< _IVEC_LEN_ ;i++) {
      if (recip_rxyz2[i] > recip_dmax2) {
	recip_dmax2 = recip_rxyz2[i];
      }
    }
    /*
      Store recip_rxyz2 for dmax computation in mobcal_dljpot.
    */
    /*
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      recip_rxyz2_00[i] = recip_rxyz2[i];
    }
    recip_rxyz2_00 += _IVEC_LEN_; // Caution address arithmetic.
    */
    /*
      recip_rxyz = sqrt(recip_rxyz2).
    */
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      recip_rxyz[i] = sqrt(recip_rxyz2[i]);
    }
    /* 
      Load qpolv_s0.
    */
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      rv[i] = qpolv_s0[i];
    }
    /*
      qpolv_s0 += recip_rxyz .* pc_scale_pcharge.
    */
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      rv[i] += (recip_rxyz[i] * lpc_scale_pcharge[i]);
    }
    /*
      Store qpolv_s0
    */
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      qpolv_s0[i] = rv[i];
    }
    /*
      recip_rxyz3 = recip_rxyz2 .* recip_rxyz;
      Here we take advantage of &recip_rxyz3 = recip_rxyz.
    */
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      recip_rxyz3[i] = recip_rxyz3[i] * recip_rxyz2[i];
    }
    /*
      vptemp2 <- recip_rxyz3 .* pc_scale_pcharge.
    */
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      vptemp2[i] = recip_rxyz3[i] * lpc_scale_pcharge[i];
    }

    /*
      Load dqpolx_s0
    */
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      rv[i] = dqpolx_s0[i];
    }
    /*
      dqpolx_s0 -= vptemp2 .* xx
    */
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      rv[i] -= (vptemp2[i] * xx[i]);
    }
    /*
      Store dqpolx_s0
    */
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      dqpolx_s0[i] = rv[i];
    }

    /*
      Load dqpoly_s0
    */
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      rv[i] = dqpoly_s0[i];
    }
    /*
      dqpoly_s0 -= vptemp2 .* yy
    */
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      rv[i] -= (vptemp2[i] * yy_center[i]);
    }
    /*
      Store dqpoly_s0
    */
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      dqpoly_s0[i] = rv[i];
    }
    
    /*
      Load dqpolz_s0
    */
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      rv[i] = dqpolz_s0[i];
    }
    /*
      dqpolz_s0 -= vptemp2 .* zz
    */
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      rv[i] -= (vptemp2[i] * zz_center[i]);
    }
    /*
      Store dqpolz_s0
    */
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      dqpolz_s0[i] = rv[i];
    }
    /*
      recip_rxyz6 <- recip_rxyz3 .* recip_rxyz3
      Here we use the fact that &recip_rxyz6 = &recip_rxyz3
    */
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      recip_rxyz6[i] = recip_rxyz6[i] * recip_rxyz6[i];
    }
    /*
      recip_rxyz12 = recip_rxyz6 .* recip_rxyz6
    */
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      recip_rxyz12[i] = recip_rxyz6[i] * recip_rxyz6[i];
    }
    /*
      recip_rxyz12 <- recip_rxyz12 .* ro12lj
    */
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      recip_rxyz12[i] = recip_rxyz12[i] * lro12lj[i];
    }
    /*
      recip_rxyz12 -= recip_rxyz6 .* ro6lj
    */
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      recip_rxyz12[i] -= (recip_rxyz6[i] * lro6lj[i]);
    }
    /*
      Load e00v_s0.
    */
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      rv[i] = e00v_s0[i];
    }
    /*
      e00v_s0 += eox4 .* (recip_rxyz12 .* ro12lj - recip_rxyz .* ro6.j)
      where the expresion in parenthesis has now ben stored in recip_rxyz12.
    */
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      rv[i] += (leox4[i] * recip_rxyz12[i]);
    }
    /*
      Store e00v_s0
    */
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      e00v_s0[i] = rv[i];
    }
    /*
      recip_rxyz8 = recip_rxyz6 .* recip_rxyz2.
      Here we make use of &recip_rxyz8 = &recip_rxyz2
    */
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      recip_rxyz8[i] = recip_rxyz8[i] * recip_rxyz6[i];
    }
    /*
      recip_rxyz14 <- recip_rxyz6 .* recip_rxyz8.
      Here we make use of &recip_rxyz14 = &recip_rxyz6
    */
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      recip_rxyz14[i] = recip_rxyz14[i] * recip_rxyz8[i];
    }
    /*
      de00g <- dro6 .* recip_rxyz8
      Now we make use of  &de00g = &recip_rxyz8
    */
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      de00g[i] = de00g[i] * ldro6[i];
    }
    /*
      de00g -= dro12 .* recip_rxyz14
    */
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      de00g[i] -= (recip_rxyz14[i] * ldro12[i]);
    }
    /*
      de00g = eox4 .* (dro6 .* recip_rxyz8 - dro12 .* recip_rxyz14).
    */
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      de00g[i] = de00g[i] * leox4[i];
    }
    /*
      Load de00x_s0
    */
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      rv[i] = de00x_s0[i];
    }
    /*
      de00x_s0 += de00g .* xx
    */
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      rv[i] += (de00g[i] * xx[i]);
    }
    /*
      Store de00x_s0
    */
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      de00x_s0[i] = rv[i];
    }

    /*
      Load de00y_s0
    */
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      rv[i] = de00y_s0[i];
    }
    /*
      de00y_s0 += de00g .* yy
    */
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      rv[i] += (de00g[i] * yy_center[i]);
    }
    /*
      Store de00y_s0
    */
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      de00y_s0[i] = rv[i];
    }

    /*
      Load de00z_s0
    */
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      rv[i] = de00z_s0[i];
    }
    /*
      de00z_s0 += de00g .* zz
    */
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      rv[i] += (de00g[i] * zz_center[i]);
    }
    /*
      Store de00z_s0
    */
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      de00z_s0[i] = rv[i];
    }
    /*
      isamp = 0, ibatom = 1. 
      xx = xx_center + bndc, yy = yy_center, zz = zz_center
    */
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      xx[i] = xx_center[i] + lbndc[i];
    }
    /*
      rxyz2 <- xx .* xx
    */
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      recip_rxyz2[i] = xx[i] * xx[i];
    }
    /*
      rxyz2 += yy .* yy
    */
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      recip_rxyz2[i] += (yy_center[i] * yy_center[i]);
    }
    /*
      rxyz2 += zz .* zz
    */
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      recip_rxyz2[i] += (zz_center[i] * zz_center[i]);
    }
    /*
      recip_rxyz2 = 1.0/recip_rxyz2.
    */
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      recip_rxyz2[i] = 1.0/recip_rxyz2[i];
    }
    /*
      track recip_dmax2.
    recip_dmax2 = mobcal_amax(recip_rxyz2,recip_dmax2);
    */
    for (i=0;i< _IVEC_LEN_ ;i++) {
      if (recip_rxyz2[i] > recip_dmax2) {
	recip_dmax2 = recip_rxyz2[i];
      }
    }
    /*
      Store recip_rxyz2 for dmax computation in mobcal_dljpot.
    */
    /*
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      //recip_rxyz2_01[i] = recip_rxyz2[i];
      recip_rxyz2_00[i] = recip_rxyz2[i];
    }
    recip_rxyz2_00 += _IVEC_LEN_; // Caution address arithmetic.
    */
    /*
      recip_rxyz = sqrt(recip_rxyz2).
    */
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      recip_rxyz[i] = sqrt(recip_rxyz2[i]);
    }
    /* 
      Load qpolv_s0.
    */
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      rv[i] = qpolv_s0[i];
    }
    /*
      qpolv_s0 += recip_rxyz .* pc_scale_pcharge.
    */
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      rv[i] += (recip_rxyz[i] * lpc_scale_pcharge[i]);
    }
    /*
      Store qpolv_s0
    */
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      qpolv_s0[i] = rv[i];
    }
    /*
      recip_rxyz3 = recip_rxyz2 .* recip_rxyz;
      Here we take advantage of &recip_rxyz3 = recip_rxyz.
    */
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      recip_rxyz3[i] = recip_rxyz3[i] * recip_rxyz2[i];
    }
    /*
      vptemp2 <- recip_rxyz3 .* pc_scale_pcharge.
    */
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      vptemp2[i] = recip_rxyz3[i] * lpc_scale_pcharge[i];
    }

    /*
      Load dqpolx_s0
    */
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      rv[i] = dqpolx_s0[i];
    }
    /*
      dqpolx_s0 -= vptemp2 .* xx
    */
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      rv[i] -= (vptemp2[i] * xx[i]);
    }
    /*
      Store dqpolx_s0
    */
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      dqpolx_s0[i] = rv[i];
    }


    /*
      Load dqpoly_s0
    */
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      rv[i] = dqpoly_s0[i];
    }
    /*
      dqpoly_s0 -= vptemp2 .* yy
    */
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      rv[i] -= (vptemp2[i] * yy_center[i]);
    }
    /*
      Store dqpoly_s0
    */
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      dqpoly_s0[i] = rv[i];
    }

    /*
      Load dqpolz_s0
    */
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      rv[i] = dqpolz_s0[i];
    }
    /*
      dqpolz_s0 -= vptemp2 .* zz
    */
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      rv[i] -= (vptemp2[i] * zz_center[i]);
    }
    /*
      Store dqpolz_s0
    */
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      dqpolz_s0[i] = rv[i];
    }
    /*
      recip_rxyz6 <- recip_rxyz3 .* recip_rxyz3
      Here we use the fact that &recip_rxyz6 = &recip_rxyz3
    */
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      recip_rxyz6[i] = recip_rxyz6[i] * recip_rxyz6[i];
    }
    /*
      recip_rxyz12 = recip_rxyz6 .* recip_rxyz6
    */
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      recip_rxyz12[i] = recip_rxyz6[i] * recip_rxyz6[i];
    }
    /*
      recip_rxyz12 <- recip_rxyz12 .* ro12lj
    */
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      recip_rxyz12[i] = recip_rxyz12[i] * lro12lj[i];
    }
    /*
      recip_rxyz12 -= recip_rxyz6 .* ro6lj
    */
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      recip_rxyz12[i] -= (recip_rxyz6[i] * lro6lj[i]);
    }
    /*
      Load e00v_s0.
    */
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      rv[i] = e00v_s0[i];
    }
    /*
      e00v_s0 += eox4 .* (recip_rxyz12 .* ro12lj - recip_rxyz .* ro6.j)
      where the expresion in parenthesis has now ben stored in recip_rxyz12.
    */
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      rv[i] += (leox4[i] * recip_rxyz12[i]);
    }
    /*
      Store e00v_s0
    */
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      e00v_s0[i] = rv[i];
    }
    /*
      recip_rxyz8 = recip_rxyz6 .* recip_rxyz2.
      Here we make use of &recip_rxyz8 = &recip_rxyz2
    */
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      recip_rxyz8[i] = recip_rxyz8[i] * recip_rxyz6[i];
    }
    /*
      recip_rxyz14 <- recip_rxyz6 .* recip_rxyz8.
      Here we make use of &recip_rxyz14 = &recip_rxyz6
    */
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      recip_rxyz14[i] = recip_rxyz14[i] * recip_rxyz8[i];
    }
    /*
      de00g <- dro6 .* recip_rxyz8
      Now we make use of  &de00g = &recip_rxyz8
    */
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      de00g[i] = de00g[i] * ldro6[i];
    }
    /*
      de00g -= dro12 .* recip_rxyz14
    */
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      de00g[i] -= (recip_rxyz14[i] * ldro12[i]);
    }
    /*
      de00g = eox4 .* (dro6 .* recip_rxyz8 - dro12 .* recip_rxyz14).
    */
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      de00g[i] = de00g[i] * leox4[i];
    }

    /*
      Load de00x_s0
    */
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      rv[i] = de00x_s0[i];
    }
    /*
      de00x_s0 += de00g .* xx
    */
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      rv[i] += (de00g[i] * xx[i]);
    }
    /*
      Store de00x_s0
    */
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      de00x_s0[i] = rv[i];
    }

    /*
      Load de00y_s0
    */
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      rv[i] = de00y_s0[i];
    }
    /*
      de00y_s0 += de00g .* yy
    */
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      rv[i] += (de00g[i] * yy_center[i]);
    }
    /*
      Store de00y_s0
    */
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      de00y_s0[i] = rv[i];
    }

    /*
      Load de00z_s0
    */
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      rv[i] = de00z_s0[i];
    }
    /*
      de00z_s0 += de00g .* zz
    */
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      rv[i] += (de00g[i] * zz_center[i]);
    }
    /*
      Store de00z_s0
    */
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      de00z_s0[i] = rv[i];
    }
    /*
      isamp = 1, ibatom = 0. 
      xx = xx_center, yy = yy_center -bndc , zz = zz_center
    */
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      yy[i] = yy_center[i] - lbndc[i];
    }
    /*
      rxyz2 <- xx .* xx
    */
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      recip_rxyz2[i] = xx_center[i] * xx_center[i];
    }
    /*
      rxyz2 += yy .* yy
    */
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      recip_rxyz2[i] += (yy[i] * yy[i]);
    }
    /*
      rxyz2 += zz .* zz
    */
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      recip_rxyz2[i] += (zz_center[i] * zz_center[i]);
    }
    /*
      recip_rxyz2 = 1.0/recip_rxyz2.
    */
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      recip_rxyz2[i] = 1.0/recip_rxyz2[i];
    }
    /*
      track recip_dmax2.
    recip_dmax2 = mobcal_amax(recip_rxyz2,recip_dmax2);
    */
    for (i=0;i< _IVEC_LEN_ ;i++) {
      if (recip_rxyz2[i] > recip_dmax2) {
	recip_dmax2 = recip_rxyz2[i];
      }
    }
    /*
      Store recip_rxyz2 for dmax computation in mobcal_dljpot.
    */
    /*
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      //recip_rxyz2_10[i] = recip_rxyz2[i];
      recip_rxyz2_00[i] = recip_rxyz2[i];
    }
    recip_rxyz2_00 += _IVEC_LEN_; // Caution address arithmetic.
    */
    /*
      recip_rxyz = sqrt(recip_rxyz2).
    */
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      recip_rxyz[i] = sqrt(recip_rxyz2[i]);
    }
    /* 
      Load qpolv_s1.
    */
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      rv[i] = qpolv_s1[i];
    }
    /*
      qpolv_s1 += recip_rxyz .* pc_scale_pcharge.
    */
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      rv[i] += (recip_rxyz[i] * lpc_scale_pcharge[i]);
    }
    /*
      Store qpolv_s1
    */
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      qpolv_s1[i] = rv[i];
    }
    /*
      recip_rxyz3 = recip_rxyz2 .* recip_rxyz;
      Here we take advantage of &recip_rxyz3 = recip_rxyz.
    */
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      recip_rxyz3[i] = recip_rxyz3[i] * recip_rxyz2[i];
    }
    /*
      vptemp2 <- recip_rxyz3 .* pc_scale_pcharge.
    */
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      vptemp2[i] = recip_rxyz3[i] * lpc_scale_pcharge[i];
    }

    /*
      Load dqpolx_s1
    */
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      rv[i] = dqpolx_s1[i];
    }
    /*
      dqpolx_s1 -= vptemp2 .* xx
    */
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      rv[i] -= (vptemp2[i] * xx_center[i]);
    }
    /*
      Store dqpolx_s1
    */
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      dqpolx_s1[i] = rv[i];
    }

    /*
      Load dqpoly_s1
    */
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      rv[i] = dqpoly_s1[i];
    }
    /*
      dqpoly_s1 -= vptemp2 .* yy
    */
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      rv[i] -= (vptemp2[i] * yy[i]);
    }
    /*
      Store dqpoly_s1
    */
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      dqpoly_s1[i] = rv[i];
    }

    /*
      Load dqpolz_s1
    */
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      rv[i] = dqpolz_s1[i];
    }
    /*
      dqpolz_s1 -= vptemp2 .* zz
    */
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      rv[i] -= (vptemp2[i] * zz_center[i]);
    }
    /*
      Store dqpolz_s1
    */
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      dqpolz_s1[i] = rv[i];
    }
    /*
      recip_rxyz6 <- recip_rxyz3 .* recip_rxyz3
      Here we use the fact that &recip_rxyz6 = &recip_rxyz3
    */
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      recip_rxyz6[i] = recip_rxyz6[i] * recip_rxyz6[i];
    }
    /*
      recip_rxyz12 = recip_rxyz6 .* recip_rxyz6
    */
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      recip_rxyz12[i] = recip_rxyz6[i] * recip_rxyz6[i];
    }
    /*
      recip_rxyz12 <- recip_rxyz12 .* ro12lj
    */
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      recip_rxyz12[i] = recip_rxyz12[i] * lro12lj[i];
    }
    /*
      recip_rxyz12 -= recip_rxyz6 .* ro6lj
    */
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      recip_rxyz12[i] -= (recip_rxyz6[i] * lro6lj[i]);
    }
    /*
      Load e00v_s1.
    */
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      rv[i] = e00v_s1[i];
    }
    /*
      e00v_s1 += eox4 .* (recip_rxyz12 .* ro12lj - recip_rxyz .* ro6.j)
      where the expresion in parenthesis has now ben stored in recip_rxyz12.
    */
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      rv[i] += (leox4[i] * recip_rxyz12[i]);
    }
    /*
      Store e00v_s1
    */
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      e00v_s1[i] = rv[i];
    }
    /*
      recip_rxyz8 = recip_rxyz6 .* recip_rxyz2.
      Here we make use of &recip_rxyz8 = &recip_rxyz2
    */
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      recip_rxyz8[i] = recip_rxyz8[i] * recip_rxyz6[i];
    }
    /*
      recip_rxyz14 <- recip_rxyz6 .* recip_rxyz8.
      Here we make use of &recip_rxyz14 = &recip_rxyz6
    */
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      recip_rxyz14[i] = recip_rxyz14[i] * recip_rxyz8[i];
    }
    /*
      de00g <- dro6 .* recip_rxyz8
      Now we make use of  &de00g = &recip_rxyz8
    */
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      de00g[i] = de00g[i] * ldro6[i];
    }
    /*
      de00g -= dro12 .* recip_rxyz14
    */
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      de00g[i] -= (recip_rxyz14[i] * ldro12[i]);
    }
    /*
      de00g = eox4 .* (dro6 .* recip_rxyz8 - dro12 .* recip_rxyz14).
    */
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      de00g[i] = de00g[i] * leox4[i];
    }

    /*
      Load de00x_s1
    */
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      rv[i] = de00x_s1[i];
    }
    /*
      de00x_s1 += de00g .* xx
    */
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      rv[i] += (de00g[i] * xx_center[i]);
    }
    /*
      Store de00x_s1
    */
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      de00x_s1[i] = rv[i];
    }

    /*
      Load de00y_s1
    */
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      rv[i] = de00y_s1[i];
    }
    /*
      de00y_s1 += de00g .* yy
    */
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      rv[i] += (de00g[i] * yy[i]);
    }
    /*
      Store de00y_s1
    */
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      de00y_s1[i] = rv[i];
    }

    /*
      Load de00z_s1
    */
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      rv[i] = de00z_s1[i];
    }
    /*
      de00z_s1 += de00g .* zz
    */
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      rv[i] += (de00g[i] * zz_center[i]);
    }
    /*
      Store de00z_s1
    */
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      de00z_s1[i] = rv[i];
    }

    /*
      isamp = 1, ibatom = 1. 
      xx = xx_center, yy = yy_center +bndc , zz = zz_center
    */
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      yy[i] = yy_center[i] + lbndc[i];
    }
    /*
      rxyz2 <- xx .* xx
    */
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      recip_rxyz2[i] = xx_center[i] * xx_center[i];
    }
    /*
      rxyz2 += yy .* yy
    */
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      recip_rxyz2[i] += (yy[i] * yy[i]);
    }
    /*
      rxyz2 += zz .* zz
    */
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      recip_rxyz2[i] += (zz_center[i] * zz_center[i]);
    }
    /*
      recip_rxyz2 = 1.0/recip_rxyz2.
    */
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      recip_rxyz2[i] = 1.0/recip_rxyz2[i];
    }
    /*
      track recip_dmax2.
    recip_dmax2 = mobcal_amax(recip_rxyz2,recip_dmax2);
    */
    for (i=0;i< _IVEC_LEN_ ;i++) {
      if (recip_rxyz2[i] > recip_dmax2) {
	recip_dmax2 = recip_rxyz2[i];
      }
    }
    /*
      Store recip_rxyz2 for dmax computation in mobcal_dljpot.
    */
    /*
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      //recip_rxyz2_11[i] = recip_rxyz2[i];
      recip_rxyz2_00[i] = recip_rxyz2[i];
    }
    recip_rxyz2_00 += _IVEC_LEN_; // Caution address arithmetic.
    */
    /*
      recip_rxyz = sqrt(recip_rxyz2).
    */
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      recip_rxyz[i] = sqrt(recip_rxyz2[i]);
    }
    /* 
      Load qpolv_s1.
    */
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      rv[i] = qpolv_s1[i];
    }
    /*
      qpolv_s1 += recip_rxyz .* pc_scale_pcharge.
    */
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      rv[i] += (recip_rxyz[i] * lpc_scale_pcharge[i]);
    }
    /*
      Store qpolv_s1
    */
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      qpolv_s1[i] = rv[i];
    }
    /*
      recip_rxyz3 = recip_rxyz2 .* recip_rxyz;
      Here we take advantage of &recip_rxyz3 = recip_rxyz.
    */
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      recip_rxyz3[i] = recip_rxyz3[i] * recip_rxyz2[i];
    }
    /*
      vptemp2 <- recip_rxyz3 .* pc_scale_pcharge.
    */
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      vptemp2[i] = recip_rxyz3[i] * lpc_scale_pcharge[i];
    }

    /*
      Load dqpolx_s1
    */
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      rv[i] = dqpolx_s1[i];
    }
    /*
      dqpolx_s1 -= vptemp2 .* xx
    */
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      rv[i] -= (vptemp2[i] * xx_center[i]);
    }
    /*
      Store dqpolx_s1
    */
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      dqpolx_s1[i] = rv[i];
    }

    /*
      Load dqpoly_s1
    */
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      rv[i] = dqpoly_s1[i];
    }
    /*
      dqpoly_s1 -= vptemp2 .* yy
    */
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      rv[i] -= (vptemp2[i] * yy[i]);
    }
    /*
      Store dqpoly_s1
    */
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      dqpoly_s1[i] = rv[i];
    }
    
    /*
      Load dqpolz_s1
    */
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      rv[i] = dqpolz_s1[i];
    }
    /*
      dqpolz_s1 -= vptemp2 .* zz
    */
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      rv[i] -= (vptemp2[i] * zz_center[i]);
    }
    /*
      Store dqpolz_s1
    */
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      dqpolz_s1[i] = rv[i];
    }
    /*
      recip_rxyz6 <- recip_rxyz3 .* recip_rxyz3
      Here we use the fact that &recip_rxyz6 = &recip_rxyz3
    */
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      recip_rxyz6[i] = recip_rxyz6[i] * recip_rxyz6[i];
    }
    /*
      recip_rxyz12 = recip_rxyz6 .* recip_rxyz6
    */
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      recip_rxyz12[i] = recip_rxyz6[i] * recip_rxyz6[i];
    }
    /*
      recip_rxyz12 <- recip_rxyz12 .* ro12lj
    */
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      recip_rxyz12[i] = recip_rxyz12[i] * lro12lj[i];
    }
    /*
      recip_rxyz12 -= recip_rxyz6 .* ro6lj
    */
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      recip_rxyz12[i] -= (recip_rxyz6[i] * lro6lj[i]);
    }
    /*
      Load e00v_s1.
    */
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      rv[i] = e00v_s1[i];
    }
    /*
      e00v_s1 += eox4 .* (recip_rxyz12 .* ro12lj - recip_rxyz .* ro6.j)
      where the expresion in parenthesis has now ben stored in recip_rxyz12.
    */
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      rv[i] += (leox4[i] * recip_rxyz12[i]);
    }
    /*
      Store e00v_s1
    */
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      e00v_s1[i] = rv[i];
    }
    /*
      recip_rxyz8 = recip_rxyz6 .* recip_rxyz2.
      Here we make use of &recip_rxyz8 = &recip_rxyz2
    */
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      recip_rxyz8[i] = recip_rxyz8[i] * recip_rxyz6[i];
    }
    /*
      recip_rxyz14 <- recip_rxyz6 .* recip_rxyz8.
      Here we make use of &recip_rxyz14 = &recip_rxyz6
    */
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      recip_rxyz14[i] = recip_rxyz14[i] * recip_rxyz8[i];
    }
    /*
      de00g <- dro6 .* recip_rxyz8
      Now we make use of  &de00g = &recip_rxyz8
    */
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      de00g[i] = de00g[i] * ldro6[i];
    }
    /*
      de00g -= dro12 .* recip_rxyz14
    */
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      de00g[i] -= (recip_rxyz14[i] * ldro12[i]);
    }
    /*
      de00g = eox4 .* (dro6 .* recip_rxyz8 - dro12 .* recip_rxyz14).
    */
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      de00g[i] = de00g[i] * leox4[i];
    }

    /*
      Load de00x_s1
    */
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      rv[i] = de00x_s1[i];
    }
    /*
      de00x_s1 += de00g .* xx
    */
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      rv[i] += (de00g[i] * xx_center[i]);
    }
    /*
      Store de00x_s1
    */
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      de00x_s1[i] = rv[i];
    }

    /*
      Load de00y_s1
    */
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      rv[i] = de00y_s1[i];
    }
    /*
      de00y_s1 += de00g .* yy
    */
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      rv[i] += (de00g[i] * yy[i]);
    }
    /*
      Store de00y_s1
    */
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      de00y_s1[i] = rv[i];
    }

    /*
      Load de00z_s1
    */
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      rv[i] = de00z_s1[i];
    }
    /*
      de00z_s1 += de00g .* zz
    */
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      rv[i] += (de00g[i] * zz_center[i]);
    }
    /*
      Store de00z_s1
    */
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      de00z_s1[i] = rv[i];
    }

    /*
      isamp = 2, ibatom = 0. 
      xx = xx_center, yy = yy_center , zz = zz_center - bndc
    */
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      zz[i] = zz_center[i] - lbndc[i];
    }
    /*
      rxyz2 <- xx .* xx
    */
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      recip_rxyz2[i] = xx_center[i] * xx_center[i];
    }
    /*
      rxyz2 += yy .* yy
    */
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      recip_rxyz2[i] += (yy_center[i] * yy_center[i]);
    }
    /*
      rxyz2 += zz .* zz
    */
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      recip_rxyz2[i] += (zz[i] * zz[i]);
    }
    /*
      recip_rxyz2 = 1.0/recip_rxyz2.
    */
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      recip_rxyz2[i] = 1.0/recip_rxyz2[i];
    }
    /*
      track recip_dmax2.
    recip_dmax2 = mobcal_amax(recip_rxyz2,recip_dmax2);
    */
    for (i=0;i< _IVEC_LEN_ ;i++) {
      if (recip_rxyz2[i] > recip_dmax2) {
	recip_dmax2 = recip_rxyz2[i];
      }
    }
    /*
      Store recip_rxyz2 for dmax computation in mobcal_dljpot.
    */
    /*
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      //recip_rxyz2_20[i] = recip_rxyz2[i];
      recip_rxyz2_00[i] = recip_rxyz2[i];
    }
    recip_rxyz2_00 += _IVEC_LEN_; // Caution address arithmetic.
    */
    /*
      recip_rxyz = sqrt(recip_rxyz2).
    */
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      recip_rxyz[i] = sqrt(recip_rxyz2[i]);
    }
    /* 
      Load qpolv_s2.
    */
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      rv[i] = qpolv_s2[i];
    }
    /*
      qpolv_s2 += recip_rxyz .* pc_scale_pcharge.
    */
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      rv[i] += (recip_rxyz[i] * lpc_scale_pcharge[i]);
    }
    /*
      Store qpolv_s2
    */
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      qpolv_s2[i] = rv[i];
    }
    /*
      recip_rxyz3 = recip_rxyz2 .* recip_rxyz;
      Here we take advantage of &recip_rxyz3 = recip_rxyz.
    */
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      recip_rxyz3[i] = recip_rxyz3[i] * recip_rxyz2[i];
    }
    /*
      vptemp2 <- recip_rxyz3 .* pc_scale_pcharge.
    */
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      vptemp2[i] = recip_rxyz3[i] * lpc_scale_pcharge[i];
    }

    /*
      Load dqpolx_s2
    */
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      rv[i] = dqpolx_s2[i];
    }
    /*
      dqpolx_s1 -= vptemp2 .* xx
    */
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      rv[i] -= (vptemp2[i] * xx_center[i]);
    }
    /*
      Store dqpolx_s2
    */
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      dqpolx_s2[i] = rv[i];
    }

    /*
      Load dqpoly_s2
    */
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      rv[i] = dqpoly_s2[i];
    }
    /*
      dqpoly_s1 -= vptemp2 .* yy
    */
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      rv[i] -= (vptemp2[i] * yy_center[i]);
    }
    /*
      Store dqpoly_s2
    */
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      dqpoly_s2[i] = rv[i];
    }
    
    /*
      Load dqpolz_s2
    */
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      rv[i] = dqpolz_s2[i];
    }
    /*
      dqpolz_s1 -= vptemp2 .* zz
    */
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      rv[i] -= (vptemp2[i] * zz[i]);
    }
    /*
      Store dqpolz_s2
    */
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      dqpolz_s2[i] = rv[i];
    }
    /*
      recip_rxyz6 <- recip_rxyz3 .* recip_rxyz3
      Here we use the fact that &recip_rxyz6 = &recip_rxyz3
    */
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      recip_rxyz6[i] = recip_rxyz6[i] * recip_rxyz6[i];
    }
    /*
      recip_rxyz12 = recip_rxyz6 .* recip_rxyz6
    */
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      recip_rxyz12[i] = recip_rxyz6[i] * recip_rxyz6[i];
    }
    /*
      recip_rxyz12 <- recip_rxyz12 .* ro12lj
    */
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      recip_rxyz12[i] = recip_rxyz12[i] * lro12lj[i];
    }
    /*
      recip_rxyz12 -= recip_rxyz6 .* ro6lj
    */
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      recip_rxyz12[i] -= (recip_rxyz6[i] * lro6lj[i]);
    }
    /*
      Load e00v_s2.
    */
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      rv[i] = e00v_s2[i];
    }
    /*
      e00v_s1 += eox4 .* (recip_rxyz12 .* ro12lj - recip_rxyz .* ro6.j)
      where the expresion in parenthesis has now ben stored in recip_rxyz12.
    */
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      rv[i] += (leox4[i] * recip_rxyz12[i]);
    }
    /*
      Store e00v_s2
    */
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      e00v_s2[i] = rv[i];
    }
    /*
      recip_rxyz8 = recip_rxyz6 .* recip_rxyz2.
      Here we make use of &recip_rxyz8 = &recip_rxyz2
    */
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      recip_rxyz8[i] = recip_rxyz8[i] * recip_rxyz6[i];
    }
    /*
      recip_rxyz14 <- recip_rxyz6 .* recip_rxyz8.
      Here we make use of &recip_rxyz14 = &recip_rxyz6
    */
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      recip_rxyz14[i] = recip_rxyz14[i] * recip_rxyz8[i];
    }
    /*
      de00g <- dro6 .* recip_rxyz8
      Now we make use of  &de00g = &recip_rxyz8
    */
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      de00g[i] = de00g[i] * ldro6[i];
    }
    /*
      de00g -= dro12 .* recip_rxyz14
    */
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      de00g[i] -= (recip_rxyz14[i] * ldro12[i]);
    }
    /*
      de00g = eox4 .* (dro6 .* recip_rxyz8 - dro12 .* recip_rxyz14).
    */
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      de00g[i] = de00g[i] * leox4[i];
    }

    /*
      Load de00x_s2
    */
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      rv[i] = de00x_s2[i];
    }
    /*
      de00x_s2 += de00g .* xx
    */
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      rv[i] += (de00g[i] * xx_center[i]);
    }
    /*
      Store de00x_s2
    */
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      de00x_s2[i] = rv[i];
    }

    /*
      Load de00y_s2
    */
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      rv[i] = de00y_s2[i];
    }
    /*
      de00y_s2 += de00g .* yy
    */
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      rv[i] += (de00g[i] * yy_center[i]);
    }
    /*
      Store de00y_s2
    */
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      de00y_s2[i] = rv[i];
    }

    /*
      Load de00z_s2
    */
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      rv[i] = de00z_s2[i];
    }
    /*
      de00z_s2 += de00g .* zz
    */
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      rv[i] += (de00g[i] * zz[i]);
    }
    /*
      Store de00z_s2
    */
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      de00z_s2[i] = rv[i];
    }

    /*
      isamp = 2, ibatom = 1. 
      xx = xx_center, yy = yy_center , zz = zz_center + bndc
    */
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      zz[i] = zz_center[i] + lbndc[i];
    }
    /*
      rxyz2 <- xx .* xx
    */
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      recip_rxyz2[i] = xx_center[i] * xx_center[i];
    }
    /*
      rxyz2 += yy .* yy
    */
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      recip_rxyz2[i] += (yy_center[i] * yy_center[i]);
    }
    /*
      rxyz2 += zz .* zz
    */
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      recip_rxyz2[i] += (zz[i] * zz[i]);
    }
    /*
      recip_rxyz2 = 1.0/recip_rxyz2.
    */
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      recip_rxyz2[i] = 1.0/recip_rxyz2[i];
    }
    /*
      track recip_dmax2.
    recip_dmax2 = mobcal_amax(recip_rxyz2,recip_dmax2);
    */
    for (i=0;i< _IVEC_LEN_ ;i++) {
      if (recip_rxyz2[i] > recip_dmax2) {
	recip_dmax2 = recip_rxyz2[i];
      }
    }
    /*
      Store recip_rxyz2 for dmax computation in mobcal_dljpot.
    */
    /*
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      //recip_rxyz2_21[i] = recip_rxyz2[i];
      recip_rxyz2_00[i] = recip_rxyz2[i];
    }
    recip_rxyz2_00 += _IVEC_LEN_; // Caution address arithmetic.
    */
    /*
      recip_rxyz = sqrt(recip_rxyz2).
    */
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      recip_rxyz[i] = sqrt(recip_rxyz2[i]);
    }
    /* 
      Load qpolv_s2.
    */
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      rv[i] = qpolv_s2[i];
    }
    /*
      qpolv_s2 += recip_rxyz .* pc_scale_pcharge.
    */
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      rv[i] += (recip_rxyz[i] * lpc_scale_pcharge[i]);
    }
    /*
      Store qpolv_s2
    */
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      qpolv_s2[i] = rv[i];
    }
    /*
      recip_rxyz3 = recip_rxyz2 .* recip_rxyz;
      Here we take advantage of &recip_rxyz3 = recip_rxyz.
    */
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      recip_rxyz3[i] = recip_rxyz3[i] * recip_rxyz2[i];
    }
    /*
      vptemp2 <- recip_rxyz3 .* pc_scale_pcharge.
    */
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      vptemp2[i] = recip_rxyz3[i] * lpc_scale_pcharge[i];
    }

    /*
      Load dqpolx_s2
    */
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      rv[i] = dqpolx_s2[i];
    }
    /*
      dqpolx_s1 -= vptemp2 .* xx
    */
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      rv[i] -= (vptemp2[i] * xx_center[i]);
    }
    /*
      Store dqpolx_s2
    */
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      dqpolx_s2[i] = rv[i];
    }

    /*
      Load dqpoly_s2
    */
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      rv[i] = dqpoly_s2[i];
    }
    /*
      dqpoly_s1 -= vptemp2 .* yy
    */
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      rv[i] -= (vptemp2[i] * yy_center[i]);
    }
    /*
      Store dqpoly_s2
    */
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      dqpoly_s2[i] = rv[i];
    }
    
    /*
      Load dqpolz_s2
    */
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      rv[i] = dqpolz_s2[i];
    }
    /*
      dqpolz_s1 -= vptemp2 .* zz
    */
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      rv[i] -= (vptemp2[i] * zz[i]);
    }
    /*
      Store dqpolz_s2
    */
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      dqpolz_s2[i] = rv[i];
    }
    /*
      recip_rxyz6 <- recip_rxyz3 .* recip_rxyz3
      Here we use the fact that &recip_rxyz6 = &recip_rxyz3
    */
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      recip_rxyz6[i] = recip_rxyz6[i] * recip_rxyz6[i];
    }
    /*
      recip_rxyz12 = recip_rxyz6 .* recip_rxyz6
    */
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      recip_rxyz12[i] = recip_rxyz6[i] * recip_rxyz6[i];
    }
    /*
      recip_rxyz12 <- recip_rxyz12 .* ro12lj
    */
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      recip_rxyz12[i] = recip_rxyz12[i] * lro12lj[i];
    }
    /*
      recip_rxyz12 -= recip_rxyz6 .* ro6lj
    */
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      recip_rxyz12[i] -= (recip_rxyz6[i] * lro6lj[i]);
    }
    /*
      Load e00v_s2.
    */
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      rv[i] = e00v_s2[i];
    }
    /*
      e00v_s1 += eox4 .* (recip_rxyz12 .* ro12lj - recip_rxyz .* ro6.j)
      where the expresion in parenthesis has now ben stored in recip_rxyz12.
    */
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      rv[i] += (leox4[i] * recip_rxyz12[i]);
    }
    /*
      Store e00v_s2
    */
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      e00v_s2[i] = rv[i];
    }
    /*
      recip_rxyz8 = recip_rxyz6 .* recip_rxyz2.
      Here we make use of &recip_rxyz8 = &recip_rxyz2
    */
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      recip_rxyz8[i] = recip_rxyz8[i] * recip_rxyz6[i];
    }
    /*
      recip_rxyz14 <- recip_rxyz6 .* recip_rxyz8.
      Here we make use of &recip_rxyz14 = &recip_rxyz6
    */
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      recip_rxyz14[i] = recip_rxyz14[i] * recip_rxyz8[i];
    }
    /*
      de00g <- dro6 .* recip_rxyz8
      Now we make use of  &de00g = &recip_rxyz8
    */
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      de00g[i] = de00g[i] * ldro6[i];
    }
    /*
      de00g -= dro12 .* recip_rxyz14
    */
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      de00g[i] -= (recip_rxyz14[i] * ldro12[i]);
    }
    /*
      de00g = eox4 .* (dro6 .* recip_rxyz8 - dro12 .* recip_rxyz14).
    */
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      de00g[i] = de00g[i] * leox4[i];
    }

    /*
      Load de00x_s2
    */
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      rv[i] = de00x_s2[i];
    }
    /*
      de00x_s2 += de00g .* xx
    */
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      rv[i] += (de00g[i] * xx_center[i]);
    }
    /*
      Store de00x_s2
    */
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      de00x_s2[i] = rv[i];
    }

    /*
      Load de00y_s2
    */
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      rv[i] = de00y_s2[i];
    }
    /*
      de00y_s2 += de00g .* yy
    */
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      rv[i] += (de00g[i] * yy_center[i]);
    }
    /*
      Store de00y_s2
    */
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      de00y_s2[i] = rv[i];
    }

    /*
      Load de00z_s2
    */
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      rv[i] = de00z_s2[i];
    }
    /*
      de00z_s2 += de00g .* zz
    */
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      rv[i] += (de00g[i] * zz[i]);
    }
    /*
      Store de00z_s2
    */
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      de00z_s2[i] = rv[i];
    }
    /*
    fx += ivec_len;
    fy += ivec_len;
    fz += ivec_len;
    pcharge += ivec_len;
    eox4    += ivec_len;
    ro6lj   += ivec_len;
    ro12lj  += ivec_len;
    dro6    += ivec_len;
    dro12   += ivec_len;
    mthree_pcharge += ivec_len;
    pc_scale_pcharge += ivec_len;
    */
    fx = dro12 + _IVEC_LEN_;
    fy = fx + _IVEC_LEN_;
    fz = fy + _IVEC_LEN_;
    pc_scale_pcharge = fz + _IVEC_LEN_;
    pcharge          = pc_scale_pcharge + _IVEC_LEN_;
    mthree_pcharge   = pcharge + _IVEC_LEN_;
    eox4             = mthree_pcharge + _IVEC_LEN_;
    ro6lj            = eox4 + _IVEC_LEN_;
    ro12lj           = ro6lj + _IVEC_LEN_;
    dro6             = ro12lj + _IVEC_LEN_;
    dro12            = dro6 + _IVEC_LEN_;
    /*
    recip_rxyz2_00 += ivec_len;
    recip_rxyz2_01 += ivec_len;
    recip_rxyz2_10 += ivec_len;
    recip_rxyz2_11 += ivec_len;
    recip_rxyz2_20 += ivec_len;
    recip_rxyz2_21 += ivec_len;
    */
  }
  *recip_dmax2_p = recip_dmax2;
}
