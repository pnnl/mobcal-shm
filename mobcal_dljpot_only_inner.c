#include "system_includes.h"
#include "mobcal_state_struct.h"
#include "mobcal_dljpot_only_inner.h"
void mobcal_dljpot_only_inner(struct mobcal_state_struct *state) {
  /*
    Workspace needs to be 328 + 6* vec_len  elements long.
    the px, py, pz, and bndc vectors are assumed to be set the calling routine.
    and are int the state->dljpot_workspace vector.
    Called by: mobcal_dljpot_only_dgt
  */
  double *fx;
  double *fy;
  double *fz;
  double *pcharge;
  double *eox4;
  double *ro6lj;
  double *ro12lj;
  double *pc_scale_pcharge;

  /*
    Workspace vectors.
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
  double *lpcharge;
  double *rxyz3i;
  double *lrx;
  double *lry;
  double *lrz;
  double *lro12lj;
  double *lro6lj;
  double *leox4;
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
  double *rv;
  double *workspace;
  int64_t vec_len;

  int avec_len;
  int iatom;

  int ivec_len;
  int i;

  int dmax_offset;
  int dljpot_input_stream_spacing;
  /*
    Input vectors.
  fx      = state->fx;
  fy      = state->fy;
  fz      = state->fz;
  pcharge = state->pcharge;
  eox4    = state->eox4;
  ro6lj   = state->ro6lj;
  ro12lj  = state->ro12lj;
  pc_scale_pcharge = state->pc_scale_pcharge;
  */
  dmax_offset = (int)state->dmax_offset;
  vec_len  = state->vec_len;
  workspace = state->dljpot_workspace;
  avec_len = state->avec_len;
  ivec_len = _IVEC_LEN_;
  dljpot_input_stream_spacing = 11 * _IVEC_LEN_;
  fx = state->dljpot_input_stream;
  /*
    Caution, the following pointer assignments use address arithmetic.
  */
  fy 		   = fx + _IVEC_LEN_;
  fz 		   = fy + _IVEC_LEN_;
  fz 		   = fy + _IVEC_LEN_;
  pc_scale_pcharge = fz + _IVEC_LEN_;
  pcharge          = pc_scale_pcharge + _IVEC_LEN_;
  eox4             = pcharge + _IVEC_LEN_ + _IVEC_LEN_;
  ro6lj            = eox4 + _IVEC_LEN_;
  ro12lj           = ro6lj + _IVEC_LEN_;
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
  /* ymm8 */
  yy_center2 = &xx_center2[avec_len];
  /* ymm9 */
  zz_center2 = &yy_center2[avec_len];
  lbndc       = zz_center2;
  /* ymm10 */
  rxyz_center2 = &zz_center2[avec_len];
  recip_rxyz_center2  = rxyz_center2;
  recip_rxyz12        = rxyz_center2;
  /* ymm11 */
  recip_rxyz_center = &rxyz_center2[avec_len];
  rxyz3i             = recip_rxyz_center;
  xx                 = recip_rxyz_center;
  yy                 = recip_rxyz_center;
  zz                 = recip_rxyz_center;
  /* ymm12 */
  lpc_scale_pcharge = &recip_rxyz_center[avec_len];
  /* ymm13 */
  recip_rxyz_center3 = &lpc_scale_pcharge[avec_len];
  lqpolvi            = recip_rxyz_center3;
  rxyz2              = recip_rxyz_center3;
  recip_rxyz2        = recip_rxyz_center3;
  /* ymm14 */
  lpcharge           = &recip_rxyz_center3[avec_len];
  recip_rxyz         = lpcharge;
  recip_rxyz3        = lpcharge;
  recip_rxyz6        = lpcharge;
  /* ymm15 */
  rv                 = &lpcharge[avec_len];
  lrx                = rv;
  lry                = rv;
  lrz                = rv;

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
    fx += ivec_len;
    fy += ivec_len;
    fz += ivec_len;
    pcharge += ivec_len;
    eox4    += ivec_len;
    ro6lj   += ivec_len;
    ro12lj  += ivec_len;
    pc_scale_pcharge += ivec_len;
    */
    fx += dljpot_input_stream_spacing;
    fy = fx + _IVEC_LEN_;
    fz = fy + _IVEC_LEN_;
    pc_scale_pcharge = fz + _IVEC_LEN_;
    pcharge          = pc_scale_pcharge + _IVEC_LEN_;
    eox4             = pcharge + _IVEC_LEN_ + _IVEC_LEN_;
    ro6lj            = eox4 + _IVEC_LEN_;
    ro12lj           = ro6lj + _IVEC_LEN_;
  }
}
