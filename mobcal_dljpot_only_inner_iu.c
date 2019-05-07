#include "system_includes.h"
#include "mobcal_state_struct.h"
#include "mobcal_dljpot_only_inner_iu.h"
void mobcal_dljpot_only_inner_iu(struct mobcal_state_struct *state) {
  /*
    Build the rx,ry,rz, and e00 components of the L-J + ion-dipole potential.
  
    Workspace needs to be 328 + 6* vec_len  elements long.
    the px, py, pz, and bndc vectors are assumed to be set by the 
    calling routine. and are int the state->dljpot_workspace vector.
    Called by: mobcal_dljpot_only_iu
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

  double *rx;
  double *ry;
  double *rz;

  double *e00v;


  double yregister[128];
  double *xx;
  double *yy;
  double *zz;
  double *rxyz2;
  double *recip_rxyz2;
  double *recip_rxyz;
  double *recip_rxyz3;
  double *lpcharge;
  double *rxyz3i;
  double *lrx;
  double *lry;
  double *lrz;
  double *lro12lj;
  double *lro6lj;
  double *leox4;
  double *le00v;
  double *recip_rxyz6;
  double *recip_rxyz12;
  double *workspace;
  int64_t vec_len;

  int avec_len;
  int iatom;

  int ivec_len;
  int i;

  int dljpot_input_stream_spacing;
  int padi;
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
    local input vectors:
  */
  px = workspace;
  py = &px[avec_len];
  pz = &py[avec_len];
  /*
    local output vectors: (assumed to be initialized to 0).
  */
  rx      = &pz[avec_len];
  ry      = &rx[avec_len];
  rz      = &ry[avec_len];

  e00v    = &rz[avec_len];

  /*
  qpolv_s1  = &e00v_s0[avec_len];
  e00v_s1   = &qpolv_s1[avec_len];

  qpolv_s2  = &e00v_s1[avec_len];
  e00v_s2   = &qpolv_s2[avec_len];
  */
  /*
    Start of local ymm registers.
  */
  /* ymm0 */
  lro12lj   = &yregister[0];

  /* ymm1 */
  lro6lj    = &lro12lj[avec_len];
  /* ymm2 */
  leox4     = &lro6lj[avec_len];
  /* ymm3 */		  
  xx = &leox4[avec_len];
  /* ymm4 */
  yy = &xx[avec_len];
  /* ymm5 */
  zz = &yy[avec_len];
  /* ymm6 */
  rxyz2 = &zz[avec_len];
  recip_rxyz2 = rxyz2;
  /* ymm7 */
  recip_rxyz   = &rxyz2[avec_len];
  recip_rxyz3  = recip_rxyz;
  recip_rxyz6  = recip_rxyz3;
  recip_rxyz12 = recip_rxyz6;
  /* ymm8 */
  rxyz3i      = &recip_rxyz[avec_len];
  lpcharge    = rxyz3i;

  /* ymm9 */
  lrx         = &rxyz3i[avec_len];
  /* ymm10 */
  lry         = &lrx[avec_len];
  /* ymm11 */
  lrz         = &lry[avec_len];

  /* ymm12 */
  le00v        = &lrz[avec_len];
  
  for (i=0;i<_IVEC_LEN_;i++) {
    lrx[i]    = rx[i];
    lry[i]    = ry[i];
    lrz[i]    = rz[i];
    le00v[i]  = e00v[i];
  }
  for (iatom=0;iatom < vec_len;iatom += _IVEC_LEN_) {
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_;i++) {
      xx[i] = px[i] - fx[i];
    }
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_;i++) {
      yy[i] = py[i] - fy[i];
    }
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_;i++) {
      zz[i] = pz[i] - fz[i];
    }
    /*
      Now we accumulate the sum of the coordinates squared in recip_rxyz2
    */
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_;i++) {
      rxyz2[i] = xx[i] * xx[i];
    }
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_;i++) {
      rxyz2[i] += (yy[i] * yy[i]);
    }
    /*
      rxyz2 = xx_center2 + yy_center2.
    */
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_;i++) {
      rxyz2[i] += (zz[i]*zz[i]);
    }
    /*
      recip_rxyz2 = 1.0 ./ rxyz2;
      Here we take advantage of the face that recip_rxyz2 points to rxyz2.
      We use recip_rxyz2 on both sized of the equation so as to
      to be obvious to the compiler that they have the same address.
      mobcal_pg1dq(ivec_len,recip_rxyz2,recip_rxyz2);
    */
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      recip_rxyz2[i] = 1.0/recip_rxyz2[i];
    }
    /*
      recip_rxyz <- sqrt(recip_rxyz2)
    */
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      recip_rxyz[i] = sqrt(recip_rxyz2[i]);
    }
    /*
      Now we overwrite recip_rxyz with recip_rxyz3 taking advantage of
      the fact they point to the same vector, but again taking care
      not to confuse the compiler.
    */
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      recip_rxyz3[i] = recip_rxyz3[i] * recip_rxyz2[i];
    }
    /*
      Load pcharge.
    */
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      rxyz3i[i] = pcharge[i];
    }
    /*
      rxyz3i <- pcharge .* recip_rxyz3
      Here we use the fact that rxyz3i and pcharge are the same vector.
    */
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      rxyz3i[i] = rxyz3i[i] * recip_rxyz3[i];
    }
    /*
      rx += rxyz3i .* xx_center.
    */
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      lrx[i] += rxyz3i[i] * xx[i];
    }
    /*
      ry += rxyz3i .* yy_center.
    */
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      lry[i] += rxyz3i[i] * yy[i];
    }
    /*
      rz += rxyz3i .* zz_center.
    */
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      lrz[i] += rxyz3i[i] * zz[i];
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
      Load ro6lj
    */
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      lro6lj[i]  = ro6lj[i];
    }
    /*
      Now put ro6lj * recip_rxyz6 in ro5lj.
    */
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      lro6lj[i]  = lro6lj[i] * recip_rxyz6[i];
    }
    /*
      recip_rxyz12 <- recip_rxyz6 .* recip_rxyz6
      Here we use the fact that &recip_rxyz12 = &recip_rxyz3
    */
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      recip_rxyz12[i] = recip_rxyz12[i] * recip_rxyz12[i];
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
      Now form ro12lj * recip_rxyz12 - ro6lj * recip_rxyz6 in ro12lj.
    */
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      lro12lj[i] = lro12lj[i]*recip_rxyz12[i] - lro6lj[i];
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
      Now form e00 = e00 + (eox4*(ro12lj*recip_rxyz12 - r06lj*recip_rxyz6))
    */
#pragma ivdep
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      le00v[i] = le00v[i] + (leox4[i] * lro12lj[i]);
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
  } /* end for (iatom...) */
  /*
    Now need to store lrx, lry, lrz, and le00 in the workspace.
  */
  for (i=0;i<_IVEC_LEN_;i++) {
    rx[i] = lrx[i];
    ry[i] = lry[i];
    rz[i] = lrz[i];
    e00v[i] = le00v[i];
  }
}
