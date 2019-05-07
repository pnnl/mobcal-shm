#include "system_includes.h"
#include "mobcal_state_struct.h"
#include "mobcal_dljpot_inner_iu.h"
void mobcal_dljpot_inner_iu(struct mobcal_state_struct *state, 
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

  double yregister[128];
  double *xx;
  double *yy;
  double *zz;
  double *xx2;
  double *yy2;
  double *zz2;
  double *rxyz2;
  double *recip_rxyz2;
  double *recip_rxyz;
  double *recip_rxyz3;
  double *xx_yy;
  double *xx_zz;
  double *yy_zz;
  
  double *lpcharge;
  double *rxyz3i;
  double *lmthree_pcharge;
  double *rxyz5i;
  double *lsum1v;
  double *lsum2v;
  double *lsum3v;
  double *lsum4v;
  double *lsum5v;
  double *lsum6v;
  double *lro12lj;
  double *lro6lj;
  double *leox4;
  double *ldro6;
  double *ldro12;
  double *le00v;
  double *recip_rxyz6;
  double *recip_rxyz12;
  double *recip_rxyz8;
  double *recip_rxyz14;
  double *ro6lj_r6;
  double *r12_m_r6;
  double *dr12_r14;
  double *dr6_m_dr12;

  double *de00;
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

  int dljpot_input_stream_spacing;
  int padi;

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

  vec_len  = state->vec_len;
  workspace = state->dljpot_workspace;
  avec_len =  state->avec_len;
  ivec_len = _IVEC_LEN_;
  /*
    local vectors:
  */
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

  /*
    Start of local ymm registers.
  */
  /* ymm0 */
  xx = &yregister[0];

  /* ymm1 */
  yy = &xx[avec_len];
  /* ymm2 */
  zz = &yy[avec_len];
  /* ymm3 */
  xx2 = &zz[avec_len];
  xx_yy = xx2;
  lro6lj = xx2;
  ro6lj_r6 = xx2;
  ldro12  = xx2;
  dr12_r14 = xx2;
  dr6_m_dr12 = xx2;
  de00       = xx2;
  /* ymm4 */
  yy2 = &xx2[avec_len];
  xx_zz = yy2;
  lro12lj = yy2;
  r12_m_r6 = yy2;
  ldro6    = yy2;
  /* ymm5 */
  zz2 = &yy2[avec_len];
  yy_zz = zz2;
  leox4 = zz2;
  /* ymm6 */
  rxyz2 = &zz2[avec_len];
  recip_rxyz2 = rxyz2;
  /* ymm7 */
  recip_rxyz  = &rxyz2[avec_len];
  recip_rxyz3 = recip_rxyz;
  recip_rxyz6 = recip_rxyz3;
  recip_rxyz8 = recip_rxyz6;
  /* ymm8 */
  lpcharge    = &recip_rxyz[avec_len];
  rxyz3i      = lpcharge;
  lmthree_pcharge = rxyz3i;
  rxyz5i      = lmthree_pcharge;
  recip_rxyz12 = lpcharge;
  recip_rxyz14 = lpcharge;
  /* ymm9 */
  lsum1v       = &lpcharge[avec_len];
  /* ymm10 */
  lsum2v       = &lsum1v[avec_len];
  /* ymm11 */
  lsum3v       = &lsum2v[avec_len];
  /* ymm12 */
  lsum4v       = &lsum3v[avec_len];
  /* ymm13 */
  lsum5v       = &lsum4v[avec_len];
  /* ymm14 */
  lsum6v       = &lsum5v[avec_len];
  /* ymm15 */
  rv          = &lsum6v[avec_len];
  le00v       = rv;
  for (i=0;i< _IVEC_LEN_;i++) {
    lsum1v[i] = 0.0;
  }
  for (i=0;i< _IVEC_LEN_;i++) {
    lsum2v[i] = 0.0;
  }
  for (i=0;i< _IVEC_LEN_;i++) {
    lsum3v[i] = 0.0;
  }
  for (i=0;i< _IVEC_LEN_;i++) {
    lsum4v[i] = 0.0;
  }
  for (i=0;i< _IVEC_LEN_;i++) {
    lsum5v[i] = 0.0;
  }
  for (i=0;i< _IVEC_LEN_;i++) {
    lsum6v[i] = 0.0;
  }

  for (iatom=0;iatom < vec_len;iatom += _IVEC_LEN_) {
#pragma simd
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_;i++) {
      xx[i] = px[i] - fx[i];
    }
#pragma simd
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_;i++) {
      yy[i] = py[i] - fy[i];
    }
#pragma simd
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_;i++) {
      zz[i] = pz[i] - fz[i];
    }
#pragma simd
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_;i++) {
      xx2[i] = xx[i] * xx[i];
    }
#pragma simd
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_;i++) {
      yy2[i] = yy[i] * yy[i];
    }
    /*
      rxyz2 = xx2 + yy2.
    */
/*
  y0      y1      y2      y3      y4      y5      y6        y7
  xx      yy      zz      xx2     yy2            rxyz2
*/
#pragma simd
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_;i++) {
      rxyz2[i] = xx2[i] + yy2[i];
    }
    /*
      zz2 = zz .* zz
    */
/*
  y0      y1      y2      y3      y4      y5      y6        y7
  xx      yy      zz      xx2     yy2     zz2    rxyz2
*/
#pragma simd
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_;i++) {
      zz2[i] = zz[i] * zz[i];
    }
    /*
      rxyz2 += zz2
    */
/*
  y0      y1      y2      y3      y4      y5      y6        y7
  xx      yy      zz      xx2     yy2     zz2    rxyz2
*/
#pragma simd
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_;i++) {
      rxyz2[i] = rxyz2[i] + zz2[i];
    }
    /*
      recip_rxyz2 = 1.0 ./ rxyz2;
      Here we take advantage of the face that recip_rxyz2 points rxyz2.
      mobcal_pg1dq(ivec_len,recip_rxyz2,recip_rxyz2);
    */
/*
  y0      y1      y2      y3      y4      y5      y6        y7
  xx      yy      zz      xx2     yy2     zz2    1/rxyz2
*/
#pragma simd
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      recip_rxyz2[i] = 1.0/recip_rxyz2[i];
    }
    /*
      recip_rxyz <- sqrt(recip_rxyz2)
    */
/*
  y0      y1      y2      y3      y4      y5      y6        y7
  xx      yy      zz      xx2     yy2     zz2    1/rxyz2  1/rxyz
*/
#pragma simd
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      recip_rxyz[i] = sqrt(recip_rxyz2[i]);
    }
/*
  y0      y1      y2      y3      y4      y5      y6        y7
  xx      yy      zz      xx2     yy2     zz2    1/rxyz2  1/rxyz3
*/
#pragma simd
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      recip_rxyz3[i] = recip_rxyz[i] * recip_rxyz2[i];
    }
    /*
      Load pcharge.
    */
/*
  y0       y1       y2 	     y3       y4       y5       y6     	   y7
  xx       yy       zz 	     xx2      yy2      zz2      1/rxyz2    1/rxyz3

  y8       y9       y10	     y11      y12      y13      y14    	   y15
  lpcharge
*/
#pragma simd
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      lpcharge[i] = pcharge[i];
    }
    /*
      rxyz3i <- pcharge .* recip_rxyz3.
    */
/*
  y0       y1       y2 	     y3       y4       y5       y6     	   y7
  xx       yy       zz 	     xx2      yy2      zz2      1/rxyz2    1/rxyz3

  y8       y9       y10	     y11      y12      y13      y14    	   y15
  rxyz3i                                                           
*/
#pragma simd
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      rxyz3i[i] = rxyz3i[i] * recip_rxyz3[i];
    }
    /*
      Load rxv.
    */
/*
  y0       y1       y2 	     y3       y4       y5       y6     	   y7
  xx       yy       zz 	     xx2      yy2      zz2      1/rxyz2    1/rxyz3

  y8       y9       y10	     y11      y12      y13      y14    	   y15
  rxyz3i                                                           rxv(rv)
*/
#pragma simd
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      rv[i] = rxv[i];
    }
    /*
      rx += rxyz3i .* xx_center.
    */
#pragma simd
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      rv[i] += rxyz3i[i] * xx[i];
    }
    /*
      Store rxv.
    */
#pragma simd
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      rxv[i] = rv[i];
    }

    /*
      Load ryv.
    */
/*
  y0       y1       y2 	     y3       y4       y5       y6     	   y7
  xx       yy       zz 	     xx2      yy2      zz2      1/rxyz2    1/rxyz3

  y8       y9       y10	     y11      y12      y13      y14    	   y15
  rxyz3i                                                           ryv(rv)
*/
#pragma simd
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      rv[i] = ryv[i];
    }
    /*
      ry += rxyz3i .* yy_center.
    */
#pragma simd
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      rv[i] += rxyz3i[i] * yy[i];
    }
    /*
      Store ryv.
    */
#pragma simd
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      ryv[i] = rv[i];
    }

    /*
      Load rzv.
    */
/*
  y0       y1       y2 	     y3       y4       y5       y6     	   y7
  xx       yy       zz 	     xx2      yy2      zz2      1/rxyz2    1/rxyz3

  y8       y9       y10	     y11      y12      y13      y14    	   y15
  rxyz3i   lsum1v   lsum2v   lsum3v   lsum4v   lsum5v   lsum6v     rzv(rv)
*/
#pragma simd
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      rv[i] = rzv[i];
    }
    /*
      rzv += rxyz3i .* zz_center.
    */
#pragma simd
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      rv[i] += rxyz3i[i] * zz[i];
    }
    /*
      Store rzv.
    */
#pragma simd
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      rzv[i] = rv[i];
    }
    /*
      Add rxyz3i to sum1, sum4 and sum6.
    */
    /*
      lsum1v += rxyz3i;
    */
#pragma simd
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      lsum1v[i] += rxyz3i[i];
    }
    /*
      lsum4 += rxyz3i;
    */
#pragma simd
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      lsum4v[i] += rxyz3i[i];
    }
    /*
      lsum6 += rxyz3i;
    */
#pragma simd
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      lsum6v[i] += rxyz3i[i];
    }
    /*
      Load mthree_pcharge over writing rxyz3i.
    */
/*
  y0        y1       y2      y3       y4       y5       y6     	   y7
  xx        yy       zz      xx2      yy2      zz2      1/rxyz2    1/rxyz3

  y8        y9       y10     y11      y12      y13      y14    	   y15
  -3pcharge lsum1v   lsum2v  lsum3v   lsum4v   lsum5v   lsum6v     rzv(rv)
*/
#pragma simd
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      lmthree_pcharge[i] = mthree_pcharge[i];
    }
    /*
      rxyz5i = mthree_pccharge .* recip_rxyz3 .* recip_rxyz2.
      Here we take advantage of &rxyz5i = &mthree_pcharge.
    */
#pragma simd
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      rxyz5i[i] = rxyz5i[i] * recip_rxyz3[i];
    }
/*
  y0        y1       y2      y3       y4       y5       y6     	   y7
  xx        yy       zz      xx2      yy2      zz2      1/rxyz2    1/rxyz3

  y8        y9       y10     y11      y12      y13      y14    	   y15
  rxyz5i    lsum1v   lsum2v  lsum3v   lsum4v   lsum5v   lsum6v     rzv(rv)
*/
#pragma simd
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      rxyz5i[i] = rxyz5i[i] * recip_rxyz2[i];
    }
    /*
      lsum1v += rxyz5i .* xx2;
    */
#pragma simd
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      lsum1v[i] += (rxyz5i[i] * xx2[i]);
    }
    /*
      lsum4 += rxyz5i .* yy2;
    */
#pragma simd
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      lsum4v[i] += (rxyz5i[i] * yy2[i]);
    }
    /*
      sum6v += rxyz5i * zz2;
    */
#pragma simd
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      lsum6v[i] += (rxyz5i[i] * zz2[i]);
    }
    /*
      Now we are finished with xx2, yy2 and zz2 and may overwrite them
      xx_yy <- xx.* yy
    */
/*
  y0        y1       y2      y3       y4       y5       y6     	   y7
  xx        yy       zz      xx_yy    yy2      zz2      1/rxyz2    1/rxyz3

  y8        y9       y10     y11      y12      y13      y14    	   y15
  rxyz5i    lsum1v   lsum2v  lsum3v   lsum4v   lsum5v   lsum6v      rzv(rv)
*/
#pragma simd
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      xx_yy[i] = xx[i] * yy[i];
    }
    /*
      lsum2v += xx_yy .* rxyz5i
    */
#pragma simd
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      lsum2v[i] += (xx_yy[i] * rxyz5i[i]);
    }
/*
  y0        y1       y2      y3       y4       y5       y6     	   y7
  xx        yy       zz      xx_yy    xx_zz    zz2      1/rxyz2    1/rxyz3

  y8        y9       y10     y11      y12      y13      y14    	   y15
  rxyz5i    lsum1v   lsum2v  lsum3v   lsum4v   lsum5v   lsum6v     rzv(rv)
*/

    /*
      xx_zz <- xx .* zz
    */
#pragma simd
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      xx_zz[i] = xx[i] * zz[i];
    }
    /*
      lsum3 += xx_zz .* rxyz5i
    */
#pragma simd
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      lsum3v[i] += (xx_zz[i] * rxyz5i[i]);
    }
/*
  y0        y1       y2      y3       y4       y5       y6     	   y7
  xx        yy       zz      xx_yy    xx_zz    yy_zz    1/rxyz2    1/rxyz3

  y8        y9       y10     y11      y12      y13      y14    	   y15
  rxyz5i    lsum1v   lsum2v  lsum3v   lsum4v   lsum5v   lsum6v     rzv(rv)
*/
    /*
      yy_zz <- yy .* zz
    */
#pragma simd
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      yy_zz[i] = yy[i] * zz[i];
    }
    /*
      lsum5v += yy_zz .* rxyz5i
    */
#pragma simd
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      lsum5v[i] += (yy_zz[i] * rxyz5i[i]);
    }
    /*
      recip_rxyz6 <- recip_rxyz3 .* recip_rxyz3
      Here we use the fact that &recip_rxyz6 = &recip_rxyz3
    */
/*
  y0        y1       y2      y3       y4       y5       y6     	   y7
  xx        yy       zz      xx_yy    xx_zz    yy_zz    1/rxyz2    1/rxyz6

  y8        y9       y10     y11      y12      y13      y14    	   y15
  rxyz5i    lsum1v   lsum2v  lsum3v   lsum4v   lsum5v   lsum6v     rzv(rv)
*/
#pragma simd
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      recip_rxyz6[i] = recip_rxyz6[i] * recip_rxyz6[i];
    }
    /*
      Load ro6lj
    */
/*
  y0        y1       y2      y3       y4       y5       y6     	   y7
  xx        yy       zz      lro6lj   xx_zz    yy_zz    1/rxyz2    1/rxyz6

  y8        y9       y10     y11      y12      y13      y14    	   y15
  rxyz5i    lsum1v   lsum2v  lsum3v   lsum4v   lsum5v   lsum6v     rzv(rv)
*/
#pragma simd
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      lro6lj[i]  = ro6lj[i];
    }
    /*
      ro6lj_r6 = ro6lj * recip_rxyz6 (using &ro6lj = &ro6lj_r6)
    */
/*
  y0        y1       y2      y3       y4       y5       y6     	   y7
  xx        yy       zz      ro6lj_r6  xx_zz    yy_zz    1/rxyz2    1/rxyz6

  y8        y9       y10     y11      y12      y13      y14    	   y15
  rxyz5i    lsum1v   lsum2v  lsum3v   lsum4v   lsum5v   lsum6v     rzv(rv)
*/
#pragma simd
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      ro6lj_r6[i]  = ro6lj_r6[i] * recip_rxyz6[i];
    }
/*
  y0        y1       y2      y3       y4       y5       y6     	   y7
  xx        yy       zz      ro6lj_r6 xx_zz    yy_zz    1/rxyz2    1/rxyz6

  y8        y9       y10     y11      y12      y13      y14    	   y15
  1/rxyz12  lsum1v   lsum2v  lsum3v   lsum4v   lsum5v   lsum6v     rzv(rv)
*/
#pragma simd
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      recip_rxyz12[i] = recip_rxyz6[i] * recip_rxyz6[i];
    }
    /*
      Load ro12lj
    */
/*
  y0        y1       y2      y3       y4       y5       y6     	   y7
  xx        yy       zz      ro6lj_r6 lro12lj  yy_zz    1/rxyz2    1/rxyz6

  y8        y9       y10     y11      y12      y13      y14    	   y15
  1/rxyz12  lsum1v   lsum2v  lsum3v   lsum4v   lsum5v   lsum6v     rzv(rv)
*/
#pragma simd
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      lro12lj[i] = ro12lj[i];
    }
/*
  form ro12lj*recip_rxyz12 - ro6lj*recip_rxyz6 in r12_m_r6 
  using the fact that &r12_m_r6 = &rol12lj.
  y0        y1       y2      y3       y4       y5       y6     	   y7
  xx        yy       zz      ro6lj_r6 r12_m_r6 yy_zz    1/rxyz2    1/rxyz6

  y8        y9       y10     y11      y12      y13      y14    	   y15
  1/rxyz12  lsum1v   lsum2v  lsum3v   lsum4v   lsum5v   lsum6v     rzv(rv)
*/
#pragma simd
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      r12_m_r6[i] = (r12_m_r6[i] * recip_rxyz12[i]) - ro6lj_r6[i];
    }
    /*
      Load eox4
    */
/*
  y0        y1       y2      y3       y4       y5       y6     	   y7
  xx        yy       zz      ro6lj_r6 r12_m_r6 le0x4    1/rxyz2    1/rxyz6

  y8        y9       y10     y11      y12      y13      y14    	   y15
  1/rxyz12  lsum1v   lsum2v  lsum3v   lsum4v   lsum5v   lsum6v     rzv(rv)
*/
#pragma simd
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      leox4[i]   = eox4[i];
    }
    /*
      Load e00;
    */
/*
  y0        y1       y2      y3       y4       y5       y6     	   y7
  xx        yy       zz      ro6lj_r6 r12_m_r6 le0x4    1/rxyz2    1/rxyz6

  y8        y9       y10     y11      y12      y13      y14    	   y15
  1/rxyz12  lsum1v   lsum2v  lsum3v   lsum4v   lsum5v   lsum6v     le00v
*/
#pragma simd
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ; i++) {
      le00v[i] = e00v[i];
    }
    /*
      Form e00 += eox4 .* r12_m_r6
    */
#pragma simd
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ; i++) {
      le00v[i] = le00v[i] + (leox4[i] * r12_m_r6[i]);
    }
    /*
      Store e00.
    */
#pragma simd
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ; i++) {
      e00v[i] = le00v[i];
    }
    /*
      Compute recip_rxyz14 overwriting recip_rxyz12
    */
/*
  y0        y1       y2      y3       y4       y5       y6     	   y7
  xx        yy       zz      ro6lj_r6 r12_m_r6 le0x4    1/rxyz2    1/rxyz6

  y8        y9       y10     y11      y12      y13      y14    	   y15
  1/rxyz14  lsum1v   lsum2v  lsum3v   lsum4v   lsum5v   lsum6v     le00v
*/
#pragma simd
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ; i++) {
      recip_rxyz14[i] = recip_rxyz14[i] * recip_rxyz2[i];
    }
    /*
      Load dro12
    */
/*
  y0        y1       y2      y3       y4       y5       y6     	   y7
  xx        yy       zz      ldro12   r12_m_r6 le0x4    1/rxyz2    1/rxyz6

  y8        y9       y10     y11      y12      y13      y14    	   y15
  1/rxyz14  lsum1v   lsum2v  lsum3v   lsum4v   lsum5v   lsum6v     le00v
*/
#pragma simd
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      ldro12[i]  = dro12[i];
    }
    /*
      Form ldro12 * recip_rxyz14 storing in dr12_r14 which overwrites ldro12
    */
/*
  y0        y1       y2      y3       y4       y5       y6     	   y7
  xx        yy       zz      dr12_r14 r12_m_r6 le0x4    1/rxyz2    1/rxyz6

  y8        y9       y10     y11      y12      y13      y14    	   y15
  1/rxyz14  lsum1v   lsum2v  lsum3v   lsum4v   lsum5v   lsum6v     le00v
*/
#pragma simd
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      dr12_r14[i]  = dr12_r14[i] * recip_rxyz14[i];
    }
    /*
      Load dro6
    */
/*
  y0        y1       y2      y3       y4       y5       y6     	   y7
  xx        yy       zz      dr12_r14 ldro6    le0x4    1/rxyz2    1/rxyz6

  y8        y9       y10     y11      y12      y13      y14    	   y15
  1/rxyz14  lsum1v   lsum2v  lsum3    lsum4v   lsum5v   lsum6v     le00v
*/
#pragma simd
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ;i++) {
      ldro6[i]   = dro6[i];
    }
    /*
      Compute recip_rxyz8 overwriting recip_rxyz6.
    */
/*
  y0        y1       y2      y3       y4       y5       y6     	   y7
  xx        yy       zz      dr12_r14 ldro6    le0x4    1/rxyz2    1/rxyz8

  y8        y9       y10     y11      y12      y13      y14    	   y15
  1/rxyz14  lsum1v   lsum2v  lsum3v   lsum4v   lsum5v   lsum6v     le00v
*/
#pragma simd
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ; i++) {
      recip_rxyz8[i] = recip_rxyz8[i] * recip_rxyz2[i];
    }
    /*
      Now overwrite dr12_r14 with dr6_m_dr12 = dro6 *recip_rxyz8 - dr12_r14 
    */
/*
  y0        y1       y2      y3         y4       y5       y6         y7
  xx        yy       zz      dr6_m_dr12 ldro6    le0x4    1/rxyz2    1/rxyz8

  y8        y9       y10     y11        y12      y13      y14        y15
  1/rxyz14  lsum1v   lsum2v  lsum3v     lsum4v   lsum5v   lsum6v     le00v
*/
#pragma simd
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ; i++) {
      dr6_m_dr12[i] = (ldro6[i] * recip_rxyz8[i]) - dr6_m_dr12[i];
    }
    /*
      Now compute de00 = leox4 .* dr6_m_dr12
    */
/*
  y0        y1       y2      y3         y4       y5       y6         y7
  xx        yy       zz      de00       ldro6    le0x4    1/rxyz2    1/rxyz8

  y8        y9       y10     y11        y12      y13      y14        y15
  1/rxyz14  lsum1v   lsum2v  lsum3v     lsum4v   lsum5v   lsum6v     le00v
*/
#pragma simd
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ; i++) {
      de00[i] = leox4[i]*dr6_m_dr12[i];
    }
    /*
      Now load,update and sotre the de00xv, de00yv, and de00zv vectors.
    */
#pragma simd
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ; i++) {
      rv[i] = de00xv[i];
    }
    /*
      Update de00xv
    */
#pragma simd
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ; i++) {
      rv[i] = rv[i] + (de00[i]*xx[i]);
    }
    /*
      Store de00xv
    */
#pragma simd
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ; i++) {
      de00xv[i] = rv[i];
    }
    /*
      Load de00yv
    */
#pragma simd
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ; i++) {
      rv[i] = de00yv[i];
    }
    /*
      Update de00yv
    */
#pragma simd
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ; i++) {
      rv[i] = rv[i] + (de00[i]*yy[i]);
    }
    /*
      Store de00yv
    */
#pragma simd
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ; i++) {
      de00yv[i] = rv[i];
    }
    /*
      Load de00zv
    */
#pragma simd
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ; i++) {
      rv[i] = de00zv[i];
    }
    /*
      Update de00zv
    */
#pragma simd
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ; i++) {
      rv[i] = rv[i] + (de00[i]*zz[i]);
    }
    /*
      Store de00zv
    */
#pragma simd
#pragma vector aligned
    for (i=0;i< _IVEC_LEN_ ; i++) {
      de00zv[i] = rv[i];
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
  } /* end for (iatom ... ) */
  /*
    Need to recover the sum1v ... sum6v vectors.
  */
  for (i=0;i< _IVEC_LEN_;i++) {
    sum1v[i] = lsum1v[i];
  }
  for (i=0;i< _IVEC_LEN_;i++) {
    sum2v[i] = lsum2v[i];
  }
  for (i=0;i< _IVEC_LEN_;i++) {
    sum3v[i] = lsum3v[i];
  }
  for (i=0;i< _IVEC_LEN_;i++) {
    sum4v[i] = lsum4v[i];
  }
  for (i=0;i< _IVEC_LEN_;i++) {
    sum5v[i] = lsum5v[i];
  }
  for (i=0;i< _IVEC_LEN_;i++) {
    sum6v[i] = lsum6v[i];
  }
  *recip_dmax2_p = recip_dmax2;
}
