#include "system_includes.h"
#include "blas.h"
#include "mobcal_state_struct.h"
#include "mobcal_compute_centroid.h"
int mobcal_compute_centroid(struct mobcal_state_struct *state) {
  /*
    Compute the molecule center of mass, and set the
    fxo, fyo and fzo fields of state.
    Called by: mobcal_ncoord
    Calls      ddot, fprintf, fflush
  */
  double *xmass;
  double *fx;
  double *fy;
  double *fz;
  double fxo;
  double fyo;
  double fzo;
  double m2;
  double recip_m2;

  int success;
  int inatom;

  int incx;
  int incy;

  int i;
  int iatom;
  
  int ivec_len;
  int dljpot_input_stream_spacing;


  FILE *ofp;
  FILE *efp;
  success = 1;
  incx    = 1;
  incy    = 1;
  inatom  = state->inatom;
  m2      = state->m2;
  ofp     = state->ofp;
  xmass   = state->xmass;
  /*
  fx      = state->fx;
  fy      = state->fy;
  fz      = state->fz;
  fxo     = ddot_(&inatom,fx,&incx,xmass,&incy);
  fyo     = ddot_(&inatom,fy,&incx,xmass,&incy);
  fzo     = ddot_(&inatom,fz,&incx,xmass,&incy);
  */
  ivec_len  = _IVEC_LEN_;
  dljpot_input_stream_spacing = 11 * _IVEC_LEN_;
  fx = state->dljpot_input_stream;
  fy = fx + _IVEC_LEN_;
  fz = fy + _IVEC_LEN_;
  fxo = 0.0;
  fyo = 0.0;
  fzo = 0.0;
  for (i=0;i<inatom;i+=_IVEC_LEN_) {
    fxo += ddot_(&ivec_len,fx,&incx,xmass,&incy);
    fyo += ddot_(&ivec_len,fy,&incx,xmass,&incy);
    fzo += ddot_(&ivec_len,fz,&incx,xmass,&incy);
    /*
      Caution, the following four pointer assignments use address arithmetic.
    */
    fx += dljpot_input_stream_spacing;
    fy = fx + _IVEC_LEN_;
    fz = fy + _IVEC_LEN_;
    xmass += _IVEC_LEN_;
  }
  if (m2 > 0.0) {
    recip_m2 = 1.0/m2;
  } else {
    success = 0;
    if (state->thread_id == 0) {
      if (ofp) {
	fprintf(ofp,"mobcal_compute_centroid: Error m2 was <= 0.0\n");
	fflush(ofp);
      }
    }
  }
  if (success) {
    fxo = fxo * recip_m2;
    fyo = fyo * recip_m2;
    fzo = fzo * recip_m2;
    state->fxo = fxo;
    state->fyo = fyo;
    state->fzo = fzo;
    if (state->thread_id == 0) {
      if (ofp) {
	fprintf(ofp,"center of mass_coordinates = %le, %le, %le\n",
		fxo,fyo,fzo);
	fflush(ofp);
      }
    }
  }
  return(success);
}
			    
