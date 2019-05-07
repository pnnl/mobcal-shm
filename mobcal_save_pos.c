#include "system_includes.h"
#include "mobcal_state_struct.h"
#include "blas.h"
#include "mobcal_save_pos.h"
int mobcal_save_pos(struct mobcal_state_struct *state) {
  /*
    Copy the coordinate fields from fx,fy,fz,  to ox,oy,oz
    relies on fx,fy,fz being stored contiguosly and
    ox,oy,oz being stored contiguosly as vectors of length
    vec_len.
    Called by: mobcal_ncoord
    Calls:     dcopy_
  */
  double *fx;
  double *fy;
  double *fz;
  double *ox;
  double *oy;
  double *oz;
  int incx;
  int inatom;

  int success;
  int copy_len;

  int vec_len;
  int iatom;

  int dljpot_input_stream_spacing;
  int ivec_len;

  success = 1;
  incx    = 1;
  inatom = state->inatom;
  vec_len = state->vec_len;
  dljpot_input_stream_spacing = 11 * _IVEC_LEN_;
  ivec_len = _IVEC_LEN_;
  ox     = state->ox;
  oy     = state->oy;
  oz     = state->oz;
  fx     = state->dljpot_input_stream;
  fy     = fx + _IVEC_LEN_;
  fz     = fy + _IVEC_LEN_;
  /*
  fx     = state->fx;
  copy_len = vec_len * 3;
  dcopy_(&copy_len,fx,&incx,ox,&incx);
  */
  for (iatom=0;iatom<vec_len;iatom+=_IVEC_LEN_) {
    dcopy_(&ivec_len,fx,&incx,ox,&incx);
    dcopy_(&ivec_len,fy,&incx,oy,&incx);
    dcopy_(&ivec_len,fz,&incx,oz,&incx);
    /*
      Caution the next six pointer assignments use address arithmetic.
    */
    ox += _IVEC_LEN_;
    oy += _IVEC_LEN_;
    oz += _IVEC_LEN_;
    fx += dljpot_input_stream_spacing;
    fy =  fx + _IVEC_LEN_;
    fz =  fy + _IVEC_LEN_;
  }
  return(success);
}
  
