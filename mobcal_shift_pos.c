#include "system_includes.h"
#include "mobcal_state_struct.h"
#include "mobcal_shift_pos.h"
int mobcal_shift_pos(struct mobcal_state_struct *state) {
  /*
    Shift the position coordinate by center of mass and scale
    by correction factor * 1e-10.
    Called by: mobcal_ncoord.
    Uses fx. fy, fz, fxo, fyo, fzo, and correct fields of state.
    Sets the fx, fy, fz fields to shifted and scaled coordinates.
  */
  double *fx;
  double *fy;
  double *fz;
  double fxo;
  double fyo;
  double fzo;
  double correct;
  double ten_m10;
  double scale;
  int inatom;
  int iatom;
  int success;
  int dljpot_input_stream_spacing;
  int i;

  success = 1;
  inatom  = state->inatom;
  /*
  fx      = state->fx;
  fy      = state->fy;
  fz      = state->fz;
  */
  fxo     = state->fxo;
  fyo     = state->fyo;
  fzo     = state->fzo;
  correct = state->correct;
  ten_m10 = 1.0e-10;
  scale   = correct*ten_m10;
  dljpot_input_stream_spacing = 11 * _IVEC_LEN_;
  fx      = state->dljpot_input_stream;
  fy      = fx + _IVEC_LEN_;
  fz      = fy + _IVEC_LEN_;
  for (iatom = 0;iatom < inatom;iatom+=_IVEC_LEN_) {
    for (i=0;i<_IVEC_LEN_;i++) {
      fx[i] = (fx[i] - fxo) * scale;
      fy[i] = (fy[i] - fyo) * scale;
      fz[i] = (fz[i] - fzo) * scale;
    }
    fx += dljpot_input_stream_spacing;
    fy = fx + _IVEC_LEN_;
    fz = fy + _IVEC_LEN_;
  }
  return(success);
}
