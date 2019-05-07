#include "system_includes.h"
#include "mobcal_state_struct.h"
#include "mobcal_print_rotated_coords.h"
int mobcal_print_rotated_coords(struct mobcal_state_struct *state) {
  /*
    Print to the output file a rotated set of coordinates.
    Called by: mobcal_rotate and mobcal_struct_asym
    Calls:     fprintf,fflush
  */
  double *ox;
  double *oy;
  double *oz;
  double *fx;
  double *fy;
  double *fz;
  int inatom;
  int iatom;

  int dljpot_input_stream_spacing;
  int i;

  int success;
  int padi;
  FILE *ofp;
  FILE *efp;
  success = 1;
  ofp     = state->ofp;
  inatom  = state->inatom;
  ox      = state->ox;
  oy      = state->oy;
  oz      = state->oz;
  /*
  fx      = state->fx;
  fy      = state->fy;
  fz      = state->fz;
  */
  dljpot_input_stream_spacing = 11 * _IVEC_LEN_;
  fx      = state->dljpot_input_stream;
  /*
    Caution the following two pointer assigments do address arithmetic.
  */
  fy      = fx + _IVEC_LEN_;
  fz      = fy + _IVEC_LEN_;
  if (ofp) {
    fprintf(ofp,"        inititial coordinates                        new coordinates\n");
    for (iatom=0;iatom<inatom;iatom+=_IVEC_LEN_) {
      for (i=0;i<_IVEC_LEN_;i++) {
	fprintf(ofp,"%le %le %le     %le %le %le\n",
		ox[iatom+i],oy[iatom+i],oz[iatom+i],
		fx[i],fy[i],fz[i]);
      }
      /*
	Caution the following three pointer assigments do address arithmetic.
      */
      fx      += dljpot_input_stream_spacing;
      fy      = fx + _IVEC_LEN_;
      fz      = fy + _IVEC_LEN_;
    }
    fflush(ofp);
  }
  return(success);
}
      
