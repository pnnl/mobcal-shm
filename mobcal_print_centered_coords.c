#include "system_includes.h"
#include "mobcal_state_struct.h"
#include "mobcal_print_centered_coords.h"
int mobcal_print_centered_coords(struct mobcal_state_struct *state) {
  /*
    Print the scaled and shifted (to the centroid) coordinates of the
    mass.
    Called by: mobcal_ncoord
    Calls:     fprintf, fflush
    Uses fields ofp, fx,fy,fz,imass,pcharge,eolj,xe,rolj of state.
    Sets no fields.
  */
  double *fx;
  double *fy;
  double *fz;
  double *eolj;
  double *rolj;
  double *imass;
  double *pcharge;
  double xe;
  double xeolj;
  double xrolj;
  double ten_p10;
  double recip_xe;

  int    success;
  int    inatom;

  int    iatom;
  int    ximass;

  int    dljpot_input_stream_spacing;
  int    i;

  FILE *ofp;
  FILE *efp;
  success    = 1;
  ten_p10    = 1.0e10;
  inatom     = state->inatom;
  xe         = state->xe;
  /*
  fx         = state->fx;
  fy         = state->fy;
  fz         = state->fz;
  pcharge    = state->pcharge;
  */
  fx         = state->dljpot_input_stream;
  fy         = fx + _IVEC_LEN_;
  fz         = fy + _IVEC_LEN_;
  pcharge    = fz + _IVEC_LEN_ + _IVEC_LEN_;
  imass      = state->imass;
  eolj       = state->eolj;
  rolj       = state->rolj;
  ofp        = state->ofp;
  dljpot_input_stream_spacing = 11 * _IVEC_LEN_;
  if (xe == 0.0) {
    success = 0;
    fprintf(ofp,"mobcal_print_centerd_coords: Error xe was 0\n");
    fflush(ofp);
  }
  if (success) {
    fprintf(ofp,"initial_coordinates         mass   charge         LJ parameters\n");
    recip_xe = 1.0/xe;
    for (iatom = 0;iatom < inatom;iatom += _IVEC_LEN_) {
      for (i=0;i<_IVEC_LEN_;i++) {
	ximass = (int)imass[iatom+i];
	xeolj = eolj[iatom+i]*recip_xe;
	xrolj = rolj[iatom+i]*ten_p10;
	fprintf(ofp,"%le %le %le %d %le %le %le\n",
		fx[i],fy[i],fz[i],ximass,
		pcharge[i],xeolj,xrolj);
      }
      fx += dljpot_input_stream_spacing;
      fy = fx + _IVEC_LEN_;
      fz = fy + _IVEC_LEN_;
      pcharge = fz + _IVEC_LEN_ + _IVEC_LEN_;
    }
    fflush(ofp);
  }
  return(success);
}
