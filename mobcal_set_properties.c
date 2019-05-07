#include "system_includes.h"
#include "mobcal_state_struct.h"
#include "mobcal_set_properties.h"
int mobcal_set_properties(struct mobcal_state_struct *state,
		          double *mass_sum_p) {
  /*
    Set the 
    xmass, eolj, rolj, and rhs vectors from the imass vector,
    field. Also set the eox4, ro6lj, ro12lj, dro6, dro12, rhs, and rhs2 
    and compute the mass of the ion, *mass_sum_p, and romax the
    maximum rolj
    Called by: mobcal_ncoord
    Calls:     fprintf, fflush
    Uses state fields:
      inatom,
      vec_len,
      romax_adj,
      gprops,
      imass,
      xmass,
      eolj,
      eox4,
      rolj,
      ro6lj,
      ro12lj,
      dro6,
      dro12,
      rhs,
      rhs2

    Sets state_fields:
      *xmass
      *eolj,
      *eox4,
      *rolj,
      *ro6lj,
      *ro12lj,
      *dro6,
      *dro12,
      *rhs,
      *rhs2,
      romax
    and *mass_sum_p

  */
  double *gprops;
  double *imass;
  double *xmass;
  double *eolj;
  double *eox4;
  double *rolj;
  double *ro6lj;
  double *ro12lj;
  double *dro6;
  double *dro12;
  double *rhs;
  double *rhs2;
  double *dljpot_input_stream;
  double m2;
  double cxmass;
  double xeolj;
  double xrolj;
  double xro2lj;
  double xeox4;
  double xro6lj;
  double xro12lj;
  double xdro6;
  double xdro12;
  double xrhs;
  double romax;
  double romax_adj;
  
  int inatom;
  int ixmass;

  int iatom;
  int success;

  int ipos;
  int iatom_p1;

  int vec_len;
  int i;

  int eox4_pos;
  int ro6lj_pos;

  int ro12lj_pos;
  int dro6_pos;

  int dro12_pos;
  int interval;

  int dljpot_input_stream_spacer;
  int ivec_len;
  
  int max_atom_imass;
  int padi;

  FILE *ofp;
  FILE *efp;
  success    = 1;
  inatom     = state->inatom;
  vec_len    = (int)state->vec_len;
  ofp        = state->ofp;
  gprops     = state->gprops;
  romax_adj  = state->romax_adj;
  imass      = state->imass;
  xmass      = state->xmass;
  eolj       = state->eolj;
  eox4       = state->eox4;
  rolj       = state->rolj;
  ro6lj      = state->ro6lj;
  ro12lj     = state->ro12lj;
  dro6       = state->dro6;
  dro12      = state->dro12;
  rhs        = state->rhs;
  rhs2       = state->rhs2;
  max_atom_imass = state->max_atom_imass;
  dljpot_input_stream = state->dljpot_input_stream;
  m2 = 0.0;
  ivec_len   = _IVEC_LEN_;

  romax = 0.0;
  interval = 0;
  eox4_pos = 6*ivec_len;
  ro6lj_pos = eox4_pos + ivec_len;
  ro12lj_pos = ro6lj_pos + ivec_len;
  dro6_pos = ro12lj_pos + ivec_len;
  dro12_pos = dro6_pos + ivec_len;
  dljpot_input_stream_spacer = 10*ivec_len;
  interval = 0;
  for (iatom = 0; ((iatom < inatom) && success); iatom++) {
    ixmass = (int) imass[iatom];
    ipos = ixmass << 2;
    if ((ixmass < 1) || (ixmass >= max_atom_imass)) {
      success = 0;
      iatom_p1 = iatom + 1;
      if (ofp) { 
	fprintf(ofp,"type not defined for atom number %d\n",iatom_p1);
      }
    } else {
      cxmass = gprops[ipos];
      if (cxmass == 0.0) {
	success = 0;
	iatom_p1 = iatom + 1;
	if (ofp) { 
	  fprintf(ofp,"type not defined for atom number %d\n",iatom_p1);
	}
      }
    }
    if (success) {
      xeolj         = gprops[ipos + 1];
      xrolj         = gprops[ipos + 2];
      xrhs          = gprops[ipos + 3];
      m2           += cxmass;
      if (xrolj > romax) {
	romax = xrolj;
      }
      xmass[iatom]        = cxmass;
      eolj[iatom]         = xeolj;
      rhs[iatom]          = xrhs;
      rolj[iatom]         = xrolj;

      /*
	The following fields should really be set by dscal and 
	other vectorized colls.
      */

      xeox4               = 4.0 * xeolj;
      xro2lj              = xrolj * xrolj;
      xro6lj              = xro2lj * xro2lj * xro2lj;
      xro12lj             = xro6lj * xro6lj;
      xdro6               = 6.0*xro6lj;
      xdro12              = 12.0*xro12lj;
      /*
      eox4[iatom]         = xeox4;
      ro6lj[iatom]        = xro6lj;
      ro12lj[iatom]       = xro12lj;
      dro6[iatom]         = xdro6;
      dro12[iatom]        = xdro12;
      */
      rhs2[iatom]         = xrhs * xrhs;
      dljpot_input_stream[eox4_pos] = xeox4;
      dljpot_input_stream[ro6lj_pos] = xro6lj;
      dljpot_input_stream[ro12lj_pos] = xro12lj;
      dljpot_input_stream[dro6_pos] = xdro6;
      dljpot_input_stream[dro12_pos] = xdro12;
      eox4_pos += 1;
      ro6lj_pos += 1;
      ro12lj_pos += 1;
      dro6_pos += 1;
      dro12_pos += 1;
      interval += 1;
      if (interval == ivec_len) {
	interval = 0;
	eox4_pos += dljpot_input_stream_spacer;
	ro6lj_pos += dljpot_input_stream_spacer;
	ro12lj_pos += dljpot_input_stream_spacer;
	dro6_pos += dljpot_input_stream_spacer;
	dro12_pos += dljpot_input_stream_spacer;
      }
    }
  } /* end for (iatom...) */
  /*
    Zero fill xmass, eolj,rhs,rolj, eox4, ro6lj,ro12lj,dro6,dro12,rhs2 vectors.
  */
  for (i=inatom;i<vec_len;i++) {
    eolj[i]   = 0.0;
    rhs[i]    = 0.0;
    rolj[i]   = 0.0;
    rhs2[i]   = 0.0;

    dljpot_input_stream[eox4_pos] = 0.0;
    dljpot_input_stream[ro6lj_pos] = 0.0;
    dljpot_input_stream[ro12lj_pos] = 0.0;
    dljpot_input_stream[dro6_pos] = 0.0;
    dljpot_input_stream[dro12_pos] = 0.0;
    eox4_pos += 1;
    ro6lj_pos += 1;
    ro12lj_pos += 1;
    dro6_pos += 1;
    dro12_pos += 1;
    interval += 1;
    if (interval == ivec_len) {
      interval = 0;
      eox4_pos += dljpot_input_stream_spacer;
      ro6lj_pos += dljpot_input_stream_spacer;
      ro12lj_pos += dljpot_input_stream_spacer;
      dro6_pos += dljpot_input_stream_spacer;
      dro12_pos += dljpot_input_stream_spacer;
    }
    /*
    eox4[i]   = 0.0;
    ro6lj[i]  = 0.0;
    ro12lj[i] = 0.0;
    dro6[i]   = 0.0;
    dro12[i]  = 0.0;
    */
  }
  if (success) {
    *mass_sum_p = m2;
    if (state->use_dgt) {
      romax = romax + romax_adj;
    }
    state->romax = romax;
  }
  return (success);
}
