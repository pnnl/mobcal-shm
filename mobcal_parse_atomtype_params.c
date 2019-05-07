#include "system_includes.h"
#include "mobcal_state_struct.h"
#include "mobcal_parse_atomtype_params.h"
int mobcal_parse_atomtype_params(struct mobcal_state_struct *state) {
  /*
    read the mass, eolj, rolj and rhs parameters for each atom type
    from the atomtype_parameters.in file.
    Called by: mobcal_ncoord
    Calls:     strcpy,fopen,fclose,fgets, fprintf, fflush, sscanf

    As the xmass, eolj, rolj, and rhs fields of state 
    are all based on the
    same calculation dependent on the imass for each atom
    we will have a gprops vextor with
    the xmass, eolj, rolj, and rhs values for each integer mass stored 
    consecutively. A zero xmass filed in the gprops array is taken as 
    an invalid integer weight.
    We fill these vectors for their values for imasses = 
    1(H), 12(C), 14(N), 16(O), 19(F), 23(Na), 28(Si), 31(P),
    32(S), 35(Cl), 56(Fe)
  */
  double *gprops;
  double xgmass;
  double xgrhs;
  double xgeolj;
  double xgrolj;
  double xe;
  double eogas;
  double rogas;
  double conve;
  double convr;
  double conve_xe;
  double convrs;
  double ten_m10;
  /*
    NB do match the single precison constants used in the fortran code,
    We declare the xgeolj and xgrolj vars as floats as they don't have
    d0's after them in the fortran mobcal, nope this was just a translation
    mistake. Original IU code has d0's.
  */
  char atp_file_b[128];
  char *atp_file;
  char *read_p;
  char *input_buffer;
  char *atom_chars;
  char atom_chars_b[8];

  int success;
  int i;

  int input_buffer_len;
  int ipos;

  int nr;
  int gprops_len;

  int gprops_len_m3;
  int max_atom_imass;

  int use_dgt;
  int padi;

  FILE *atp_fp;
  FILE *efp;

  success  = 1;
  atp_file = (char*)&atp_file_b[0];
  xe       = state->xe;
  use_dgt  = state->use_dgt;
  
  atp_fp = fopen(state->at_param_file,"r");
  if (atp_fp == NULL) {
    success = 0;
    fprintf(stderr,"mobcal_parse_atomtype_params: Error unable to open %s\n",
	    state->at_param_file);
    fflush(stderr);
  }
  if (success) {
    eogas = 0.06900;
    rogas = 3.6600;
    conve = (4.2*0.01036427);
    convr = 0.890898718;
    conve_xe = conve * xe;
    ten_m10 = 1.0e-10;
    convrs = convr*ten_m10;
    input_buffer     = state->input_buffer;
    input_buffer_len = state->input_buffer_len;
    atom_chars       = (char*)&atom_chars_b[0];
    gprops           = state->gprops;
    max_atom_imass   = state->max_atom_imass;
    gprops_len = max_atom_imass << 2;
    for (i =0;i < gprops_len;i++) {
      gprops[i] = 0.0;
    }
    /*
      skip the first header line.
    */
    read_p = fgets(input_buffer,input_buffer_len,atp_fp);
    if (read_p != input_buffer) {
      success = 0;
      fprintf(stderr,"mobcal_parse_atomtype_params: Error skipping header line\n");
      fflush(stderr);
    }
  }
  if (success) {
    read_p = fgets(input_buffer,input_buffer_len,atp_fp);
    gprops_len_m3 = gprops_len - 3;
    while (!feof(atp_fp)) {
      nr = sscanf(input_buffer,"%s %le %le %le %le",
		  atom_chars,&xgmass,&xgeolj,&xgrolj,&xgrhs);
      if (nr != 5) {
	fprintf(stderr,"mobcal_parse_atomtype_params: Error nalformed"
		" atom line was\n%s\n",input_buffer);
	fflush(stderr);
      } else {
	ipos = (int)(xgmass + .5);
	ipos = ipos << 2;
	if ((ipos > 0) && (ipos < gprops_len_m3)) {
	  gprops[ipos]   = xgmass;
	  /*
	  gprops[ipos+1] = sqrt(eogas * xgeolj) * conve_xe; 
	  gprops[ipos+2] = sqrt(rogas * xgrolj) * convrs; 
	  */
	  if (use_dgt) {
	    gprops[ipos+1] = sqrt(eogas * xgeolj) * conve_xe; 
	    gprops[ipos+2] = sqrt(rogas * xgrolj) * convrs; 
	  } else {
	    gprops[ipos+1] = xgeolj * xe; 
	    gprops[ipos+2] = xgrolj * ten_m10;
	  }
	  gprops[ipos+3] = xgrhs*ten_m10; /*rhs */
	}
      }
      read_p = fgets(input_buffer,input_buffer_len,atp_fp);
    }
    fclose(atp_fp);
  }
  return(success);
}
