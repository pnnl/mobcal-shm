#include "system_includes.h"
#include "blas.h"
#include "mobcal_state_struct.h"
#include "mobcal_read_coords_masses_charges.h"
int mobcal_read_coords_masses_charges(struct mobcal_state_struct *state) {
  /*
    Read from the input file the coordinates, masses, and charges of the
    atoms.
    Called by: mobcal_ncoord.
    Calls:     fgets, sscanf, fprintf, fflush
    Uses state_fields:
       inatom,
       vec_len, 
       xe,
       xeo,
       four_pi,
       pc,
       fx,
       fy,
       fz,
       pcharge,
       imass,
       iunit,
       charge_dist,
       mthree_pcharge,
       pc_scale_pcharge
    Sets state fields
    *fx,
    *fy,
    *fz,
    *pcharge,
    *imass
  */
  double *fx;
  double *fy;
  double *fz;
  double *pcharge;
  double *imass;
  double *mthree_pcharge;
  double *pc_scale_pcharge;
  double *dljpot_input_stream;
  double xe;
  double xeo;
  double four_pi;
  double pc;
  double xe2_o_4pi_xeo;
  double pc_scale;
  double mthree;
  double xpcharge;
  double ximass;
  double xfx;
  double xfy;
  double xfz;
  double tcharge;
  double acharge;
  double unit_multiplier;
  double fxi;
  double fyi;
  double fzi;
    
  char *fgets_p;
  char *line_buffer;

  int success;
  int inatom;

  int iatom;
  int iatom_p1;

  int nr;
  int iunit;

  int charge_dist;
  int line_len;

  int vec_len;
  int i; 

  int incx;
  int interval;


  int fx_pos;
  int fy_pos;

  int fz_pos;
  int pc_scale_pcharge_pos;

  int pcharge_pos;
  int mthree_pcharge_pos;

  int dljpot_input_stream_spacer;
  int ivec_len;

  FILE *ifp;
  FILE *ofp;
  
  success     	   = 1;
  ivec_len         = _IVEC_LEN_;
  line_buffer 	   = state->input_buffer;
  line_len    	   = state->input_buffer_len;
  inatom      	   = state->inatom;
  vec_len     	   = (int)state->vec_len;
  xe          	   = state->xe;
  xeo         	   = state->xeo;
  four_pi     	   = state->four_pi;
  pc          	   = state->pc;
  ifp         	   = state->ifp;
  ofp         	   = state->ofp;
  fx          	   = state->fx;
  fy          	   = state->fy;
  fz          	   = state->fz;
  imass       	   = state->imass;
  pcharge     	   = state->pcharge;
  mthree_pcharge   = state->mthree_pcharge;
  pc_scale_pcharge = state->pc_scale_pcharge;
  dljpot_input_stream = state->dljpot_input_stream;
  tcharge          = 0.0;
  acharge          = 0.0;
  iunit       = state->iunit;
  charge_dist = state->charge_dist;
  mthree = -3.0;
  xe2_o_4pi_xeo = (xe * xe) / (four_pi * xeo);
  pc_scale = pc * xe2_o_4pi_xeo;
  incx = 1;
  if (iunit == 0) {
    unit_multiplier = 0.52917706;
  } else {
    unit_multiplier = 1.0;
  }
  dljpot_input_stream_spacer = 10 * _IVEC_LEN_;
  fx_pos = 0;
  fy_pos = fx_pos + _IVEC_LEN_;
  fz_pos = fy_pos + _IVEC_LEN_;
  pc_scale_pcharge_pos = fz_pos + _IVEC_LEN_;
  pcharge_pos = pc_scale_pcharge_pos + _IVEC_LEN_;
  mthree_pcharge_pos = pcharge_pos + _IVEC_LEN_;
  interval = 0;
  for (iatom = 0;iatom < inatom;iatom++) {
    fgets_p = fgets(line_buffer,line_len,ifp);
    if (fgets_p == NULL) {
      success = 0;
      if (ofp) {
	iatom_p1 = iatom + 1;
	fprintf(ofp,"mobcal_read_coords_masses_charges: Error reading coordinate line %d.\n",iatom_p1);
	fflush(ofp);
      }
      success = 0;
    }
    if (success) {
      nr = sscanf(line_buffer,"%le %le %le %le %le",
		  &xfx,&xfy,&xfz,&ximass,&xpcharge);
      if (nr != 5) {
	success = 0;
	if (ofp) {
	  iatom_p1 = iatom + 1;
	
	  fprintf(ofp,"mobcal_read_coords_masses_charges: Error scannind coordinate line %d\n",iatom_p1);
	  fflush(ofp);
	}
	break;
      }
    }
    if (success) {
      fxi          = xfx * unit_multiplier;
      fx[iatom]    = fxi;
      fyi          = xfy * unit_multiplier;
      fy[iatom]    = fyi;
      fzi          = xfz * unit_multiplier;
      fz[iatom]    = fzi;
      imass[iatom] = (double)((int)(ximass + .5));
      pcharge[iatom] = xpcharge;
      tcharge += xpcharge;
      acharge += fabs(xpcharge);
      dljpot_input_stream[fx_pos] = fxi;
      dljpot_input_stream[fy_pos] = fyi;
      dljpot_input_stream[fz_pos] = fzi;
      dljpot_input_stream[pc_scale_pcharge_pos] = pc_scale * xpcharge;
      dljpot_input_stream[pcharge_pos] = xpcharge;
      dljpot_input_stream[mthree_pcharge_pos] = mthree * xpcharge;
    }
    interval += 1;
    fx_pos += 1;
    fy_pos += 1;
    fz_pos += 1;
    pc_scale_pcharge_pos += 1;
    pcharge_pos += 1;
    mthree_pcharge_pos += 1;
    if (interval == _IVEC_LEN_) {
      interval = 0;
      fx_pos += dljpot_input_stream_spacer;
      fy_pos += dljpot_input_stream_spacer;
      fz_pos += dljpot_input_stream_spacer;
      pc_scale_pcharge_pos += dljpot_input_stream_spacer;
      pcharge_pos += dljpot_input_stream_spacer;
      mthree_pcharge_pos += dljpot_input_stream_spacer;
    }
  } /* end for iatom */
  /*
    Reset particle charges if equal or none options were chosen.
  */
  if (charge_dist == 1) {
    /*
      Equal charge distribution.
    */
    xpcharge = tcharge/((double)inatom);
    pc_scale_pcharge_pos = 3 * _IVEC_LEN_;
    pcharge_pos = pc_scale_pcharge_pos + _IVEC_LEN_;
    mthree_pcharge_pos = pcharge_pos + _IVEC_LEN_;
    interval = 0;
    for (iatom=0;iatom<inatom;iatom++) {
      pcharge[iatom] = xpcharge;
      dljpot_input_stream[pc_scale_pcharge_pos] = pc_scale * xpcharge;
      dljpot_input_stream[pcharge_pos]          = xpcharge;
      dljpot_input_stream[mthree_pcharge_pos]   = mthree * xpcharge;
      interval += 1;
      pc_scale_pcharge_pos += 1;
      pcharge_pos += 1;
      mthree_pcharge_pos += 1;
      if (interval == _IVEC_LEN_) {
	interval = 0;
	pc_scale_pcharge_pos += dljpot_input_stream_spacer;
	pcharge_pos += dljpot_input_stream_spacer;
	mthree_pcharge_pos += dljpot_input_stream_spacer;
      }
    }
  } else {
    if (charge_dist == 0) {
      /*
	No charge distributed.
      */
      for (iatom=0;iatom<inatom;iatom++) {
	pcharge[iatom] = 0.0;
	dljpot_input_stream[pc_scale_pcharge_pos] = 0.0;
	dljpot_input_stream[pcharge_pos]          = 0.0;
	dljpot_input_stream[mthree_pcharge_pos]   = 0.0;
	interval += 1;
	pc_scale_pcharge_pos += 1;
	pcharge_pos += 1;
	mthree_pcharge_pos += 1;
	if (interval == _IVEC_LEN_) {
	  interval = 0;
	  pc_scale_pcharge_pos += dljpot_input_stream_spacer;
	  pcharge_pos += dljpot_input_stream_spacer;
	  mthree_pcharge_pos += dljpot_input_stream_spacer;
	}
      }
    }
  }
  /*
    Zero fill fx,fy, fz, imass, and pcharge vectors.
  */
  
  for (i=inatom;i<vec_len;i++) {
    fx[i] = 0.0;
    fy[i] = 0.0;
    fz[i] = 0.0;
    imass[i] = 0.0;
    pcharge[i] = 0.0;
    dljpot_input_stream[fx_pos] = 0.0;
    dljpot_input_stream[fy_pos] = 0.0;
    dljpot_input_stream[fz_pos] = 0.0;
    dljpot_input_stream[pc_scale_pcharge_pos] = 0.0;
    dljpot_input_stream[pcharge_pos] = 0.0;
    dljpot_input_stream[mthree_pcharge_pos] = 0.0;
    interval += 1;
    fx_pos += 1;
    fy_pos += 1;
    fz_pos += 1;
    pc_scale_pcharge_pos += 1;
    pcharge_pos += 1;
    mthree_pcharge_pos += 1;
    if (interval == _IVEC_LEN_) {
      interval = 0;
      fx_pos += dljpot_input_stream_spacer;
      fy_pos += dljpot_input_stream_spacer;
      fz_pos += dljpot_input_stream_spacer;
      pc_scale_pcharge_pos += dljpot_input_stream_spacer;
      pcharge_pos += dljpot_input_stream_spacer;
      mthree_pcharge_pos += dljpot_input_stream_spacer;
    }
  }
  /*
    Generate the scaled vectors of pcharge, mthree_p_charge,
    pc_scale_pcharge.
  */
  dcopy_(&vec_len,pcharge,&incx,mthree_pcharge,&incx);
  dscal_(&vec_len,&mthree,mthree_pcharge,&incx);
  dcopy_(&vec_len,pcharge,&incx,pc_scale_pcharge,&incx);
  dscal_(&vec_len,&pc_scale,pc_scale_pcharge,&incx);
  if (ofp) {
    fprintf(ofp,"total charge = %le\ntotal absolute charge = %le\n",
	    tcharge,acharge);
  }
  return(success);
}
