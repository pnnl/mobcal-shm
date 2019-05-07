#include "system_includes.h"
#include "mobcal_state_struct.h"
#include "mobcal_parse_atomtype_params.h"
#include "mobcal_read_coords_header.h"
#include "mobcal_alloc1.h"
#include "mobcal_read_coords_masses_charges.h"
#include "mobcal_set_properties.h"
#include "mobcal_compute_centroid.h"
#include "mobcal_save_pos.h"
#include "mobcal_shift_pos.h"
#include "mobcal_print_centered_coords.h"
#include "mobcal_struct_asym.h"
#include "mobcal_ncoord.h"
int mobcal_ncoord(struct mobcal_state_struct *state, int iic) {
  /*
    Read in coordinates and other parameters.
    This is a combination of fcoord and ncoord since they are very similar.
    The iic is the coordinate set, with iic = 0 being the first
    coordinate set. This makes for more consistent and managable 
    code source. The assymetry factor is put directly into the
    asympp element correspoinding to iic.
    Called by: mobcal
    Calls:     fprintf, fflush, fgets, sscanf,
               mobcal_parse_atomtype_params, (when iic = 0)
	       mobcal_read_coords_header, (when iic = 0)
               mobcal_alloc1,  (when iic = 0)
	       mobcal_read_coords_masses_charges,
	       mobcal_set_xmass_eolj_rolj_rhs,
	       mobcal_set_properties,
	       mobcal_compute_centroid,
	       mobcal_shift_pos,
	       mobcal_print_centered_coords, (when iic = 0)
	       mobcal_struct_asym	       
  */
  double m2;
  double mx;
  double *asympp;
  
  int  success;
  int  padi;

  FILE *ofp;
  FILE *efp;
  ofp    = state->ofp;
  success = 1;
  /*
    On the first iteration,
    we need to initialise the atomtype paremters, gmass,geolj,grolj,
    grhs, and gvalid fields of state. They have been allocated in
    mobcal_alloc0, and are need by mobcal_set_xmass_eolj_rolj_rhs.
  */
  if (iic == 0) {
    success = mobcal_parse_atomtype_params(state);
    if (success) {
      /*
	Read the header from the coordinates file setting
	the contents of the xlabel, units, and dchar fields of state,
	and setting the icoord, inatom, iunit, charge_dist, and 
	correct fields     of state.
      */
      success = mobcal_read_coords_header(state);
    }
    /*
      At this point we actually know the number of atoms , inatom,
      and the number of coordinate sets,
      so we should allocate those fields call mobcal_alloc1 from
      here instead of from mobcal.
    */
    /*
      Allocate space for the vector and fields of state, that
      depend on number of atoms or number of coordinate sets.
    */
    if (success) {
      success = mobcal_alloc1(state);
    }
  }
  if (success) {
    success = mobcal_read_coords_masses_charges(state);
  }
  if (success) {
    /*
      Set the xmaxx, eolj, rolj, eox4, rolj6, rolj12, dro6, dro12 rhs, and rhs2
      fields in the property vector and determine romax, the maximum rolj.
    */
    success = mobcal_set_properties(state, &mx);
  }
  if (success) {
    if (iic == 0) {
      state->m2 = mx;
      if (state->thread_id == 0) {
	if (ofp) {
	  fprintf(ofp,"mass of ion = %le\n",mx);
	  fflush(ofp);
	}
      }
    } else {
      m2 = state->m2;
      /*
	On invocations after the first compare mx to m2 to see they
	are the same.
      */
      if (mx != m2) {
	if (state->thread_id == 0) {
	  if (ofp) {
	    fprintf(ofp,"masses do not add up\n");
	    fflush(ofp);
	  }
	}
	success = 0;
      }
    } /* end else ic > 0 */
  } /* end if success */
  if (success) {
    /*
      Find the center of mass.
    */
    success = mobcal_compute_centroid(state);
  }
  if (success) {
    /*
      Shift the position coordinate by center of mass and scale
      by correction factor * 1e-10.
    */
    success = mobcal_shift_pos(state);
  }
  if(success) {
    success = mobcal_save_pos(state);
  }
  if (success) {
    if (iic == 0) {
      if (state->iu1 == 1) {
	if (state->thread_id == 0) {
	  success = mobcal_print_centered_coords(state);
	}
      }
    }
  }
  if (success) {
    /*
      determine the structural asymmetry parameter.
    */
    asympp = state->asympp;
    success = mobcal_struct_asym(state,&asympp[iic]);
  }
  return(success);
}
