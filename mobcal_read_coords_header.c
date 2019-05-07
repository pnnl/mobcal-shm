#include "system_includes.h"
#include "mobcal_state_struct.h"
#include "mobcal_read_coords_header.h"
int mobcal_read_coords_header(struct mobcal_state_struct *state) {
  /*
    Open the coordinates file for reading and 
    read in its header lines, 
    and echo some information to the output file.
    Sets the contents of the xlabel, units, and dchar fields of state.
    Sets the icoord, inatom, iunit, charge_dist, and correct fields
    of state.
    Called by: mobcal_ncoord.
    Calls:     fopen,fgets, fprintf, fflush, sscanf
  */
  double correct;
  char *line_buffer;
  char *fgets_p;
  char *mfj_file;
  char *unit;
  char *dchar;
  char *xlabel;
  int  line_len;
  int  success;

  int  icoord;
  int  inatom;

  int  charge_dist;
  int  nr;

  FILE *ofp;
  FILE *ifp;
  success     = 1;
  ofp         = state->ofp;
  ifp         = state->ifp;
  line_buffer = state->input_buffer;
  mfj_file    = state->mfj_file;
  unit        = state->unit;
  dchar       = state->dchar;
  xlabel      = state->xlabel;
  line_len    = state->input_buffer_len;
  if (ofp) {
    fprintf(ofp,"input file name = %s\n",mfj_file);
  }
  /*
    The input file is now opened in mobcal_io_init routine.
  */
  /*
  ifp = fopen(mfj_file,"r");
  if (ifp == NULL) {
    fprintf(ofp,"mobcal_read_coords_header: unable to open input file %s\n",
	    filen1);
    fflush(ofp);
    success = 0;
  }
  */
  if (success) {
    /*
    state->ifp = ifp;
    */
    fgets_p = fgets(line_buffer,line_len,ifp);
    if (fgets_p == NULL) {
      if (ofp) {
	success = 0;
	fprintf(ofp,"mobcal_read_coords_header: Error reading xlabel line\n");
	fflush(ofp);
      }
    }
  }
  if(success) {
    nr = sscanf(line_buffer,"%s",xlabel);
    if (nr != 1) {
      if (ofp) {
	success = 0;
	fprintf(ofp,"mobcal_read_coords_header: Error scanning xlabel\n");
	fflush(ofp);
      }
    } 
  }
  if (success) {
    if (ofp) {
      fprintf(ofp,"input file label = %s\n",xlabel);
    }
    fgets_p = fgets(line_buffer,line_len,ifp);
    if (fgets_p == NULL) {
      if (ofp) {
	success = 0;
	fprintf(ofp,"mobcal_read_coords_header: Error reading icoord line.\n");
	fflush(ofp);
      }
    }
  }
  if (success) {
    nr = sscanf(line_buffer,"%d",&icoord);
    if (nr != 1) {
      if (ofp) {
	success = 0;
	fprintf(ofp,"mobcal_read_coords_header: Error scanning icoord\n");
      }
    }
  }
  if (success) {
    state->icoord = icoord;
    state->num_parameters = icoord;
    if (ofp) {
      fprintf(ofp,"number of coordinate sets = %d\n",
	      icoord);
    }
    fgets_p = fgets(line_buffer,line_len,ifp);
    if (fgets_p == NULL) {
      if (ofp) {
	success = 0;
	fprintf(ofp,"mobcal_read_coords_header: Error reading inatom line.\n");
	fflush(ofp);
      }
    }
  }
  if (success) {
    nr = sscanf(line_buffer,"%d",&inatom);
    if (nr != 1) {
      if (ofp) {
	success = 0;
	fprintf(ofp,"mobcal_read_coords_header: Error scanning inatom\n");
      }
    }
  }
  if (success) {
    state->inatom = inatom;
    state->num_data_points = inatom;
    if (ofp) {
      fprintf(ofp,"number of atoms = %d\n",inatom);
    }
    fgets_p = fgets(line_buffer,line_len,ifp);
    if (fgets_p == NULL) {
      if (ofp) {
	success = 0;
	fprintf(ofp,"mobcal_read_coords_header: Error reading unit line.\n");
	fflush(ofp);
      }
    }
  }
  if (success) {
    nr = sscanf(line_buffer,"%s",unit);
    if (nr != 1) {
      if (ofp) {
	success = 0;
	fprintf(ofp,"mobcal_read_coords_header: Error scanning unit\n");
      }
    }
  }
  if (success) {
    state->iunit = -1;
    if (strncmp(unit,"au",2) == 0) {
      state->iunit = 0;
      if (ofp) {
	fprintf(ofp,"coordinates in atomic units\n");
      }
    } else {
      if (strncmp(unit,"ang",3) == 0) {
	state->iunit = 1;
	if (ofp) {
	  fprintf(ofp,"coordinates in angstroms\n");
	}
      } else {
	fprintf(ofp,"mobcal_read_coords_header: units not specified\n");
	fflush(ofp);
	success = 0;
      }
    }
  }
  if (success) {
    fgets_p = fgets(line_buffer,line_len,ifp);
    if (fgets_p == NULL) {
      if (ofp) {
	success = 0;
	fprintf(ofp,"mobcal_read_coords_header: Error reading dchar line.\n");
	fflush(ofp);
      }
    }
  }
  if (success) {
    nr = sscanf(line_buffer,"%s",dchar);
    if (nr != 1) {
      if (ofp) {
	success = 0;
	fprintf(ofp,"mobcal_read_coords_header: Error scanning dchar\n");
      }
    }
  }
  if (success) {
    charge_dist = -1;
    if (strncmp(dchar,"equal",5) == 0) {
      charge_dist = 1;
    } else {
      if (strncmp(dchar,"calc",4) == 0) {
	charge_dist = 2;
      } else {
	if (strncmp(dchar,"none",4) == 0) {
	  charge_dist = 0;
	}
      }
    }
    if (charge_dist < 0) {
      success = 0;
      if (ofp) {
	fprintf(ofp,"mobcal_read_coords_header: charge distribution not specified\n");
	fflush(ofp);
      }
    }
  }
  if (success) {
    state->charge_dist = charge_dist;
    fgets_p = fgets(line_buffer,line_len,ifp);
    if (fgets_p == NULL) {
      if (ofp) {
	success = 0;
	fprintf(ofp,"mobcal_read_coords_header: Error reading correction factor line.\n");
	fflush(ofp);
      }
    }
  }
  if (success) {
    nr = sscanf(line_buffer,"%le",&correct);
    if (nr != 1) {
      if (ofp) {
	success = 0;
	fprintf(ofp,"mobcal_read_coords_header: Error scanning correction factor\n");
      }
    }
  }
  if (success) {
    state->correct = correct;
  }
  return(success);
}
