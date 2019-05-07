#include "system_includes.h"
#include "mobcal_state_struct.h"
#include "mobcal_io_init.h"
int mobcal_io_init (struct mobcal_state_struct *state) {
  /*
    Open input, mfj_file, and output files output_file, returning 
    a 1 on success, 0 on failure, setting fp1, and fp2 fields
    int the state structure as their stream file pointers respecitivley.
    In the fortran code unit 8 corresponds to fp2, the output file,
    and unit 9 corresponds to fp1, the input file. 
    Called by: mobcal.
    Calls:     fopen,fprintf,fflush
    Uses: mfj_file and output_file fields of state
    Sets: ifp and ofp fields of state, the file pointers.
  */
  int success;
  int padi;
  success = 1;
  FILE *ifp;
  FILE *ofp;
  ifp = fopen(state->mfj_file,"r");
  if (ifp == NULL) {
    success = 0;
    fprintf(stderr,"mobcal_io_init: Error could not open %s for reading\n",
	   state->mfj_file);
    fflush(stderr);
  }
  ofp = NULL;
  if (success) {
    if (state->thread_id == 0) {
      ofp = fopen(state->output_file,"w");
      if (ofp == NULL) {    
	success = 0;
	fprintf(stderr,"mobcal_io_init: Error could not open %s for writing\n",
		state->output_file);
	fflush(stderr);
      }
    }
  }
  if (success) {
    state->ifp = ifp;
    state->ofp = ofp;
  }
  return(success);
}
