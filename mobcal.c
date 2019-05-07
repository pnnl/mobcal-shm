#include "system_includes.h"
#include "mobcal_state_struct.h"

#include "mobcal_read_parameters.h"
#include "mobcal_alloc0.h"
#include "mobcal_unpack_params.h"
#include "mobcal_io_init.h"
#include "mobcal_init_constants_0.h"
#include "mobcal_ncoord.h"
#include "mobcal_init_constants_1.h"
/*
#include "mobcal_xrand.h"
*/
#include "mobcal_mobil2.h"
#include "mobcal_print_summary.h"

int  main (int argc, char **argv) {
  struct mobcal_state_struct *state;
  char *filenames;
  int  *iparams;
  double *dparams;

  double *asympp;
  double *ehsm;
  double *tmm;
  double *tmc;
  double *ranlist;
  double sdevpc;
  char *log_file;
  char log_file_b[2048];

  double dparams_b[8];
  char filenames_b[8192];
  int  iparams_b[32];

  int immmin;
  int success;

  int icoord;
  int iic;

  int iic_p1;
  int imm;

  int thread_id;
  int num_threads;

  int filename_limit;
  int padi;

  FILE *ofp;
  FILE *lfp;
  /*
    Allocate storage for top level state struct.
    and its character fields.
  */
  dparams   = (double*)&dparams_b[0];
  iparams   = (int*)&iparams_b[0];
  filenames = (char*)&filenames_b[0];

  success = mobcal_read_parameters(argc,argv,
				   filenames,iparams,dparams);
  filename_limit = iparams[27];
  if (success) {
    success = mobcal_alloc0(&state,filename_limit);
  }
  /*
    This is the serial driver so we set thread_id to 0, and num_threads to 1.
  */
  if (success) {
    success = mobcal_unpack_params(filenames,iparams,dparams,state);
  }
  thread_id = 0;
  num_threads = 1;
  if (success) {
    state->thread_id = 0;
    state->num_threads = 1;
    log_file = (char*)&log_file_b[0];
    sprintf(log_file,"mobcal_log.0.%d",thread_id);
    lfp = fopen(log_file,"w");
    state->lfp = lfp;
  }
  /*
    Read in parameters and set num_data_points and num_parameters fields
    in state.
  */
  /*
    Open input and output files.
  */
  if (success) {
    success = mobcal_io_init(state);
  }
  /*
    Initialize some constants used by fcoord.
  */
  if (success) {
    ofp = state->ofp;
    success = mobcal_init_constants_0(state);
  }
  /*
    Read in the coordinates,
  */
  if (success) {
    iic     = 0;
    /*
      This call to ncoord with iic = 0 replaces fcoord call.
    */
    success = mobcal_ncoord(state,iic);
  }
  if (success) {
    /*
      Needs m2 field set by mobcal_ncoord call above.
    */
    success = mobcal_init_constants_1(state);
  }
  if (success) {
    ranlist = state->ranlist;
    asympp  = state->asympp;
    /*
    hold    = mobcal_xrand(state);
    */
    ehsm    = state->ehsm;
    tmm     = state->tmm;
    tmc     = state->tmc;
    ehsm[0] = 0.0;
    tmm[0]  = 0.0;
    immmin  = state->inor;
    state->imm = immmin;
    icoord  = state->icoord;
    for (iic = 0; iic<icoord; iic += 1) {
      iic_p1 = iic+1;
      state->iic = iic;
      if (ofp) {
	fprintf(ofp,"\n");
	if (icoord > 1) { 
	  fprintf(ofp," coordinate set = %d\n",iic);
	}
	fprintf(ofp,"structural asymmetry parameter = %le\n\n",asympp[iic]);
	fflush(ofp);
      }
      mobcal_mobil2(state,&tmm[iic],&tmc[iic],&sdevpc);
      imm = state->imm;
      if (imm < state->immmin) {
        state->immmin = imm;
      }
      if (imm > state->immmax) {
        state->immmax = imm;
      }
      state->im2 = 1;
      state->im4 = 1;
      if (iic != (icoord-1)) {
	mobcal_ncoord(state,iic_p1);
      }
    } /* end for (iic) */
    mobcal_print_summary(state,sdevpc);
  }
  return(0);
}
