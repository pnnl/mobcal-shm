#include "system_includes.h"
#include "mobcal_state_struct.h"
#include "mobcal_alloc0.h"
#include "mobcal_unpack_params.h"
#include "attach_shm.h"
/*
#include "mobcal_read_parameters.h"
*/
#include "mobcal_io_init.h"
#include "mobcal_init_constants_0.h"
#include "mobcal_ncoord.h"
#include "mobcal_init_constants_1.h"
/*
#include "mobcal_xrand.h"
*/
#include "mobcal_shm_mobil2.h"
#include "mobcal_print_summary.h"
#include "mobcal_shm_main.h"
int mobcal_shm_main(char *filenames, int *iparams, double *dparams) {
  struct mobcal_state_struct *state;
  void   *shm_base;
  double *asympp;
  double hold;
  double *ehsm;
  double *tmm;
  double *tmc;
  double *ranlist;
  double sdevpc;
  double temp;

  char *log_file;
  char log_file_b[1024];

  int itn;
  int imp;

  int inp;
  int immmax;

  int immmin;
  int success;

  int icoord;
  int iic;

  int iic_p1;
  int imm;

  int shmid;
  int thread_id;

  int shm_block_size;
  int num_threads;

  int filename_limit;
  int padi;

  FILE *ofp;
  FILE *lfp;

  /*
#define DBG 1
  */

  /*
    Allocate storage for top level state struct.
    and its character fields.
  */
  filename_limit = iparams[27];
  success = mobcal_alloc0(&state,filename_limit);
  /*
    Here we realy need to open up a log file for each thread.
  */
  
  /*
    Read in parameters and set num_data_points and num_parameters fields
    in state.
  if (success) {
    success = mobcal_read_parameters(state,"mobcal.in");
  }
  */
  /*
    Unpack parameters into state.
  */
  if (success) {
    success = mobcal_unpack_params(filenames,iparams,dparams,state);
  }
  thread_id = state->thread_id;
  shmid     = state->shmid;
  shm_block_size = state->shm_block_size;
  num_threads    = state->num_threads;
  
  log_file = (char*)&log_file_b[0];
  sprintf(log_file,"mobcal_log.0.%d",thread_id);
#ifdef DBG
  lfp = fopen(log_file,"w");
  state->lfp = lfp;
  if (lfp) {
    fprintf(lfp,"Top of mobcal_shm_main thread_id = %d, num_threads = %d\n",
	    thread_id,num_threads);
    fflush(lfp);
  }
#else
  state->lfp = NULL;
#endif
  /*
    Attach to the shared memory segment.
  */
  shm_base = attach_shm(thread_id,shmid);
  state->server_area_base = (void*)shm_base;
  state->client_area_base = (void*)shm_base + 
    (shm_block_size * num_threads);
  state->my_server_area = (void*)shm_base + thread_id * shm_block_size;
  state->my_client_area = state->my_server_area + 
    (shm_block_size * num_threads);
#ifdef DBG
  if (lfp) {
    fprintf(lfp,"mobcal_shm_main: shm_base = %p, shm_block_size = %d\n",shm_base,shm_block_size);
    fprintf(lfp,"mobcal_shm_main: shm_base = %p\n",shm_base);
    fprintf(lfp,"state->server_area_base = %p\n",state->server_area_base);
    fprintf(lfp,"state->client_area_base = %p\n",state->client_area_base);
    fprintf(lfp,"state->my_server_area   = %p\n",state->my_server_area);
    fprintf(lfp,"state->my_client_area   = %p\n",state->my_client_area);
    fflush(lfp);
  }
#endif
  /*
    Open input and output files.
  */
  if (success) {
    success = mobcal_io_init(state);
  }
#ifdef DBG
  if (lfp) {
    fprintf(lfp,"mobcal_shm_main: after mocal_io_init, success = %d\n",
	    success);
    fflush(lfp);
  }
#endif
  /*
    Initialize some constants used by fcoord.
  */
  if (success) {
    ofp = state->ofp;
    success = mobcal_init_constants_0(state);
  }
#ifdef DBG
  if (lfp) {
    fprintf(lfp,"mobcal_shm_main: after mocal_init_constants_0, success = %d\n",
	    success);
    fflush(lfp);
  }
#endif
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
#ifdef DBG
  if (lfp) {
    fprintf(lfp,"mobcal_shm_main: after mocal_ncoord, success = %d\n",
	    success);
    fflush(lfp);
  }
#endif
  if (success) {
    /*
      Needs m2 field set by mobcal_ncoord call above.
    */
    success = mobcal_init_constants_1(state);
  }
#ifdef DBG
  if (lfp) {
    fprintf(lfp,
	    "mobcal_shm_main: after mocal_init_constants_1, success = %d\n",
	    success);
    fflush(lfp);
  }
#endif
  ranlist = state->ranlist;
  if (success) {
    asympp  = state->asympp;
    /*
    hold    = mobcal_xrand(state);
    */
    hold    = ranlist[0];
    ehsm    = state->ehsm;
    tmm     = state->tmm;
    tmc     = state->tmc;
    ehsm[0] = 0.0;
    tmm[0]  = 0.0;
    immmax  = 0;
    immmin  = state->inor;
    state->imm = immmin;
    itn     = state->itn;
    imp     = state->imp;
    inp     = state->inp;
    temp    = state->temp;
    icoord  = state->icoord;
    for (iic = 0; iic<icoord; iic += 1) {
      iic_p1 = iic+1;
      state->iic = iic;
      if (thread_id == 0) {
	if (ofp) {
	  fprintf(ofp,"\n");
	  if (icoord > 1) { 
	    fprintf(ofp," coordinate set = %d\n",iic);
	  }
	  fprintf(ofp,"structural asymmetry parameter = %le\n",asympp[iic]);
	  fflush(ofp);
	}
      }
#ifdef DBG
      if (lfp) {
	fprintf(lfp,"mobcal_shm_main: before mobil2 iic = %d\n",iic);
	fflush(lfp);
      }
#endif      
      mobcal_shm_mobil2(state,&tmm[iic],&tmc[iic],&sdevpc);
#ifdef DBG
      if (lfp) {
	fprintf(lfp,"mobcal_shm_main: after mobil2 iic = %d\n",iic);
	fflush(lfp);
      }
#endif      
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
#ifdef DBG
    if (lfp) {
      fprintf(lfp,"mobcal_shm_main: before mobcal_print_summary\n");
      fflush(lfp);
    }
#endif      
    if (thread_id == 0) {
      mobcal_print_summary(state,sdevpc);
    }
#ifdef DBG
    if (lfp) {
      fprintf(lfp,"mobcal_shm_main: after mobcal_print_summary\n");
      fflush(lfp);
    }
#endif      
  }
  return(0);
}
