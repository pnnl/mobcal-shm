#include "system_includes.h"
#include "mobcal_state_struct.h"
#include "mobcal_unpack_params.h"
int mobcal_unpack_params(char *filenames, int *iparams, double *dparams,
			 struct mobcal_state_struct *state) {
  /*
    Unpack the parameters from filenames, iparams, and dparams into
    the state structure.
    Called by: mobcal_shm_main, mobcal
    Calls:     strcpy
  */
  char *param_file;
  char *at_param_file;
  char *mfj_file;
  char *output_file;

  int  success;
  int  ipos;

  int filename_limit;
  int padi;

  success = 1;
  ipos = 2;
  state->num_threads    = iparams[ipos];
  ipos += 1;
  state->shm_block_size = iparams[ipos];
  ipos += 1;
  state->shm_size      	= iparams[ipos];
  ipos += 1;
  state->shmid          = iparams[ipos];
  ipos += 1;
  state->shm_prjid     	= iparams[ipos];
  ipos += 1;
  state->thread_id     	= iparams[ipos];
  ipos += 1;
  state->i2  	       	= iparams[ipos];
  ipos += 1;
  state->ipr 	       	= iparams[ipos];
  ipos += 1;
  state->itn 	       	= iparams[ipos];
  ipos += 1;
  state->inp 	       	= iparams[ipos];
  ipos += 1;
  state->imp 	       	= iparams[ipos];
  ipos += 1;
  state->irn 	       	= iparams[ipos];
  ipos += 1;
  state->ibstmax       	= iparams[ipos];
  ipos += 1;
  state->igs            = iparams[ipos];
  ipos += 1;
  state->im2            = iparams[ipos];
  ipos += 1;
  state->im4            = iparams[ipos];
  ipos += 1;
  state->ip             = iparams[ipos];
  ipos += 1;
  state->it             = iparams[ipos];
  ipos += 1;
  state->iu1            = iparams[ipos];
  ipos += 1;
  state->iu2            = iparams[ipos];
  ipos += 1;
  state->iu3            = iparams[ipos];
  ipos += 1;
  state->iv             = iparams[ipos];
  ipos += 1;
  state->use_dgt        = iparams[ipos];
  ipos += 1;
  state->buffer_gas     = iparams[ipos];
  /*
#define DBG 1
  */
#ifdef DBG
    fprintf(stdout,"mobcal_unpack_params: ipos = %d, buffer_gas = %d\n",ipos,state->buffer_gas);
    fflush(stdout);
#endif    
  ipos += 1;
  state->use_mt         = iparams[ipos];
  ipos += 1;
  filename_limit        = iparams[ipos];
  state->filename_limit = filename_limit;
  ipos += 1;
  state->use_iu_dljpot  = iparams[ipos];
  ipos += 1;
  state->avec_len       = iparams[ipos];
  state->temp           = dparams[0];
  state->m1             = dparams[1];
  state->bg_dipole_mult = dparams[2];
  param_file = (char *)&filenames[0];
  at_param_file = (char *)&param_file[filename_limit];
  mfj_file      = (char *)&at_param_file[filename_limit];
  output_file   = (char *)&mfj_file[filename_limit];
  strcpy(state->param_file,param_file);
  strcpy(state->at_param_file,at_param_file);
  strcpy(state->mfj_file,mfj_file);
  strcpy(state->output_file,output_file);
  return(success);
}
