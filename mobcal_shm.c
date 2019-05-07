#include "system_includes.h"
#include "mobcal_read_parameters.h"
#include "get_num_threads.h"
#include "mobcal_shm_block_size.h"
#include "setup_shm.h"
#include "fork_and_bind.h"
#include "mobcal_shm_main.h"
#include "remove_shm.h"
int  main (int argc, char **argv) {
  size_t shm_size;
  char *filenames;
  int  *iparams;
  double *dparams;

  double dparams_b[8];
  char filenames_b[8192];
  int  iparams_b[32];

  int  my_thread_id;
  int  success;

  int  num_threads;
  int  shm_block_size;

  int  shmmid;
  int  shm_proj_id;

  int  itn;
  int  inp;

  int  ipos;
  int  padi;

  dparams   = (double *)&dparams_b[0];
  iparams   = (int*)&iparams_b[0];
  filenames = (char*)&filenames_b[0];

  success = mobcal_read_parameters(argc,argv,
				   filenames,iparams,dparams);
  if (success) {
    /*
    iparams[2] = numthreads;
    iparams[3] = shm_block_size (large enough to hold q1st,q2st, om*st arrays)
    iparams[4] = shm_size;
    iparams[5] = shmmid
    iparams[6] = 14623516;
    iparams[7] = mythread_id
    iparams[8] = i2;
    iparams[9] = ipr;
    iparams[10] = itn;
    iparams[11] = inp;
    iparams[12] = imp;
    iparams[13] = irn;
    iparams[14] = ibstmax;
    iparams[15] = igs;
    iparams[16] = im2;
    iparams[17] = im4;
    iparams[18] = ip;
    iparams[19] = it;
    iparams[20] = iu1;
    iparams[21] = iu2;
    iparams[22] = iu3;
    iparams[23] = iv;
    iparams[24] = use_dgt;
    iparams[25] = buffer_gas;
    iparams[26] = use_mt;
    iparams[27] = filename_limit;
    iparams[28] = use_iu_dljpot
    iparams[29] = avec_len;
    dparams[0]  = temp;
    dparams[1]  = m1;
    dparams[2]  = bg_dipole_mult;
    */
    itn = iparams[10];
    inp = iparams[11];

    num_threads = iparams[2];
    mobcal_shm_block_size(itn,inp,&shm_block_size);
    shm_size   = ((int64_t)2)*num_threads*shm_block_size;
    system("touch ./mobcal_shm");

    shm_proj_id = 14623516;
    /*
      set up a shared memory segment of size iparams[9], 
      storing the shmid in iparm[10], with iparams[11] being passed as
      the projid to the ftok call to generate the key for the shmget call.
    */
    success = setup_shm("./mobcal_shm",shm_proj_id,&shmmid,shm_size);

    ipos = 2;
    iparams[ipos] = num_threads;
    ipos += 1;
    iparams[ipos] = shm_block_size;
    ipos += 1;
    iparams[ipos] = (int)shm_size;
    ipos += 1;
    iparams[ipos] = shmmid;
    ipos += 1;
    iparams[ipos] = shm_proj_id;
    ipos += 1;
  }
  if (success) {
    /*
      Now we fork with iparams, dparams, and filenames as being
      the only state that needs to replicate.
    */
    my_thread_id = fork_and_bind(iparams[2]);
    dparams = (double *)&dparams_b[0];
    iparams = (int*)&iparams_b[0];
    filenames = (char*)&filenames_b[0];
    iparams[ipos] = my_thread_id;
    success = mobcal_shm_main(filenames,iparams,dparams);
    if (my_thread_id == 0) {
      /*
	Remove the shared memory segment.
      */
      remove_shm(shmmid);
    }
  }
  return(success);
}

