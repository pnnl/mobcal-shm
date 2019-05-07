#include "system_includes.h"
int  mobcal_setup_shm(const char *pathname, int proj_id, int *shmid_p, 
		      size_t shm_size) {
  /*
    Via call to ftok generate a key and and then
    with a cal to shmget create a shmid for passing to
    shmat.
  */
  key_t shm_key;
  key_t spare_key;
  int shmid;
  int success;
  int sav_err;
  int shmflg;
  int prot;
  success = 1;
  shm_key = ftok(pathname, proj_id);
  sav_err = errno;
  if (shm_key == -1) {
    success = 0;
    fprintf(stderr,"mobcal_setup_shm: ftok call failed with %s\n",
	    strerror(sav_err));
    fflush(stderr);
  } else {
    shmflg = IPC_CREAT;
    /*
      Set permissiona to be rwx------
    */
    prot = 0x180;
    shmflg = shmflg | prot;
    shmid = shmget(shm_key,shm_size,shmflg);
    sav_err = errno;
    if (shmid == -1) {
      success = 0;
      fprintf(stderr,"mobcal_setup_shm: shmget call failed with %s\n",
	      strerror(sav_err));
      fflush(stderr);
    }
  }
  *shmid_p = shmid;
  return(success);
}
