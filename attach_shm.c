#include "system_includes.h"
#include "attach_shm.h"
void *attach_shm(int my_thread_id, int shmid) {
  /*  
#define DBG 1
  */
  void *shm_suggest;
  void *shm_base;
  int shmflg;
  int sav_err;
  shm_suggest = NULL;
  shmflg = 0;
#ifdef DBG
  fprintf(stdout,
	  "attach_shm: my_thread_id = %d, shmid = %d, shmflg = %d, shm_suggest = %p\n",
	  my_thread_id,shmid,shmflg,shm_suggest);
  fflush(stdout);
#endif
  shm_base = shmat(shmid,shm_suggest,shmflg);
  sav_err = errno;
  if (shm_base == ((void *)-1)) {
    fprintf(stderr,"attach_shm:  Error from shmat on %d was %s, shm_base = %p\n",
	    my_thread_id,strerror(sav_err),shm_base);
  } 
#ifdef DBG
  fprintf(stdout,"attach_shm: my_thread_id = %d, shm_base = %p\n",my_thread_id,shm_base);
  fflush(stdout);
#endif
  return(shm_base);
}
