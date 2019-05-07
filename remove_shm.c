#include "system_includes.h"
#include "remove_shm.h"
void remove_shm(int shmid) {
  struct shmid_ds shm_ds_buf;
  int shmctl_cmd;
  int shmctl_ret;
  int sav_err;
  int padi;
  shmctl_cmd = IPC_RMID;
  shmctl_ret = shmctl(shmid,shmctl_cmd,&shm_ds_buf);
  sav_err = errno;
  if (shmctl_ret < 0) {
    fprintf(stderr," Error return from shmctl command was %s\n",
	    strerror(sav_err));
    fflush(stderr);
  }
}

