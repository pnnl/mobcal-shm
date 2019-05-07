#ifndef _SETUP_SHM_H_
#define _SETUP_SHM_H_ 1
extern int setup_shm(const char *pathname, int proj_id, int *shmid_p,
		     size_t shm_size);
#endif
