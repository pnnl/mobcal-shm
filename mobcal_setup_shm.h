#ifndef _MOBCAL_SETUP_SHM_H_
#define _MOBCAL_SETUP_SHM_H_ 1
extern int  mobcal_setup_shm(const char *pathname, 
			     int proj_id, int *shmid_p, 
			     size_t shm_size);
#endif
