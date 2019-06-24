#include "system_includes.h"
#include "fork_and_bind.h"
#include "setup_shm.h"
#include "attach_shm.h"
#include <sched.h>
int  main (int argc, char **argv) {
  /*
    Local variables replicated by the fork.
  */
  struct shmid_ds shm_ds_buf;
  void *shm_base;
  int64_t shm_page_size;
  size_t  shm_size;
  size_t  mask_size;
  cpu_set_t *cpu_maskp;
  cpu_set_t cpu_mask;

  int max_threads;
  int num_threads;


  int my_thread_id;
  int bound_core;

  int shmid;
  int mobcal_projid;

  int shmflg;
  int sav_err;

  int shmctl_cmd;
  int shmctl_ret;
  
  int shmdt_ret;
  int i;

  int success;
  int padi;

  void *my_client_area;
  void *my_server_area;

  pid_t pid;
  pid_t epid;
  num_threads = 16;
  /*
    Now before forking we will want to set up shared memory
    on thread 0.
    The notion will be first to get the number of threads
    from the environment some how. For now we'll just have it
    return a fixed value.
    Then we will want a "server_area" and a "client_area" per thread.
    Each "area" will be a full page in size. 

    Some info that each processor will need then 
    num_threads,
    page_size,
    their_threadid,
    shmid

    We could do a separate shared memory segment for each "area"
    but that gets cumbersome and may not scale well.
    We will layout the pages in the shared memory segment with
    all the server pages consecutive followed by the client page.
    There is a separate server page for each client, and each client
    has its own client page.

    So we have 2*num_threads pages

    Client i's (i>=1) server page starts at shm_base + i*page_size and 
    ends at shm_base + page_size-1; 
    Client i's (i>=1) client page starts at 
       shm_base + (i+num_threads)* page_size and has length page_size.

    The server may only set values in the server pages 
    and each client is only allowed to set values in its own client page.
    strict adherence to the SWMR (Single Writer, Multiple Reader)  
    model is crucial

    Each client will need a pointer to 
    my_client_page = shm_base + ((my_thread_id + num_threads)) * page_size 
    my_server_page = shm_base + (my_thread_id * page_size)

  */
  shm_page_size = (int64_t)4096;
  shm_size  = (num_threads + num_threads) * shm_page_size;
  system("touch ./mobcal_shm");
  mobcal_projid = 14623516;
  success = setup_shm("./mobcal_shm",mobcal_projid,&shmid,shm_size);
  if (shmid >= 0) {
    my_thread_id = fork_and_bind(num_threads);
    max_threads = 1024;
    /*
      cpu_maskp = CPU_ALLOC(max_threads);
    */
    cpu_maskp = &cpu_mask;
    mask_size = CPU_ALLOC_SIZE(max_threads);
    pid = 0;
    if (sched_getaffinity(pid,mask_size,cpu_maskp)  == -1) {
      fprintf(stderr,"sched_getaffinity error\n");
      fflush(stderr);
    } else {
      for (i=0;i<num_threads;i++) {
	if (CPU_ISSET_S(i, mask_size, cpu_maskp)) {
	  bound_core = i;
	  break;
	}
      }
    }
    shm_base = attach_shm(my_thread_id,shmid);
    if (shm_base == ((void *)-1)) {
      fprintf(stderr," Error from shmat on %d was %s, shm_base = %p\n",my_thread_id,
	      strerror(sav_err),shm_base);
    } else {
      fprintf(stdout,"hello from %d, I think I'm bound to core %d\n",
	      my_thread_id,bound_core);
      fflush(stdout);
      /*
	Do work here.
      */
      /*
	Now we need to detach.
      */
      shmdt_ret = shmdt(shm_base);
      sav_err = errno;
      if (shmdt_ret < 0) {
	fprintf(stderr," Error return from shmdt was %s\n",
		strerror(sav_err));
      }
    }
    if (my_thread_id == 0) {
      /*
	If I'm thread 0, I should issue a call to destroy the shared
	memory segment after all have detached it.
	Can't issue this call till after all threads have attached,
	else their attach will fail with invalid argument.
	so first I need to make sure all have attached.
      */
      /*
      shmctl_cmd = IPC_RMID;
      shmctl_ret = shmctl(shmid,shmctl_cmd,&shm_ds_buf);
      sav_err = errno;
      if (shmctl_ret < 0) {
        fprintf(stderr," Error return from shmctl command was %s\n",
	        strerror(sav_err));
	fflush(stderr);
      }
      */
    }

  }
}
