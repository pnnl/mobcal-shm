#include "system_includes.h"
#include <sched.h>
#include "fork_and_bind.h"
int fork_and_bind(int num_threads) {
  /*
  returns a thread id of a bound forked process.
  */
  size_t    mask_size;
  cpu_set_t *cpu_maskp;

  int       max_threads;
  int       padi;
  int mythread_id;
  int not_done;
  int lchild_id;
  int rchild_id;

  
  pid_t pid;
  pid_t epid;

  mythread_id = 0;
  if (num_threads > 1) {
    not_done    = 1;
    while (not_done) {
      lchild_id = mythread_id + mythread_id + 1;
      if (lchild_id < num_threads) {
        pid = fork();
        if (pid == 0) {
    	/*
    	  I am the child.
    	*/
	  mythread_id = lchild_id;
        } else {
    	/*
    	  I am the parent, fork a right child if needed.
    	*/
	  rchild_id = lchild_id + 1;
	  if (rchild_id < num_threads) {
	    pid = fork();
	    if (pid == 0) {
	      /*
		I am the childe.
	      */
	      mythread_id = rchild_id;
	    } else {
	      /*
		I am the parent and I have forked two children.
	      */
	      not_done = 0;
	    }
	  } else {
	    /* 
	       I am the parent and no right child needed forking.
	    */
	    not_done = 0;
	  }
        }
      } else {
        /*
    	I am a leaf and need to fork no children.
        */
        not_done = 0;
      }
    } 
  }
  max_threads = 1024;
  cpu_maskp = CPU_ALLOC(max_threads);
  mask_size = CPU_ALLOC_SIZE(max_threads);
  CPU_ZERO_S(mask_size, cpu_maskp);
  CPU_SET_S(mythread_id,mask_size, cpu_maskp);

  pid = 0;
  sched_setaffinity(pid, mask_size, cpu_maskp);

  return(mythread_id);
}
