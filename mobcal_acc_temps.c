#include "system_includes.h"
#include "mobcal_state_struct.h"
#include "mobcal_acc_temps.h"
void mobcal_acc_temps(struct mobcal_state_struct *state) {
  void *client_base;
  volatile void *left_area;
  volatile void *right_area;
  volatile void *my_area;
  volatile int64_t *jmy_area;
  volatile int64_t *jleft_area;
  volatile int64_t *jright_area;
  volatile int64_t *master_area;
  volatile double  *dmy_area;
  volatile double  *dleft_area;
  volatile double  *dright_area;
  volatile double  *my_temps;
  volatile double  *left_temps;
  volatile double  *right_temps;
  volatile int64_t my_flag2;
  volatile int64_t *left_flag2;
  volatile int64_t *right_flag2;
  volatile int64_t *master_flag2;
  int64_t left_spin_count;
  int64_t right_spin_count;
  int64_t finish_spin_count;
  int64_t one_l;

  int itn_inp;
  int itn_inp_2;

  int num_threads;
  int my_thread_id;

  int shm_block_size;
  int i;

  int left_child;
  int right_child;
  
  client_base    = state->client_area_base;
  my_thread_id   = state->thread_id;
  shm_block_size = state->shm_block_size;
  num_threads    = state->num_threads;
  itn_inp        = state->itn_inp;
  /*
    Wait  for left child to set his accumulated flag (flag2 field) to
  */
  master_area  = (int64_t *)client_base;
  one_l        = (int64_t) 1;
  my_area      = client_base + (my_thread_id * shm_block_size);
  jmy_area = (int64_t *)my_area;
  dmy_area = (double *)my_area;
  my_flag2     = jmy_area[_FLAG2_];
  left_child = my_thread_id + my_thread_id + 1;
  left_spin_count = 0;
  right_spin_count = 0;
  finish_spin_count = 0;
  if (left_child < num_threads) {
    /* Caution address arithmetic used in the statement below */
    left_area = client_base + (left_child * shm_block_size); 
    jleft_area = (int64_t *)left_area;
    dleft_area = (double *)left_area;
    left_flag2   = (int64_t *)&jleft_area[_FLAG2_];
    while (*left_flag2 == my_flag2) {
      left_spin_count += 1;
    }
    left_temps = (double*)&dleft_area[8];
    my_temps   = (double *)&dmy_area[8];
    itn_inp_2  = itn_inp + itn_inp;
    for (i=0;i<itn_inp_2;i++) {
      my_temps[i] += left_temps[i];
    }
    right_child = left_child + 1;
    if (right_child < num_threads) {
      right_area = client_base + (right_child * shm_block_size);
      jright_area = (int64_t *)right_area;
      dright_area = (double *)right_area;
      right_flag2 = (int64_t*)&jright_area[_FLAG2_];
      while (*right_flag2 == my_flag2) {
	right_spin_count += 1;
      }
      right_temps = (double*)&dright_area[8];
      for (i=0;i<itn_inp_2;i++) {
	my_temps[i] += right_temps[i];
      }
    }
  }
  /*
  jmy_area[_FLAG2_] = my_flag2;
  */
  my_flag2 += one_l;
  jmy_area[_FLAG2_] += one_l;
  if (my_thread_id != 0) {
    master_flag2 = &master_area[_FLAG2_];
    while (*master_flag2 != my_flag2) {
      finish_spin_count += 1;
    }
  }
  state->finish_spin_count = finish_spin_count;
  state->left_spin_count   = left_spin_count;
  state->right_spin_count  = right_spin_count;
}
