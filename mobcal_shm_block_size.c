#include "system_includes.h"
#include "mobcal_shm_block_size.h"
void mobcal_shm_block_size(int itn, int inp, int *shm_block_size_p) {
  int shm_block_size;
  int block_increment;

  int need;
  int padi;
  /*
    A shared memory block will have the following fields
    int64_t flag_1;
    int64_t flag_2;
    int64_t iic;
    int64_t ic;
    int64_t ig;
    int64_t im;
    int64_t spacer_1;
    int64_t spacer_2;
    double temp1[itn,inp]
    double temp2[itn,inp]

    Called by: mobcal_shm

  */
  block_increment = 4096;
  need = 64 + (2*inp*itn)*8;
  need = (need + block_increment - 1) / block_increment;
  shm_block_size = block_increment * need;
  *shm_block_size_p = shm_block_size;
}
