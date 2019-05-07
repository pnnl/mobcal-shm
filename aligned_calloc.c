#include "system_includes.h"
#include "aligned_calloc.h"
void *aligned_calloc(int64_t nmemb, int64_t size_t, int64_t align_len) {
  int64_t mask;
  int64_t one_l;
  int64_t pad_len;
  int64_t tot_size;
  void *aligned_address;
  void *allocated_address;
  mask = align_len - 1;
  tot_size = (nmemb * size_t) + align_len;
  one_l = (int64_t)1;
  allocated_address = calloc(one_l,tot_size);
  pad_len = (align_len - (((int64_t)allocated_address) & mask)) & mask;
  aligned_address = allocated_address + pad_len;
  return(aligned_address);
}
  
