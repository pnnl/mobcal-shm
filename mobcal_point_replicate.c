#include "system_includes.h"
#include "mobcal_point_replicate.h"
int mobcal_point_replictate(double *p, int n, int lda, double *pc) {
  /*
    Replicate the x,y,z coordinates stored in vector p[0],p[1],p[2],
    n times into vector pc separated by lda indices. lda must be >=3.
  */
  double *pc_pos;
  int success;
  int i;
  success = 1;
  if (n > 0) {
    if (lda >= 3) {
      pc_p = pc;
      for (i=0;i<n;i++) {
	pc_p[0] = p[0];
	pc_p[1] = p[1];
	pc_p[2] = p[2];
	pc_p += lda; /* Caution address arithmetic here. */
      }
    } else {
      success = 0;
      fprintf(stderr,"mobcal_point_replicate: Error lda < 3\n");
      fflush(stderr);
    } else {
      success = 0;
      fprintf(stderr,"mobcal_point_replicate: Error n < 1\n");
      fflush(stderr);
    }
    return success;
  }
