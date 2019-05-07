#ifndef _DGEMV_H_
#define _DGEMV_H_ 1
extern void dgemv_(char *trans, int *m_p, int *n_p, double *alpha_p, double *a, 
		  int *lda_p, double *x, int *incx_p, double *beta_p, 
		  double *y, int *incy_p);
#endif
