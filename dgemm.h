#ifndef _DGEMM_H_
#define _DGEMM_H_ 1
extern void dgemm_(char *transa_p, 
		  char *transb_p, 
		  int *m_p, 
		  int *n_p, 
		  int *k_p,
		  double *alpha_p, 
		  double *a, 
		  int *lda_p, 
		  double *b, 
		  int *ldb_p, 
		  double *beta_p, 
		  double *c, 
		  int *ldc_p);
#endif
