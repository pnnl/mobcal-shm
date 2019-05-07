#include "system_includes.h"
#include "lsame.h"
#include "blas.h"
#include "dgemm.h"
void dgemm_(char *transa_p, char *transb_p, int *m_p, int *n_p, int *k_p,
	   double *alpha_p, double *a, int *lda_p, double *b, 
	   int *ldb_p, double *beta_p, double *c, int *ldc_p) {
  /*
    Adapted from netlib dgemm.f
    Compute 
       c = alpha * op(A) * op(b) + beta * c
       op(x) = x or x^T
       
       op (a) is m by k
       op (b) is k by n
       c is m by n
       alpha and beta are scalars.
    Calls: lsame, dscal, daxpy, ddot
  */
  double alpha;
  double beta;
  double one;
  double zero;
  double temp;
  double *a_col;
  double *b_col;
  double *c_col;
  int64_t a_pos;
  int64_t b_pos;
  int64_t c_pos;

  int  m;
  int  n;

  int  k;
  int  lda;

  int  ldb;
  int  ldc;

  int  i;
  int  j;

  int  l;
  int  info;

  int  nota;
  int  notb;

  int nrowa;
  int ncola;
  
  int nrowb;
  int inc1;

  char transa;
  char transb;
  char n_char;
  char c_char;
  char t_char;
  char padc0;
  char padc1;
  char padc2;

  transa = *transa_p;
  transb = *transb_p;
  m      = *m_p;
  n      = *n_p;
  k      = *k_p;
  lda    = *lda_p;
  ldb    = *ldb_p;
  ldc    = *ldc_p;
  alpha  = *alpha_p;
  beta   = *beta_p;

  n_char  = 'N';
  c_char  = 'C';
  t_char  = 'T';
  one     = 1.0;
  zero    = 0.0;
  inc1    = 1;

  nota = lsame_(&transa,&n_char);
  notb = lsame_(&transb,&n_char);
   
  if (nota) {
    nrowa = m;
    ncola = k;
  } else {
    nrowa = k;
    ncola = m;
  }
  if (notb) {
    nrowb = k;
  } else {
    nrowb = n;
  }
  info = 0;
  if ((nota == 0)  && 
      (lsame_(&transa,&c_char) + lsame_(&transa,&t_char) == 0)) {
    info = 1;
  } else {
    if ((notb == 0) && 
	(lsame_(&transb,&c_char) + lsame_(&transb,&t_char) == 0)) {
      info = 2;
    } else {
      if (m < 0) {
	info = 3;
      } else {
	if (n < 0) {
	  info = 4;
	} else {
	  if (k < 0) {
	    info = 5;
	  } else {
	    if ((lda < 1) || (lda < nrowa)) {
	      info = 8;
	    } else {
	      if ((ldb < 1) || (ldb < nrowb)) {
		info = 10;
	      } else {
		if ((ldc < 1) || (ldc < m)) {
		  info = 13;
		}
	      }
	    }
	  }
	}
      }
    }
  }
  if (info != 0) {
    fprintf(stderr,"dgemm: Argument error: info = %d\n",info);
    fflush(stderr);
  } else {
    /*
      Arguments ok.
    */
    if (((m > 0) && (n > 0)) &&
	(((alpha != zero) && (k > 0)) || (beta != one))) {
      /* work to be done. */
      if (alpha == zero) {
	if (beta == zero) {
	  c_pos = 0;
	  for (j=0;j<n;j++) {
	    c_col = &c[c_pos];
	    for (i=0;i<m;i++) {
	      c_col[i] = zero;
	    }
	    c_pos += ldc;
	  }
	} else {
	  c_pos = 0;
	  for (j=0;j<n;j++) {
	    c_col = &c[c_pos];
	    for (i=0;i<m;i++) {
	      c_col[i] = beta*c_col[i];
	    }
	    c_pos += ldc;
	  }
	} /* end else (beta != 0) */
	/* end if (alpha == 0.0) */
      } else {
	/*
	  alpha != 0.0
	*/
	if (notb == 1) {
	  if (nota == 1) {
	    c_pos = 0;
	    b_pos = 0;
	    for (j=0;j<n;j++) {
	      c_col = &c[c_pos];
	      b_col = &b[b_pos];
	      if (beta == zero) {
		for (i=0;i<m;i++) {
		  c_col[i] = zero;
		}
	      } else {
		if (beta != one) {
		  dscal_(&m,&beta,c_col,&inc1);
		} 
	      }
	      a_pos = 0;
	      for (l = 0; l<k;l++) {
		a_col = &a[a_pos];
		temp = alpha * b_col[l];
		daxpy_(&m,&temp,a_col,&inc1,c_col,&inc1);
		a_pos += lda;
	      } /* end for (l...) */
	      c_pos += ldc;
	      b_pos += ldb;
	    } /* end (for j ...) */
	    /* end if nota */
	  } else {
	    /*
	      Form c = alpha*a^T * b + beta * c
	    */
	    c_pos = 0;
	    b_pos = 0;
	    for (j=0;j<n;j++) {
	      c_col = &c[c_pos];
	      b_col = &b[b_pos];
	      if (beta == zero) {
		for (i=0;i<m;i++) {
		  c_col[i] = zero;
		}
	      } else {
		if (beta != one) {
		  dscal_(&m,&beta,c_col,&inc1);
		} 
	      }
	      a_pos = 0;
	      for (i=0;i<m;i++) {
		a_col = &a[a_pos];
		temp = ddot_(&k,a_col,&inc1,b_col,&inc1);
		c_col[i] += alpha * temp;
		a_pos += lda;
	      }
	      b_pos += ldb;
	      c_pos += ldc;
	    }
	  } /* end else nota == 0 (use A^T) */
	} else {
	  /*
	    Use B*T
	  */
	  if (nota) {
	    /*
	      form  c= alpha * a * b^T + beta * c
	    */
	    c_pos = 0;
	    for (j=0;j<n;j++) {
	      c_col = &c[c_pos];
	      if (beta == zero) {
		for (i=0;i<m;i++) {
		  c_col[i] = 0;
		} 
	      } else {
		if (beta != one) {
		  dscal_(&m,&beta,c_col,&inc1);
		}
	      }
	      a_pos = 0;
	      b_pos = 0;
	      for (l=0;l<k;l++) {
		a_col = &a[a_pos];
		b_col = &b[b_pos];
		temp = alpha*b_col[j];
		daxpy_(&m,&temp,a_col,&inc1,c_col,&inc1);
		a_pos += lda;
		b_pos += ldb;
	      }
	      c_pos += ldc;
	    }
	    /* end  nota == 1 */
	  } else {
	    /*
	      form  c = alaph * a^T * b^T + beta * c
	    */
	    c_pos = 0;
	    for (j=0;j<n;j++) {
	      c_col = &c[c_pos];
	      if (beta == zero) {
		for (i=0;i<m;i++) {
		  c_col[i] = 0;
		} 
	      } else {
		if (beta != one) {
		  dscal_(&m,&beta,c_col,&inc1);
		}
	      }
	      a_pos = 0;
	      for (i=0;i<m;i++) {
		a_col = &a[a_pos];
		b_pos = 0;
		temp = zero;
		for (l=0;l<k;l++) {
		  b_col = &b[b_pos];
		  temp += a_col[l] * b_col[j];
		  b_pos += ldb;
		}
		c_col[i] += alpha * temp;
		a_pos += lda;
	      }
	      c_pos += ldc;
	    }
	  } /* end else nota == 0 */
	} /* end else notb == 0 */
      } /* end else alpha != 0; */
    } /* end if non trivial case; */
  } /* end else info == 0 */
} /* end dgemm */
