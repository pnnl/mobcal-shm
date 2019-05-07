#include "system_includes.h"

#include "blas.h"
void dscal_(int *nx, double *a_p, double *x, int *inc) {
  /*
    Dscal routine that ignores inc, assuming its 1.
    x <- a * x
    Called by: ode23tb
  */
  double a;
  int    i;
  int    padi;
  a = *a_p;
  for (i=0;i<(*nx);i++) {
    x[i] = a * x[i];
  }
}

void dcopy_(int *nx, double *x, int *incx, double *y, int *incy) {
  /*
    dcopy routine that ignores incx and incy assuming them both to be 1.
    y <- x
    Called by: ode_num_jac
  */
  int i;
  int padi;
  for (i=0;i<(*nx);i++) {
    y[i] = x[i];
  }
}

double dnrm2_(int *nx, double *x, int *inc) {
  /*
    Dnrm2 routine that ignores inc, assuming its 1.
    This still needs to be rewritten for stability and 
    scaling.
    Called by: ode23tb
    Calls:     sqrt
  */
  double sum;
  double twonorm;
  int i;
  int padi;
  sum = 0.0;
  for (i=0;i<(*nx);i++) {
    sum = sum + (x[i]*x[i]);
  }
  twonorm = sqrt(sum);
  return(twonorm);
}


void dgemv_(char *trans, int *m_p, int *n_p, double *alpha_p, double *a, 
	   int *lda_p, double *x, int *incx_p, double *beta_p, 
	   double *y, int *incy_p) {
   /*
     Simplified version of dgemv that assumes trans = 0,
     incx and incy are both 1.
     a is an m x n  matrix stored contiguously,
     x is a n x 1 vector, and y is an  m x 1 vector
     We compute y = beta * y + alpha *Ax using daxpy operations.
     Called by: ode23tb
     Calls:     daxpy
   */
  double *a_col;
  double alpha;
  double xv;
  int m;
  int n;
  int lda;
  int j;
  int incx;
  int incy;
  a_col = a;
  m   = *m_p;
  n   = *n_p;
  incx = 1;
  incy = 1;
  if ((*incx_p != 1) || (*incy_p != 1)) {
    fprintf(stderr,"dgemv_: Error, assumption of unit increments violated, You need a more general dgemv_ routine \n");
    fflush(stderr);
  } else {
    lda  = *lda_p;
    alpha = *alpha_p;
    dscal_(&m,beta_p,y,&incy);
    for (j=0;j<n;j++) {
      xv = alpha * x[j];
      daxpy_(&m,&xv,a_col,&incx,y,&incy);
      a_col += lda; /* Caution address arithmetic here */
    }
  }
}

void daxpy_(int *n_p, double *alpha_p, double *x, int *incx_p, 
	   double *y, int *incy_p) {
  /*
    daxpy routine from lapack where we assume incx and incy are both 1.
    Called by dgemv.
  */
  double alpha;
  int i;
  int n;
  n     = *n_p;
  alpha = *alpha_p;
  for (i=0;i<n;i++) {
    y[i] = y[i] + alpha * x[i];
  }
}

double ddot_(int *n_p, double *x, int *incx_p, double *y, int *incy_p) {
  double dtemp;
  double *xt;
  double *yt;
  int n;
  int incx;
  int incy;
  int i;
  n = *n_p;
  incx = *incx_p;
  incy = *incy_p;
  xt = x;
  yt = y;
  dtemp = 0.0;
  if ((incx == 1) && (incy == 1)) {
    for (i=0;i<n;i++) {
      dtemp += (*xt) * (*yt);
      xt += 1; /* Address arithmetic */
      yt += 1; /* Address arithmetic */
    }
  } else {
    if (incx < 0) {
      xt = &x[(1-n)*incx];
    }
    if (incy < 0) {
      yt = &y[(1-n)*incy];
    }
    for (i=0;i<n;i++) {
      dtemp += (*xt) * (*yt);
      xt += incx; /* Address arithmetic */
      yt += incy; /* Address arithmetic */
    }
  }
  return(dtemp);
}
int idamax_(int *n_p, double *dx, int *incx_p) {
  int n;
  int incx;
  int i;
  int ix;
  int max_loc;
  double dxmax;
  max_loc = -1;
  n = *n_p;
  incx = *incx_p;
  max_loc = -1;
  if (n > 0) {
    max_loc = 0;
    dxmax = fabs(dx[0]);
    if (incx > 1) {
      ix = incx;
      for (i=1;i<n;i++) {
	if (fabs(dx[ix]) > dxmax) {
	  dxmax = fabs(dx[ix]);
	  max_loc = i;
	}
	ix += incx;
      }
    } else {
      for (i=1;i<n;i++) {
	if (fabs(dx[i]) > dxmax) {
	  dxmax = fabs(dx[i]);
	  max_loc = i;
	}
      }
    }
  }
  return (max_loc+1);
}
