#include "system_includes.h"
#include "blas.h"
#include "mobcal_state_struct.h"
#include "mobcal_deriv.h"
#include "mobcal_diffeq.h"
void mobcal_diffeq(struct mobcal_state_struct *state,
		   int *l_p,
		   double *tim_p,
		   double *dt_p,
		   double *w,
		   double *dw,
		   double *pot_p,
		   double *dmax_p) {
  /*
    Integration subroutine - uses 5th order runge-kutta-gill to 
    initiate and 5th order adams-moulton predictor-corrector to 
    propagate. Parameter l is initially set to zero and then 
    incremented to tell the subroutine when to switch between 
    integration methods. DIFFEQ calls subroutine DERIV to define 
    the equations of motion to be integrated.
    Called by: mobcal_gsang
    Calls:     mobcal_deriv, dscal, dcopy, daxpy, dgemv, memmove
  */
  struct mobcal_diffeq_state_struct *diffeq_state;
  double tim;
  double dt;
  double pot;
  double dmax;
  double *a;
  double *b;
  double *c;
  double *three_a;
  double *ampc;
  double *amcc;
  double *savw;
  double *savdw;
  double *q;
  double *r;
  double *array;
  double var;
  double cvar;
  double acst;

  double hvar;
  double hcvar;

  double alpha;
  double beta;
  double mbj;

  double *predictor;
  double *corrector;
  double *array_col_1;
  double *array_col_4;

  double half;
  double two;
  double zero;

  int64_t move_size;
  int l;
  int j;

  int k;
  int lindex;

  int ifour;
  int ifive;

  int isix;
  int incx;

  int ld_array;
  int four_ld_array;

  char trans;
  

  l    = *l_p;
  dt   = *dt_p;
  tim  = *tim_p;
  pot  = *pot_p;
  dmax = *dmax_p;
  half  = 0.5e0;
  two   = 2.0e0;
  zero  = 0.0e0;
  ifive = 5;
  ifour = 4;
  isix  = 6;
  incx  = 1;

  diffeq_state = state->diffeq_state;
  /*
    Fields of diffeq_state are allocated in mobca_alloc0, and
    initialized in mobcal_init_constants_0.
  */
  a = diffeq_state->a;
  b = diffeq_state->b;
  c = diffeq_state->c;
  three_a = diffeq_state->three_a;
  q = diffeq_state->q;
  r = diffeq_state->r;
  ampc      	= diffeq_state->ampc;
  amcc      	= diffeq_state->amcc;
  var       	= diffeq_state->var;
  cvar      	= diffeq_state->cvar;
  acst      	= diffeq_state->acst;
  array     	= diffeq_state->array;
  savw      	= diffeq_state->savw;
  savdw     	= diffeq_state->savdw;
  ld_array  	= diffeq_state->ld_array;
  four_ld_array = diffeq_state->four_ld_array;
  predictor    	= diffeq_state->predictor;
  corrector    	= diffeq_state->corrector;
  array_col_4  	= diffeq_state->col_4;
  array_col_1  	= diffeq_state->col_1;

  if (l == 0) {
    /* statement 1 */
    for (j=0;j<5;j++) {
      q[j] = zero;
    }
    diffeq_state->hvar = dt * var;
    diffeq_state->hcvar = dt * cvar;
    dt    = half * dt;
  } 
  if (l >= 0) {
    /* statement 3 */
    l = l + 1;
    /*
      This is the runge-kutta-gill part...the steps are broken up into
      half steps to improve accuracy.
    */
    for (k=0;k<2;k++) {  
/*
   15 do 7 j=1,4
      if ((-1)**j.gt.0) tim=tim+0.5*dt
c	write(8,*) 'j = ', j, ', dt = ', dt
      call deriv(w,dw,pot,dpotx,dpoty,dpotz,dmax)
      do 7 i=1,6
      dw(i)=dt*dw(i)
      r=a(j)*(dw(i)-b(j)*q(i))
      w(i)=w(i)+r
    7 q(i)=q(i)+3.0*r+c(j)*dw(i)
*/
      for (j=0;j<4;j++) {
	if (j & 1) {
	  tim = tim + .5 * dt;
	}
	mobcal_deriv(state,w,dw,&pot,&dmax);
	mbj = -b[j];
	/*
	  dw -> dt * dw
	*/
	dscal_(&isix,&dt,dw,&incx);
	/*
	  r <- dw
	*/
	dcopy_(&isix,dw,&incx,r,&incx);
	/*
	  r <- dw - bj * q
	*/
	daxpy_(&isix,&mbj,q,&incx,r,&incx);
	/*
	  w <- w + aj *r
	*/
	daxpy_(&isix,&a[j],r,&incx,w,&incx);
	/*
	  q <- q + 3*aj*r 
	*/
	daxpy_(&isix,&three_a[j],r,&incx,q,&incx);
	/*
	  q <- q + cj * dw 
	*/
	daxpy_(&isix,&c[j],dw,&incx,q,&incx);

      }
      mobcal_deriv(state,w,dw,&pot,&dmax);
    } /* end for (k ...) */
    if (l < 6) {
      /*
	Coopy dw to column l of arrau
      */
      lindex = (l-1) * ld_array;
      dcopy_(&isix,dw,&incx,&array[lindex],&incx);
    } else {
      l = -1;
      dt = two * dt;
    }
  } else {
    /*
      l < 0
      This is the adams-moulton predictor-corrector part.
    */
    /*
    4 do 10 j=1,6
      savw(j)=w(j)
      savdw(j)=dw(j)
      array(6,j)=savdw(j)
      do 9 i=1,5
    9 array(6,j)=array(6,j)+ampc(i)*array(i,j)
   10 w(j)=array(6,j)*hvar+w(j)
      tim=tim+dt
    */
    hcvar = diffeq_state->hcvar;
    hvar  = diffeq_state->hvar;
    dcopy_(&isix,w,&incx,savw,&incx);
    dcopy_(&isix,dw,&incx,savdw,&incx);

    /*
      Column 5 (the predictor) of array gets savdw.
    */
    dcopy_(&isix,savdw,&incx,predictor,&incx);
    trans = 'N';
    beta = 1.0;
    alpha = 1.0;
    /*
      predictor += array(*,0:4) * ampc(0:4)'
    */
    dgemv_(&trans,&isix,&ifive,&alpha,array,&ld_array,ampc,&incx,
	   &beta,predictor,&incx);

    /* 
      w <- w + hvar * predictor.
    */
    daxpy_(&isix,&hvar,predictor,&incx,w,&incx);

    tim += dt;
    /*
      call deriv(w,dw,pot,dpotx,dpoty,dpotz,dmax)
      do 12 j=1,6
      array(6,j)=acst*dw(j)
      do 11 i=1,4
      array(i,j)=array(i+1,j)
   11 array(6,j)=array(i,j)*amcc(i)+array(6,j)
      array(5,j)=savdw(j)
   12 w(j)=savw(j)+hcvar*(array(5,j)+array(6,j))
      call deriv(w,dw,pot,dpotx,dpoty,dpotz,dmax)
    */

    mobcal_deriv(state,w,dw,&pot,&dmax);

    corrector = predictor;
    /*
      corrector <- dw
    */
    dcopy_(&isix,dw,&incx,corrector,&incx);
    /*
      Shift columns 1 to 4 of array into columns 0 to three.
      We cannot use dcopy for this as the source and target overlap.
      We use memmove instead.
    */
    move_size = (int64_t)(four_ld_array * sizeof(double));
    memmove(array,array_col_1,move_size);
    /*
      Compute the correcotr.
    */
    beta = acst;
    alpha = 1.0;

    dgemv_(&trans,&isix,&ifour,&alpha,array,&ld_array,amcc,&incx,
	   &beta,corrector,&incx);
    /*
      Move savdw into column 4 of array.
    */
    dcopy_(&isix,savdw,&incx,array_col_4,&incx);
    /*
      move saved w into w.
    */
    dcopy_(&isix,savw,&incx,w,&incx);
    /*
      Add savdw to corrector.
    */
    daxpy_(&isix,&hcvar,savdw,&incx,w,&incx);
    /*
      Correct w.
    */
    daxpy_(&isix,&hcvar,corrector,&incx,w,&incx);
    mobcal_deriv(state,w,dw,&pot,&dmax);
  } /* end else l was < 0 */
  *l_p    = l;
  *tim_p  = tim;
  *dt_p   = dt;
  *pot_p  = pot;
  *dmax_p = dmax;
}

 
