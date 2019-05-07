#include "system_includes.h"
#include "mobcal_state_sruct.h"
#include "mobcal_ranlux.h"
int ranlux(struct mobcal_state_struct *state,double *rvec; int lenv) {
  /*
    Adapted from the fortran 77 flavor of the RANLUX code by F.James. 
    documentation from which is included below. Adaptaion by Doug Baxter 2017.
  */
  /*
C         Subtract-and-borrow random number generator proposed by
C         Marsaglia and Zaman, implemented by F. James with the name
C         RCARRY in 1991, and later improved by Martin Luescher
C         in 1993 to produce "Luxury Pseudorandom Numbers".
C     Fortran 77 coded by F. James, 1993
C          
C       references:
C  M. Luscher, Computer Physics Communications  79 (1994) 100
C  F. James, Computer Physics Communications 79 (1994) 111
C
C   LUXURY LEVELS.
C   ------ ------      The available luxury levels are:
C
C  level 0  (p=24): equivalent to the original RCARRY of Marsaglia
C           and Zaman, very long period, but fails many tests.
C  level 1  (p=48): considerable improvement in quality over level 0,
C           now passes the gap test, but still fails spectral test.
C  level 2  (p=97): passes all known tests, but theoretically still
C           defective.
C  level 3  (p=223): DEFAULT VALUE.  Any theoretically possible
C           correlations have very small chance of being observed.
C  level 4  (p=389): highest possible luxury, all 24 bits chaotic.
C
C!!! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C!!!  Calling sequences for RANLUX:                                  ++
C!!!      CALL RANLUX (RVEC, LEN)   returns a vector RVEC of LEN     ++
C!!!                   32-bit random floating point numbers between  ++
C!!!                   zero (not included) and one (also not incl.). ++
C!!!      CALL RLUXGO(LUX,INT,K1,K2) initializes the generator from  ++
C!!!               one 32-bit integer INT and sets Luxury Level LUX  ++
C!!!               which is integer between zero and MAXLEV, or if   ++
C!!!               LUX .GT. 24, it sets p=LUX directly.  K1 and K2   ++
C!!!               should be set to zero unless restarting at a break++
C!!!               point given by output of RLUXAT (see RLUXAT).     ++
C!!!      CALL RLUXAT(LUX,INT,K1,K2) gets the values of four integers++
C!!!               which can be used to restart the RANLUX generator ++
C!!!               at the current point by calling RLUXGO.  K1 and K2++
C!!!               specify how many numbers were generated since the ++
C!!!               initialization with LUX and INT.  The restarting  ++
C!!!               skips over  K1+K2*E9   numbers, so it can be long.++
C!!!   A more efficient but less convenient way of restarting is by: ++
C!!!      CALL RLUXIN(ISVEC)    restarts the generator from vector   ++
C!!!                   ISVEC of 25 32-bit integers (see RLUXUT)      ++
C!!!      CALL RLUXUT(ISVEC)    outputs the current values of the 25 ++
C!!!                 32-bit integer seeds, to be used for restarting ++
C!!!      ISVEC must be dimensioned 25 in the calling program        ++
C!!! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  */
  /*
    N.B. ISVEC most likely refers to isdext in the code.
  */
  struct ranlux_state_struct *rstate;
  double uni;
  double carry;
  double twom24;
  double twom48;
  double twom12;
  double *seeds;
  int *inext;

  int success;
  int lux;

  int jseed;
  int ivec;

  int k1;
  int k2;

  int i24;
  int j24;

  int in24;
  int nskip;

  int kount;
  int mkount;

  int isk;
  
  success = 1;
  if (state == NULL) {
    success = 0;
  } 
  if (success) {
    rstate = state->ranlux_state;
    if (rstate == NULL) {
      success = 0;
    }
  }
  if (success) {
    if (rstate->uninitialized) {
      lux = 3;
      jseed = 314159265;
      k1  = 0;
      k2  = 0;
      success = rluxgo(state,lux.jseed,k1,k2);
    }
  }
  /*
C
C          The Generator proper: "Subtract-with-borrow",
C          as proposed by Marsaglia and Zaman,
C          Florida State University, March, 1989
C
  */
  if (success) {
    twom48  = rstate->twom48;
    twom24  = rstate->twom24;
    twom12  = rstate->twom12;
    seeds   = rstate->seeds;
    inext   = rstate->inext;
    i24     = rstate->i24;
    j24     = rstate->j24;
    in24    = rstate->in24;
    kount   = rstate->kount;
    mkount  = rstate->mkount;
    nskip   = rstate->nskip;
    carry   = rstate->carry;
    for (ivec=0;ivec < lenv; ivec+=1) {
      uni = seeds[j24] - seeds[i24] - carry;
      if (uni < 0.0) {
	uni += 1.0;
	carry = twom24;
      } else {
	carry = 0.0;
      }
      seeds[i24] = uni;
      i24 = inext[i24];
      j24 = jnext[j24];
      kcount += 1;
      if (uni < twom12) {
	/*
	  small numbers (with less than 12 "significant" bits) are "padded".
	*/
	uni += twom24 * seeds[j24];
	/*
	  And sero is forbidden in case some takes a log.
	*/
	if (uni == 0.0) {
	  uni = twom48;
	}
      }
      recv[ivec] = uni;
      in24 += 1;
      if (in24 == 24) {
	in24 = 0;
	/*
	  The following 5 lines are a bug fix, if 
	  if lenv * nskip > 2 * igiga you could 
	  wrap kcount.
	*/
	kount = kount + nskip;
	if (kount >= igiga) {
	  mkount = mkount + 1;
	  kount = kount - igiga;
	}
	for (isk = 0; isk < nskip; isk += 1) {
	  uni = seeds[j24] - seeds[i24] - carry;
	  if (uni < 0.) {
	    uni += 1.0;
	    carry = twom24;
	  } else {
	    carry = 0;
	  }
	  seeds[i24] = uni;
	  i24 = next[i24];
	  j24 = next[j24];
	}
      } /* end for (isk... ) */
    } /* end for (ivec ...) */
    /*
      The following five lines are no longer added as we
      added a 1 in the inside increment of kcount.
    kount = kount + lenv;
    if (kount >= igiga) {
      mkount = mkount + 1;
      kount = kount -igiga;
    }
    */
    rstate->carry  = carry;
    rstate->i24    = i24;
    rstate->j24    = j24;
    rstate->in24   = in24;
    rstate->kount  = kount;
    rstate->mkount = mkount;
  }
  return(success);
}
