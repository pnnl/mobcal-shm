#include "system_includes.h"
#include "mobcal_state_struct.h"
#include "mobcal_init_rn.h"
#include "mobcal_gen_rn.h"
#include "mobcal_init_constants_1.h"
int mobcal_init_constants_1 (struct mobcal_state_struct *state) {
  /*
    Initialize more physics and math constants, and the random
    number generators.
    Called by: mobcal, mobcal_shm_main
    Calls:     fprintf,fflush mobcal_init_rn, mobcal_gen_rn
    Uses: pi, xe, xn, xeo, m2, xmv, i2, ofp
    Sets  e0, ro, ro2, dipol, m1, mu, temp (t in fortran code),
          i1, i3, i4, i5
          
  */
  double *ranlist;
  double pi; 
  double xe;
  double xn;
  double xmv;
  double m2;
  double eo;
  double xeo;
  double ro;
  double ro2;
  double m1;
  double mu;
  double xk;
  double recip_mconst;
  double bg_dipole_mult;
  
  double unscaled_ro;
  double unscaled_eo;
  double dipol;
  double dens;
  double mconst;
  double t;

  int    i1;
  int    i2;

  int    i3;
  int    i4;

  int    i5;
  int    success;

  int    icoord;
  int    itn;

  int    inp;
  int    imp;

  int    ranlist_size;
  int    buffer_gas;

  int    use_mt;
  int    padi;

  FILE *ofp;
  FILE *lfp;

  success = 1;
  pi  = state->pi;
  xe  = state->xe;
  xn  = state->xn;
  xeo = state->xeo;
  m2  = state->m2;
  xmv = state->xmv;
  ofp = state->ofp;
  xk  = state->xk;
  icoord = state->icoord;
  itn    = state->itn;
  inp    = state->inp;
  imp    = state->imp;
  buffer_gas = state->buffer_gas;
  ranlist      = state->ranlist;
  ranlist_size = state->ranlist_size;
  /*
    Lennard-Jones scaling parameters
  */
  unscaled_eo = .00134;
  eo          = unscaled_eo*xe;
  state->eo   = eo;
  /*
    ro          = 3.043*1.0e-10;
  */
  unscaled_ro = 3.043;
  ro          = unscaled_ro * 1.0e-10;
  state->ro   = ro;
  ro2         = ro * ro;
  state->ro2  = ro2;
  state->pi_ro2 = pi * ro2;
  /*
    Constant for ion-induced dipole potential
     Polarizability of helium = 0.204956d-30 m3
     xeo is permitivity of vacuum, 8.854187817d-12 F.m-1
     dipol=1.641d-30/(2.d0*4.d0*pi*xeo)
  */
  dipol = 0.0;
  /*
#define DBG 1
  */
#ifdef DBG
  fprintf(stdout,"mobcal_init_constants_1: %d, buffer_gas = %d\n",
	  state->thread_id,buffer_gas);
  fflush(stdout);
#endif
  if (buffer_gas == 14) {
    /*
      Buffer gas is nitrogen
      Mass constants
      Original mobcal from iu had a different mass constant.
    */
    m1=28.0134;
    /*
      mobcal orignal source code had a different constant,
    */
    /*
      New value from DGT, but not used in his versions.
    bg_dipole_mult = 1.641e-30;
    */
    bg_dipole_mult = 1.710e-30;
  } else {
    if (buffer_gas == 4) {
      /*
	Buffer gas is helium.
      */
      m1 = 4.0026;
      bg_dipole_mult = 0.204956e-30;
    } else {
      m1 = state->m1;
      bg_dipole_mult = state->bg_dipole_mult;
      if (ofp) {
	fprintf(ofp,"mobcal_init_constants_1: error buffer_gas = %d, must be 14 for nitrogen or 4 for helium, using %le for m1\n",buffer_gas,m1);
	fflush(ofp);
      }
    }
  }
  dipol=bg_dipole_mult/(8.0*pi*xeo);
  state->m1             = m1;
  state->bg_dipole_mult = bg_dipole_mult;
  dipol                 = dipol*xe*xe;
  /*
    This value of dipol is now used int mobcal_dljpot and mobcal_dljpot_only
    for dipolxx, and dipolzz.
  */
  state->dipol          = dipol;
  /*
    mu=((m1*m2)/(m1+m2))/(xn*1.0e3);
  */
  mu=(m1*m2)/((m1+m2)*(xn*1000.0));
  state->mu = mu;
  state->recip_half_mu = 2.0/mu;
  /*
    Mobility constant
    
    mconst=sqrt(18.0*pi)/16.0
    = 3/16 * sqrt(2*pi)
  */
  mconst=0.1875*sqrt(pi+pi);
  /*
    mconst=mconst*dsqrt(xn*1.0d3)*dsqrt((1.d0/m1)+(1.d0/m2))
  */
  mconst = mconst * sqrt(1.0/mu);
  mconst = mconst*xe/sqrt(xk);
  dens   = xn/xmv;
  mconst = mconst/dens;
  state->mconst = mconst;
  recip_mconst  = 1.0/mconst;
  state->recip_mconst = recip_mconst;
  /*
    temp field of state now set in mobcal_read_pareameters.
    t = 301.0;
    state->temp = t;
  */
  /*
    Define parameters for random number generator

    If i5=1 RANLUX is used otherwise RAND is used. If RAND is used 
    i2 is the seed integer. If RANLUX is used i2, i3, and i4 are seed
    integers. i3 and i4 (which are used to start RANLUX in a particular
    place) are usually set to zero. i1 contains the "luxury level" - how
    good the random number generator is. Values of 0 to 4 can be used
    and the default value in RANLUX is 3. At this level "any 
    theoretically possible correlations have a very small chance of
    being observed".
    if use_mt is set to 1 (the default) the mersenned-twister mt19937-64
    generator is used with seed i2.
    
  */
  i2        = state->i2;
  i1        = 3;
  i3        = 0;
  i4        = 0;
  state->i1 = 3;
  state->i3 = 0;
  state->i4 = 0;
  i5        = 1;
  use_mt    = state->use_mt;
  state->i5 = i5;
  /*
    Initialze random number generators.
  */
  success = mobcal_init_rn(state);
  if (ofp) {
    fprintf(ofp,"Lennard-Jones scaling parameters: eo= %le ro=%le\n",
	    unscaled_eo,unscaled_ro);
    fprintf(ofp,"dipole constant =%le\n",dipol);
    fprintf(ofp,"mobility constant = %le\n",mconst);
    fprintf(ofp,"temperature = %le\n\n",t);
    if (i5 != 1) {
      fprintf(ofp,"using RAND with seed integer = %d\n",state->i2);
    } else {
      if (use_mt) {
	fprintf(ofp,"using mt19937-64 with seed = %d\n",state->i2);
      } else {
	fprintf(ofp,"using RANLUX with seed = %d\n",state->i2);
      }
    }
  }
  /*
    Here we will build all the random integers to be used, storing them
    in the ranlist field, actually we build 4 more than we need.
  */
  success = mobcal_gen_rn(state,ranlist,ranlist_size);
  /*
    Minimum value of (1-cosX). This quantity determines the maximum
    impact parameter at each velocity. Default value is 0.0005.
  */
  state->cmin=0.0005;
  /*
    Define some parameters for trajectories: sw1 defines the potential
    energy where the trajectory starts and dtsf1 is related to the 
    time step at the start of the trajectory. As the trajectory comes
    close to a collision the time step is reduced. sw2 defines the 
    potential energy where the time step is reduced and dtsf2 defines
    the reduced time step. Default values are: 
    sw1 = 0.00005   dtsf1=0.5
    sw2 = 0.0025    dtsf2=0.05
    inwr is the number of integration steps before the program tests
    to see if the trajectory is done or lost. ifail is the number of
    failed trajectories that are permitted (a failed trajectory is one
    that does not conserve energy to within 1%. Default values are: 
    inwr = 1        ifail = 100
  */
  state->sw1=0.00005;
  state->sw2=0.005;
  state->dtsf1=0.5;
  state->dtsf2=0.1;
  state->inwr=1;
  state->ifail=100;
  state->ifailc=0;
  /*
    Define parameters for MOBIL4
    inum is the number of Monte Carlo trajectories employed and inor
    is the maximum number of successive reflections followed. Default
    values inum=600000 and inor=30
  */
  state->inum=250000;
  state->inor=30;
  /*
    Define parameters for TRAJONE
  */
  state->v=2309.9;
  state->b=1.5654e-10;
  state->ntheta=33.300;
  state->nphi=64.250;
  state->ngamma=134.30;
  return(success);
}
