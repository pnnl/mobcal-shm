void mobcal_che(struct mobcal_state_struct *state, int im, 
		double *cof_p, double *cop_p, double *yr_p,
		double *zr_p, int *kp_p) {

  /*
    Guides hard sphere scattering trajectory. Adapted from code 
    written by Alexandre Shvartsburg.
    Called by: mobil4
    Calls:     sqrt
  */
  double cof;
  double cop;
  double yr;
  double zr;
  double xl;
  double *rhs;
  double *rhs2;
  double *fx;
  double *fy;
  double *fz;
  double ras;
  double dev;
  double eps;
  double xfx;

  int kp;
  int ki;
  int in;
  int inatom;

  
  inatom = state->inatom;
  rhs    = state->rhs;
  rhs2   = state->rhs2;
  fx     = state->fx;
  fy     = state->fy;
  fz     = state->fz;
  kp = 0;
  eps = 1.0e-16;
  if (im != 1) {
    yr = 0.0;
    zr = 0.0;
  }
  xl = 1000000.0;
  ki = -1;
  for (in=0;in<inatom;in++) {
    xfx = fx[iatom];
    if ((xfx > eps) || (im == 1)) {
      yd = yr - fy[iatom];
      zd = zr - fz[iatom];
      ras = (yd*yd) + (zd * zd);
      dev = sqrt(ras);
      if (dev <= rhs[in]) {
	xc = fx- sqrt(rhs2[in] - ras);
	if (xc < xl) {
	  xl = xc;
	  ki = in;
	}
      }
    }
  }
  if (ki == -1) {
    cof = cop;
  } else {
    kp = 1;
  }
  *cof_p = cof;
  *cop_p = cop;
  *yr_p  = yr;
  *zr_p  = zr;
  *kp_p  = kp;
}
