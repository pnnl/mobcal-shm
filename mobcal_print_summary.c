#include "system_includes.h"
#include "mobcal_state_struct.h"

#include "mobcal_print_summary.h"
void mobcal_print_summary(struct mobcal_state_struct *state, double sdevpc) {
  /* 
     Print a summary of the results.
     Called by: mobcal
     Calls:     fprintf, fflush, fclose
  */
  char *mfj_file;
  char *xlabel;
  double *asympp;
  double *ehsm;
  double *pam;
  double *pac;
  double *ehsc;
  double *tmm;
  double *tmc;
  double temp;
  double sw1;
  double sw2;
  double dtsf1;
  double dtsf2;
  double recip_pam;
  double recip_ehsm;
  double recip_tmm;
  double ten_p20;
  double scaled_pac;
  double scaled_ehsc;
  double scaled_tmc;
  double tmcs;
  double scaled_asymp;
  double pacs;
  double pamob;
  double ehscs;
  double ehsmob;
  double aasymp;
  double scaled_pacs;
  double scaled_ehscs;
  double scaled_aasymp;
  double scaled_tmcs;
  double recip_pamob;
  double recip_ehsmob;
  double recip_icoord;
  double tmmob;
  double recip_tmmob;

  int  charge_dist;
  int  icoord;

  int  i5;
  int  i2;

  int  imm;
  int  total_np;

  int  itn;
  int  inp;

  int  imp;
  int  immmin;

  int  immmax;
  int  ifailc;

  int  i;
  int  ip1;

  int  inum;
  

  FILE *ofp;
  FILE *efp;
  ofp         = state->ofp;
  mfj_file    = state->mfj_file;
  xlabel      = state->xlabel;
  charge_dist = state->charge_dist;
  icoord      = state->icoord;
  i5          = state->i5;
  i2          = state->i2;
  asympp      = state->asympp;
  ehsm        = state->ehsm;
  ehsc        = state->ehsc;
  pam         = state->pam;
  pac         = state->pac;
  tmm         = state->tmm;
  temp        = state->temp;
  imm         = state->imm;
  inum        = state->inum;
  sw1         = state->sw1;
  sw2         = state->sw2;
  dtsf1       = state->dtsf1;
  dtsf2       = state->dtsf2;
  itn         = state->itn;
  inp         = state->inp;
  imp         = state->imp;
  ifailc      = state->ifailc;
  immmin      = state->immmin;
  immmax      = state->immmax;
  tmc         = state->tmc;
  ten_p20     = 1.e20;
  if (ofp) {
    fprintf(ofp,"\n\nSUMMARY\n\n program version = junkn.f\n"
	    " input file name = %s\n input file label = %s\n",
	    mfj_file,xlabel);
    if (charge_dist == 1) {
      fprintf(ofp,"using a uniform charge distribution\n");
    } else {
      if (charge_dist == 2) {
	fprintf(ofp,"using a calculated (non-uniform) charge distribution\n");
      } else {
	fprintf(ofp,"using no charge - only LJ interactions\n");
      }
    }
    fprintf(ofp,"temperature = %le\n",temp);
    if (i5 != 1) {
      fprintf(ofp,"using RAND with seed integer = %d\n",i2);
    } else {
      fprintf(ofp,"using RANLUX with seed integer = %d\n",i2);
    }
    if (icoord == 1) {
      fprintf(ofp,"structural asymmetry parameter = %le\n\n",
	      asympp[0]);
      if (ehsm[0] != 0.0) {
	fprintf(ofp,"mobility calculation by MOBIL4 (HS scattering)\n\n");
	fprintf(ofp,"number of Monte Carlo trajectories = %d\n"
		"maximun number of reflections encounter = %d\n",inum,
		imm);
	recip_pam = 1.0/pam[0];
	scaled_pac = pac[0] * ten_p20;
	fprintf(ofp,"inverse average PA mobility = %le\n",recip_pam);
	fprintf(ofp,"average PA cross section = %le\n",scaled_pac);
	recip_ehsm = 1.0/ehsm[0];
	scaled_ehsc = ehsc[0] * ten_p20;
	fprintf(ofp,"inverse average EHS mobility = %le\n",recip_ehsm);
	fprintf(ofp,"average EHS cross section = %le\n",scaled_ehsc);
      }
      if (tmm[0] != 0.0) {
	fprintf(ofp,"mobility calculation by MOBIL2 (trajectory method)\n\n");
	fprintf(ofp,"tractory parameters\n");
	fprintf(ofp,"sw1 = %le       sw2 = %le\n",sw1,sw2);
	fprintf(ofp,"dtsf1 = %le     dtsf2 = %le\n\n",dtsf1,dtsf2);
	total_np = itn * inp * imp;
	fprintf(ofp,"number of complete cycles (itn) = %d\n",itn);
	fprintf(ofp,"number of velocity points (inp) = %d\n",inp);
	fprintf(ofp,"number of random points (imp) = %d\n",imp);
	fprintf(ofp,"total number of points = %d\n\n",total_np);
	recip_tmm = 1.0/tmm[0];
	scaled_tmc = tmc[0] * ten_p20;
	fprintf(ofp,"inverse average (second order) TM mobility = %le\n",
		recip_tmm);
	fprintf(ofp,"average TM cross section = %le\n",scaled_tmc);
	fprintf(ofp,"standard deviation (percent) = %le\n",sdevpc);
	fprintf(ofp,"number of failed tracjectories = %d\n",ifailc);
	fflush(ofp);
      }
      /* iend if (icoord == 1) */ 
    } else { 
      /*
	icoord > 1
      */
      recip_icoord = 1.0/((double)icoord);
      if (ehsm[0] != 0.0) {
	fprintf(ofp,"mobility calculation by MOBIL4 (HS scattering)\n");
	fprintf(ofp,"number of Monte Carlo trajectories = %d\n",inum);
	fprintf(ofp,"minimum and maximum number of reflections = %d  %d\n",
		immmin,immmax);
	fprintf(ofp,"\n\n    set     PA CS       PA MOB^-1      EHSS CS      EHSS MOB^-1    ASYMP\n");
	for (i=0;i<icoord;i++) {
	  ip1 = i + 1;
	  scaled_pac = pac[i] * ten_p20;
	  recip_pam  = 1.0 / pam[i];
	  scaled_ehsc = ehsc[i]*ten_p20;
	  recip_ehsm  = 1.0 / ehsm[i];
	  scaled_asymp = asympp[i] * 0.1;
	  fprintf(ofp,"%d %le %le %le %le %le\n",
		  ip1,scaled_pac,recip_pam,
		  scaled_ehsc,recip_ehsm,scaled_asymp);
	}
	pacs   = 0.0;
	pamob  = 0.0;
	ehscs  = 0.0;
	ehsmob = 0.0;
	aasymp = 0.0;
	for (i=0;i<icoord;i++) {
	  pacs += pac[i];
	  pamob += pam[i];
	  ehscs += ehsc[i];
	  ehsmob += ehsm[i];
	  aasymp += asympp[i];
	}
	pacs   = pacs * recip_icoord;
	pamob  = pamob * recip_icoord;
        ehscs  = ehscs * recip_icoord;
	ehsmob = ehsmob * recip_icoord;
	aasymp = aasymp * recip_icoord;
	scaled_pacs = pacs * ten_p20;
	recip_pamob = 1.0 / pamob;
	scaled_ehscs = ehscs * ten_p20;
	recip_ehsmob = 1.0 / ehsmob;
	scaled_aasymp = aasymp * 0.1;
	fprintf(ofp,"AVGE %le %le %le %le %le\n",
		scaled_pacs,recip_pamob,
		scaled_ehscs,recip_ehsmob,scaled_aasymp);
      } /* end if (ehsm[0] != 0) */
      if (tmm[0] != 0) {
	fprintf(ofp,"mobility calculation by MOBIL2 (trajectory method)\n");
	fprintf(ofp,"tractory parameters\n");
	fprintf(ofp,"sw1 = %le       sw2 = %le\n",sw1,sw2);
	fprintf(ofp,"dtsf1 = %le     dtsf2 = %le\n",dtsf1,dtsf2);
	total_np = itn * inp * imp;
	fprintf(ofp,"number of complete cycles (itn) = %d\n",itn);
	fprintf(ofp,"number of velocity points (inp) = %d\n",inp);
	fprintf(ofp,"number of random points (imp) = %d\n",imp);
	fprintf(ofp,"total number of points = %d\n",total_np);
	fprintf(ofp,"number of failed tracjectories = %d\n",ifailc);
	fprintf(ofp,"\n\n    set     TM CS       TM MOB^-1\n");
	tmcs = 0.0;
	tmmob = 0.0;
	for (i=0;i<icoord;i++) {
	  ip1 = i+1;
	  scaled_tmc = tmc[i] * ten_p20;
	  tmcs       += tmc[i];
	  recip_tmm  = 1.0/tmm[i];
	  tmmob      += tmm[i];
	  fprintf(ofp,"%d   %le   %le\n",ip1,scaled_tmc,recip_tmm);
	}
	tmcs = tmcs * recip_icoord;
	tmmob = tmmob * recip_icoord;
	scaled_tmcs = tmcs * ten_p20;
	recip_tmmob = 1.0/tmmob;
	fprintf(ofp,"AVGE  %le   %le\n",scaled_tmcs,recip_tmmob);
      }
    }
    fclose(ofp);
  }
}
