#include "system_includes.h"
#include "mobcal_state_struct.h"
#include "mobcal_get_filenames.h"
#include "mobcal_read_parameters.h"
int mobcal_read_parameters(int argc, char **argv,
			   char *filenames, int *iparams, double *dparams) {
  /*
    Read in the parameter, atomtype_parameter, input_mfj, and output  filenames,    storing them in filenames, and the 
    i2,ipr,itn,inp,imp,irn, ibstmax,
    igs, im2, im4, ip, it, iu1, iu2, iu3, iv, use_dgt, buffer_gas, use_mt,
    parameters storing them in that order in iparams, 
    Also store filename_limit in iparams after use_mt.
    Store the temp and m1 (BUFFER_GAS_MASS) parameters in dparams.

    Called by: mobcal, mobcal_shm
    Calls:     fopen, fgets, sscanf, strcmp, fprintf, fflush, atoi
  */
  double temp;
  double m1;
  double bg_dipole_mult;
  char *infile;
  char *param_file;
  char *at_param_file;
  char *mfj_file;
  char *output_file;
  char *input_buffer;
  char *ib;
  char *filen1;
  char *filen2;
  char input_b[2048];
  char *key;
  char key_b[32];
  char *value;
  char value_b[32];
  int64_t key_upcase;
  int64_t *ikey;
  int  success;
  int  i2;

  int  nr;
  int  buff_len;

  int  ipr;
  int  itn;

  int  inp;
  int  imp;

  int  irn;
  int  ibstmax;

  int igs;
  int im2;

  int im4;
  int ip;

  int it;
  int iu1;

  int iu2;
  int iu3;

  int iv;
  int use_dgt;

  int buffer_gas;
  int i;

  int ipos;
  int use_mt;

  int num_threads;
  int filename_limit;

  int use_iu_dljpot;
  int avec_len;

  FILE *ifp;
  FILE *pfp;
  FILE *efp;
  /*
    The hardwired name of the input file mobcal.in could be
    taken from the first agument.
  */
  success = 1;
  buff_len = 64;
  filename_limit = 2048;
  param_file = &filenames[0];
  at_param_file = &param_file[filename_limit];
  mfj_file      = &at_param_file[filename_limit];
  output_file   = &mfj_file[filename_limit];
  input_buffer  = (char*)&input_b[0];
  infile        = input_buffer;
  avec_len      = _AVEC_LEN_;
  if ((argc != 2) && (argc !=5))  {
    fprintf(stdout,
	    " mobcal usage is\n"
    " mobcal paremeters.in atomtype_parameters.in input.mfj output_file\n"
	    " or\n"
	    " mobcal infile \n"
	    " where infile is the name of a file that has those four filenames on succesive lines.\n"
	    );
    fflush(stdout);
    success = 0;
  }
  if (success) {
    success = mobcal_get_filenames(argc,argv,filenames, filename_limit);
  }
  if (success) {
    pfp = fopen(param_file,"r");
    if (pfp == NULL) {
      success = 0;
      fprintf(stderr,"mobcal_read_parameaters: Error opening parameter file %s\n",param_file);
      fflush(stderr);
    }
  }
  if (success) {
    /*
      Default random number seed.
    */
    i2 = 5013489;
    /*
      Define parameters for POTENT
      Number of rotations in average potential calculation. If ipr=0
      an unrotated configuration is used. Otherwise ipr random rotations
      are employed.
    */
    ipr = 1000;
    /*
      Define parameters for MOBIL2
      itn = Number of complete cycles for average mobility calculation in 
      MOBIL2. Default value is 10.
    */
    itn = 1;
    /*
      Inp = Number of points in velocity integration in MOBIL2. Default 
      value is 20.
    */
    inp=40;
    /*
      Number of points in Monte Carlo integrations of impact parameter
      and orientation in MOBIL2. Default value is 500.
    */
    imp=1000;
    /*
      max length of cosx vector int mobcal_mobil2
    */
    ibstmax = 500;
    /*
      irn = reciprocal of step size weh computing lj potentials.
    */
    irn = 1000;
    temp = 301.0;
    /*
      Allow user to keyword specify these fields at the end of the 
      paremeter input file.
    */
    igs = 0;
    im2 = 0;
    im4 = 0;
    ip = 0;
    it = 0;
    iu1 = 0;
    iu2 = 0;
    iu3 = 0;
    iv = 0;
    num_threads = 16;
    use_dgt = 1;
    use_iu_dljpot = 0;
    buffer_gas = 14;
    bg_dipole_mult = 1.710e-30;
    use_mt  = 1;
    key = (char*)&key_b[0];
    value = (char*)&value_b[0];
    ikey = (int64_t *)key;
    key_upcase = 0x5f5f5f5f5f5f5f5f;
    ib = fgets(input_buffer,buff_len,pfp);
    while (!feof(pfp)) {
      nr = sscanf(input_buffer,"%s %s",key,value);
      
      if (nr == 2) {
	/*
	  Only parse lines with two "strings"
	*/
	/*
	  upper case the key value.
	*/
	for (i=0;i<strlen(key);i++) { 
	  if ((key[i] >= 97) && key[i] <= 122) {
	    key[i] = key[i]-32;
	  }
	}
	if (strcmp(key,"I2") == 0) {
	  i2 = atoi(value);
	} else if (strcmp(key,"IBSTMAX") == 0) {
	  ibstmax = atoi(value);
	} else if (strcmp(key,"IGS") == 0) {
	  igs = atoi(value);
	} else if (strcmp(key,"IMP") == 0) {
	  imp = atoi(value);
	} else if (strcmp(key,"IM2") == 0) {
	  im2 = atoi(value);
	} else if (strcmp(key,"IM4") == 0) {
	  im4 = atoi(value);
	} else if (strcmp(key,"INP") == 0) {
	  inp = atoi(value);
	} else if (strcmp(key,"IPR")  == 0) {
	  ipr = atoi(value);
	} else if (strcmp(key,"IP") == 0) {
	  ip = atoi(value);
	} else if (strcmp(key,"ITN") == 0) {
	  itn = atoi(value);
	} else if (strcmp(key,"IRN") == 0) {
	  irn = atoi(value);
	} else if (strcmp(key,"IT") == 0) {
	  it = atoi(value);
	} else if (strcmp(key,"IU1") == 0) {
	  iu1 = atoi(value);
	} else if (strcmp(key,"IU2") == 0) {
	  iu2 = atoi(value);
	} else if (strcmp(key,"IU3") == 0) {
	  iu3 = atoi(value);
	} else if (strcmp(key,"IV") == 0) {
	  iv = atoi(value);
	} else if (strcmp(key,"TEMP") == 0) {
	  nr = sscanf(value,"%le",&temp);
	} else if (strcmp(key,"USE_DGT") == 0) {
	  use_dgt = atoi(value);
	} else if (strcmp(key,"USE_MT") == 0) {
	  use_mt = atoi(value);
	} else if (strcmp(key,"BUFFER_GAS") == 0) {
	  if (strcmp(value,"NITROGEN") == 0) {
	    buffer_gas = 14;
	    use_iu_dljpot = 0;
	  } else if (strcmp(value,"HELIUM") == 0) {
	    buffer_gas = 4;
	    use_iu_dljpot = 1;
	  } else {
	    buffer_gas = atoi(value);
	    if (buffer_gas == 4) {
	      use_iu_dljpot = 1;
	    } else {
	      use_iu_dljpot = 0;
	    }
	  }
	} else if (strcmp(key,"BUFFER_GAS_MASS") == 0) {
	  nr = sscanf(value,"%le",&m1);
	} else if (strcmp(key,"BG_DIPOLE_MULT") == 0) {
	  nr = sscanf(value,"%le",&bg_dipole_mult);
	} else if (strcmp(key,"NUM_THREADS") == 0) {
	  num_threads = atoi(value);
	  if (num_threads < 2) {
	    num_threads = 2;
	  }
	} else if (strcmp(key,"FILENAME_LIMIT") == 0) {
	  filename_limit = atoi(value);
	  if (filename_limit < 256) {
	    filename_limit = 256;
	  }
	  /* We do not let the user specifiy this any more.
             if buffer gas is helium, 
	} else if (strcmp(key,"USE_IU_DLJPOT") == 0) {
	  use_iu_dljpot = atoi(value);
	  */
	} else if (strcmp(key,"AVEC_LEN") == 0) {
	  avec_len = atoi(value);
	}
      } else {
	fprintf(stderr,
		"mobcal_read_paremeters ignored the followint line\n%s\n",
		input_buffer);
	fflush(stderr);
      }
      ib = fgets(input_buffer,buff_len,pfp);
    } /* end while (!eof(pfp)) */
    if (pfp) {
      fclose(pfp);
    }
  }
  if (success) {
    /*
      Slots 0 through 7 in iparams are used for 
      num_procs,  my_rank,
      numthreads, shm_block_size, shm_size, shmmid, and shm_projid,
      and mythread_id.
    */
    ipos = 8;
    iparams[ipos] = i2;
    ipos += 1;
    iparams[ipos] = ipr;
    ipos += 1;
    iparams[ipos] = itn;
    ipos += 1;
    iparams[ipos] = inp;
    ipos += 1;
    iparams[ipos] = imp;
    ipos += 1;
    iparams[ipos] = irn;
    ipos += 1; 
    iparams[ipos] = ibstmax;
    ipos += 1; 
    iparams[ipos] = igs;
    ipos += 1;
    iparams[ipos] = im2;
    ipos += 1;
    iparams[ipos] = im4;
    ipos += 1;
    iparams[ipos] = ip;
    ipos += 1;
    iparams[ipos] = it;
    ipos += 1;
    iparams[ipos] = iu1;
    ipos += 1;
    iparams[ipos] = iu2;
    ipos += 1;
    iparams[ipos] = iu3;
    ipos += 1;
    iparams[ipos] = iv;
    ipos += 1;
    iparams[ipos] = use_dgt;
    ipos += 1;
    /*
#define DBG 1
    */
#ifdef DBG
    fprintf(stdout,"mobcal_read_params: ipos = %d, buffer_gas = %d\n",ipos,buffer_gas);
    fflush(stdout);
#endif    
    iparams[ipos] = buffer_gas;
    ipos += 1;
    iparams[ipos] = use_mt;
    ipos += 1;
    iparams[ipos] = filename_limit;
    ipos += 1;
    iparams[ipos] = use_iu_dljpot;
    ipos += 1;
    iparams[ipos] = avec_len;
    iparams[2]    = num_threads;
    dparams[0]  = temp;
    dparams[1]  = m1;
    dparams[2]  = bg_dipole_mult;
  }
  return(success);
}
