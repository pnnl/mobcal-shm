#include "system_includes.h"
#include "mobcal_state_struct.h"
#include "mobcal_get_filenames.h"
int mobcal_get_filenames(int argc, char **argv,
			 char *filenames, int filename_limit) {
  /*
    Get the names of the parameters file (parameters.in),
    atomtype parameters file (atomtype_parameters.in),
    mfj file (*.mfj), and output file.
    Called by: mobcal_read_parameters
    Calls:     fopen,fgets,strcpy,sscanf,fprintf,fflush,fclose
  */
  char *param_file;
  char *at_param_file;
  char *mfj_file;
  char *output_file;
  char *input_buffer;
  char *infile;
  char *ib;
  char input_b[2048];
  int success;
  int nr;
  int buff_len;
  int padi;
  FILE *ifp;
  FILE *efp;
  param_file = &filenames[0];
  at_param_file = &param_file[filename_limit];
  mfj_file      = &at_param_file[filename_limit];
  output_file   = &mfj_file[filename_limit];
  input_buffer  = (char*)&input_b[0];
  infile        = input_buffer;
  buff_len      = filename_limit;
  success       = 1;
  if (argc == 2) {
    strcpy(infile,(char*)argv[1]);
    ifp = fopen(infile,"r");
    if (ifp == NULL) {
      success = 0;
      fprintf(stderr, "mobcal_get_filenames: Error unable to open %s\n",infile);
      fflush(stderr);
    } 
    if (success) {
      ib = fgets(input_buffer,buff_len,ifp);
      if (ib != input_buffer) {
	success = 0;
	fprintf(stderr,"mobcal_get_filenames: Error reading parameter file name from %s\n",infile);
	fflush(stderr);
      }
    }
    if (success) {
      nr = sscanf(input_buffer,"%s",param_file);
      if (nr < 1) {
	success = 0;
	fprintf(stderr,"mobcal_get_filenames: Error scanning paremeter file name\n");
	fflush(stderr);
      }
    }
    if (success) {
      ib = fgets(input_buffer,buff_len,ifp);
      if (ib != input_buffer) {
	success = 0;
	fprintf(stderr,"mobcal_read_get_filenames: Error reading atomtypes parameter file name from %s\n",infile);
	fflush(stderr);
      }
    }
    if (success) {
      nr = sscanf(input_buffer,"%s",at_param_file);
      if (nr < 1) {
	success = 0;
	fprintf(stderr,"mobcal_get_filenames: Error scanning atomtype paremeter file name\n");
	fflush(stderr);
      }
    }
    if (success) {
      ib = fgets(input_buffer,buff_len,ifp);
      if (ib != input_buffer) {
	success = 0;
	fprintf(stderr,"mobcal_get_filenames: Error reading mfj file name from %s\n",infile);
	fflush(stderr);
      }
    }
    if (success) {
      nr = sscanf(input_buffer,"%s",mfj_file);
      if (nr < 1) {
	success = 0;
	fprintf(stderr,"mobcal_get_filenames: Error scanning mfj input file name\n");
	fflush(stderr);
      }
    }
    if (success) {
      ib = fgets(input_buffer,buff_len,ifp);
      if (ib != input_buffer) {
	success = 0;
	fprintf(stderr,"mobcal_get_filenames: Error reading output file name from %s\n",infile);
	fflush(stderr);
      }
    }
    if (success) {
      nr = sscanf(input_buffer,"%s",output_file);
      if (nr < 1) {
	success = 0;
	fprintf(stderr,"mobcal_get_filenames: Error scanning output file name\n");
	fflush(stderr);
      }
    }
    if (ifp) {
      fclose(ifp);
    }
  } else {
    strcpy(param_file,argv[1]);
    strcpy(at_param_file,argv[2]);
    /* mfj_file was formerly filen1. */
    strcpy(mfj_file,argv[3]);
    /* output_file was formerly filen2 */
    strcpy(output_file,argv[4]);
  }
  return(success);
}
