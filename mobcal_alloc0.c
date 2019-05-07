#include "system_includes.h"
#include "mobcal_state_struct.h"
#include "ranlux_state_struct.h"
#include "aligned_calloc.h"
#include "mobcal_alloc0.h"
int mobcal_alloc0(struct mobcal_state_struct **state_p, int filename_limit) {
  /*
    Allocate the state structure from the heap, and its
    character variables which are targets in the
    read_parameters routine.
    Allocate space for the random number generator state
    Allocate space for differential equation solver.
    ranlux_state_struct and its vector fields, seeds inext and ndskip,
    Allocate space for the atomtype_parameters table
    (fields gmass,geolj,grolj,grhs)
    Allocate space for the input line buffer (2048 bytes for now).
    Called by: mobcal, mobcal_shm_main
    Calls:     aligned_calloc, fprintf,fflush
  */
  struct mobcal_state_struct state_instance;
  struct mobcal_state_struct *state;
  struct ranlux_state_struct rstate_instance;
  struct ranlux_state_struct *rstate;
  struct mobcal_diffeq_state_struct diffeq_instance;
  struct mobcal_diffeq_state_struct *diffeq_state;
  struct mt_state_struct mt_state_instance;
  struct mt_state_struct *mt_state;
  double *seeds;
  double *gprops;
  /*
  double *gmass;
  double *geolj;
  double *grolj;
  double *grhs;
  */
  double *diffeq_a;
  double *array;
  double *predictor;
  char   *param_file;
  char   *at_param_file;
  char   *mfj_file;
  char   *output_file;
  /*
  char *filen1;
  char *filen2;
  */
  char *unit;
  char *dchar;
  char *xlabel;
  char *input_buffer;
  int64_t char_var_len;
  int64_t num_char_vars;
  int64_t one_l;
  int64_t ask_for;
  int64_t l_64;
  int    *inext;
  int success;
  int max_atom_imass;
  int input_buffer_len;
  int ld_array;
  int four_ld_array;
  int five_ld_array;
  success = 1;
  l_64    = (int64_t)64;
  one_l   = (int64_t)1;
  ask_for = (int64_t)sizeof(state_instance);
  state = (struct mobcal_state_struct *)aligned_calloc(one_l,ask_for,l_64);
  if (state == NULL) {
    success = 0;
    fprintf(stderr,
	    "mobcal_alloc0: Could not allocate %ld bytes for state\n",
	    ask_for);
    fflush(stderr);
  }
  if (success) {
    state->usage = ask_for;
    *state_p = state;
    char_var_len = (int64_t)filename_limit;
    num_char_vars = (int64_t)8;
    num_char_vars += (num_char_vars & one_l);
    ask_for = char_var_len * num_char_vars;
    state->usage  += ask_for;
    param_file  = (char*)aligned_calloc(one_l,ask_for,l_64);
    if (param_file == NULL) {
      success = 0;
      fprintf(stderr,"mobcal_alloc0: Error could not allocate %ld bytes for character variables\n",ask_for);
      fflush(stderr);
    }
  }
  if (success) {
    at_param_file = (char*)&param_file[char_var_len];
    mfj_file      = (char*)&at_param_file[char_var_len];
    output_file   = (char*)&mfj_file[char_var_len];
    unit          = (char*)&output_file[char_var_len];
    dchar         = (char*)&unit[char_var_len];
    xlabel        = (char*)&dchar[char_var_len];
    state->param_file = param_file;
    state->at_param_file = at_param_file;
    state->mfj_file = mfj_file;
    state->output_file = output_file;
    state->unit   = unit;
    state->dchar  = dchar;
    state->xlabel = xlabel;
  }
  if (success) {
    ask_for = (int64_t)sizeof(rstate_instance);
    state->usage  += ask_for;
    rstate  = (struct ranlux_state_struct *)aligned_calloc(one_l,ask_for,l_64);
    if (rstate == NULL) {
      success = 0;
      fprintf(stderr,
	     "mobcal_alloc0: Could not allocate %ld bytes for ranlux state\n",
	      ask_for);
      fflush(stderr);
    }
  }
  if (success) {
    ask_for = (int64_t)(24*sizeof(double));
    state->usage += ask_for;
    seeds = (double *)aligned_calloc(one_l,ask_for,l_64);
    if (seeds == NULL) {
      success = 0;
      fprintf(stderr,
	     "mobcal_alloc0: Could not allocate %ld bytes for ranlux state seeds field\n",
	      ask_for);
      fflush(stderr);
    }
  }
  if (success) {
    ask_for = (int64_t)(64*sizeof(int));
    state->usage += ask_for;
    inext   = (int *)aligned_calloc(one_l,ask_for,l_64);
    if (inext == NULL) {
      success = 0;
      fprintf(stderr,
	     "mobcal_alloc0: Could not allocate %ld bytes for ranlux state inext and ndskip fields\n",
	      ask_for);
      fflush(stderr);
    }
  }
  if (success) {
    state->ranlux_state   = rstate;
    rstate->seeds         = seeds;
    rstate->inext         = inext;
    rstate->ndskip        = (int*)&inext[24];
    rstate->isdext        = (int*)&inext[32];
    rstate->uninitialized = 1;
    ask_for = (int64_t)sizeof(mt_state_instance);
    state->usage  += ask_for;
    mt_state  = (struct mt_state_struct *)aligned_calloc(one_l,ask_for,l_64);
    if (mt_state == NULL) {
      success = 0;
      fprintf(stderr,
	     "mobcal_alloc0: Could not allocate %ld bytes for mt state\n",
	      ask_for);
      fflush(stderr);
    }
  }
  if (success) {
    state->mt_state = mt_state;
    max_atom_imass        = 256;
    ask_for = (((int64_t)4) * max_atom_imass) * sizeof(double);
    state->usage += ask_for;
    gprops = (double *)aligned_calloc(one_l,ask_for,l_64);
    if (gprops == NULL) {
      success = 0;
      fprintf(stderr,
	      "mobcal_alloc0: Could not allocate %ld bytes for atomtype_parameters gmass, geolj, grolj, and grhs\n",ask_for);
      fflush(stderr);
    }
  }
  if (success) {
    state->max_atom_imass = max_atom_imass;
    /*
    gmass  = gprops;
    geolj  = &gmass[max_atom_imass];
    grolj  = &geolj[max_atom_imass];
    grhs   = &grolj[max_atom_imass];
    */
    state->gprops = gprops;
    input_buffer_len = filename_limit;
    state->input_buffer_len = filename_limit;
    ask_for         = input_buffer_len;
    input_buffer     = (char*)aligned_calloc(one_l,ask_for,l_64);
    if (input_buffer == NULL) {
      success = 0;
      fprintf(stderr,
	      "mobcal_alloc0: Could not allocate %ld bytes for input_buffer\n",ask_for);
      fflush(stderr);
    }
  }
  if (success) {
    state->input_buffer = input_buffer;
    ask_for = (int64_t)sizeof(diffeq_instance);
    diffeq_state = (struct mobcal_diffeq_state_struct *)aligned_calloc(one_l,ask_for,l_64);
    if (diffeq_state == NULL) {
      success = 0;
      fprintf(stderr,
	      "mobcal_alloc0: Could not allocate %ld bytes for diffeq_state\n",ask_for);
      fflush(stderr);
    }
  }
  if (success) {
    state->diffeq_state = diffeq_state;
    /*
      Now allocate the fields of diffeq_state:
      double *a;       length 4 
      double *b;       length 4 
      double *c;       length 4 
      double *three_a  length 4
      double *ampc;    length 5 
      double *amcc;    length 4 
      double *savw;    length 6 
      double *savdw;   length 6 
      double *r        length 6
      double *q;       length 6 
      double *array;   length 36

      We have 11 entities. To keep things cache aligned, we will allocate
      them on 8 double word boundaries.
      So if we ask for 144 words of space that will be enough.
    */
    ask_for = ((int64_t)144) * sizeof(double);
    diffeq_a = (double *)aligned_calloc(one_l,ask_for,l_64);
    if (diffeq_a == NULL) {
      success = 0;
      fprintf(stderr,
	      "mobcal_alloc0: Could not allocate %ld bytes for diffeq_state fields a,b,c,ampc, amcc, q, savw, savdw, q_scale, dw_scale and array\n",ask_for);
      fflush(stderr);
    }
  }
  if (success) {
    diffeq_state->a       = diffeq_a;
    diffeq_state->b       = &diffeq_a[8];
    diffeq_state->c       = &diffeq_a[16];
    diffeq_state->three_a = &diffeq_a[24];
    diffeq_state->ampc    = &diffeq_a[32];
    diffeq_state->amcc    = &diffeq_a[40];
    diffeq_state->savw    = &diffeq_a[48];
    diffeq_state->savdw   = &diffeq_a[56];
    diffeq_state->q       = &diffeq_a[64];
    diffeq_state->r       = &diffeq_a[72];
    array                 = &diffeq_a[80];
    diffeq_state->array    = array;
    /*
      mobcall_diffeq needs a 6 by 6 array for array, 
      For aligment purposes we will give it an 8 by 8 array.
    */
    ld_array = 8;
    four_ld_array = 4 * ld_array;
    five_ld_array = four_ld_array + ld_array;
    predictor = (double *)&array[five_ld_array];
    diffeq_state->ld_array = ld_array;
    diffeq_state->four_ld_array = four_ld_array;
    diffeq_state->predictor = predictor;
    diffeq_state->corrector = predictor;
    diffeq_state->col_4     = (double *)&array[four_ld_array];
    diffeq_state->col_1     = (double *)&array[ld_array];
  }
  return(success);
}

