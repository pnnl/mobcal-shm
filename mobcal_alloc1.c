#include "system_includes.h"
#include "mobcal_state_struct.h"
#include "aligned_calloc.h"
#include "mobcal_alloc1.h"
int mobcal_alloc1(struct mobcal_state_struct *state) {
  /*
    Allocate the vector and character fields in state
    Requires the num_data_points and num_parameters fields
    of state to be set and positive.
    Called by: mobcal_ncoord
    Calls:     fprint,fflush, calloc
  */
  double *ox;
  double *oy;
  double *oz;
  double *fx;
  double *fy;
  double *fz;
  double *imass; 
  double *xmass;
  double *pcharge;
  double *eolj;
  double *rolj;
  double *eox4;
  double *ro6lj;
  double *ro12lj;
  double *dro6;
  double *dro12;
  double *rhs;
  double *rhs2;
  double *ranlist;
  double *dljpot_input_stream;
  double *asym_work;
  
  double *tmc;
  double *tmm;
  double *ehsc;
  double *ehsm;
  double *pac;
  double *pam;
  double *asympp;

  double *pgst;
  double *wgst;
  double *b2max;
  double *q1st;
  double *q2st;
  double *om11st;
  double *om12st;
  double *om13st;
  double *om22st;
  double *cosx;

  double *mthree_pcharge;
  double *pc_scale_pcharge;
  double *dljpot_workspace;

  double *temp1;
  double *temp2;

  int64_t ask_for;
  int64_t one_l; 
  int64_t zero_l;
  int64_t l_64;
  int64_t ndoubles;
  
  int64_t align_size;
  int64_t align_mask;
  int64_t vec_len;
  int64_t vec_pad; 
  int64_t param_len;
  int64_t param_pad;
  int64_t num_data_points;
  int64_t num_parameters;
  int64_t num_inp_vecs;
  int64_t num_itn_vecs;
  int64_t num_ibstmax_vecs;
  int64_t dmax_offset;

  int success;
  int avec_len;

  int inp;
  int inp_align;

  int itn;
  int itn_align;

  int num_doubles;
  int ibstmax;

  int icoord;
  int imp;

  int ranlist_size;
  int itn_inp_imp_4;

  int inp_imp_4;
  int imp_4;

  int itn_inp;
  int itn_inp_imp;
  /*
    Setup to allocate vectors aligned on cache lines. 
  */
  one_l           = (int64_t)1;
  zero_l          = (int64_t)0;
  l_64            = (int64_t)64;
  align_size      = (int64_t)8;
  align_mask      = align_size - one_l;
  success         = 1;
  num_data_points = (int64_t)state->num_data_points;
  num_parameters  = (int64_t)state->num_parameters;
  icoord          = state->icoord;
  imp             = state->imp;
  itn             = state->itn;
  inp             = state->inp;
  /*
    NB we really only need (icoord * itn * inp * imp * 4) + 1, random numbers 
    but for good measure we generate an extra 7.
  */
  itn_inp         =  itn * inp;
  itn_inp_imp     =  itn_inp * imp;
  imp_4           =  imp << 2; 
  inp_imp_4       =  inp * imp_4;
  itn_inp_imp_4   =  itn * inp_imp_4;
  ranlist_size    =  (icoord * itn_inp_imp_4) + 8;
  state->ranlist_size = ranlist_size;
  state->imp_4    = imp_4;
  state->inp_imp_4 = inp_imp_4;
  state->itn_inp_imp_4 = itn_inp_imp_4;
  state->itn_inp   = itn_inp;
  state->itn_inp_imp = itn_inp_imp;
  if (num_parameters <= zero_l) {
    success = 0;
    fprintf(stderr,
	    "mobcall_alloc1: Error num_parameters must be > 0, value was %ld\n",
	    num_parameters);
    fflush(stderr);
  }
  if (success) {
    if (num_data_points <= zero_l) {
      success = 0;
      fprintf(stderr,
      "mobcall_alloc1: Error num_data_points must be > 0, value was %ld\n",
	      num_data_points);
      fflush(stderr);
    }
  }
  if (success) {
    /*
      ndoubles =  
      (number of fields with length num_data_points) * num_data_points
       	+ (number of fields with length num_parameters) * num_parameters
	+ num_data_point vectors: 
	   pcharge, fx,fy,fz,ox,oy,oz,eolj,rolj,eox4,ro6lj,
	   ro12lj,dro6,dro12,rhs,rhs2, xmass, imass = 18,
	   num_parameter vectors: tmc,tmm,ehsc,ehsm,pac,pam,asympp = 7. 
	   Even though imass is an integer we will store it as
	   a double.

	   So we will have the notion of two padding parmeters,
	   a num_data points pad vec_pad, and a 
	   parameter  padding param_pad, 

      We will require vec_len, and param_len to be multiples of 
 	align_size (which needs to be some power of 2 for the mask stuff 
	to work).

      Pad coordinate length out to an even multiple of 8.
      mobcal_mobil4 needs space for an incident ray, originally 
      stored at element inatom+1 of the fx,fy, and fz arrays.
    */
    vec_len = num_data_points + 1;
    vec_pad = (align_size - (vec_len & align_mask)) & align_mask;
    vec_len = vec_len + vec_pad;
    state->vec_len = vec_len;
    param_pad  = (align_size - (num_parameters & align_mask)) & align_mask;
    param_len = num_parameters + param_pad;
    state->param_len = param_len;
    /*   
      We neeed an original coordinate set, and
      a rotated coordinate set. corresponding to the ox, oy, oz, and
      fx, fy, fz, vectors in the fortran flavor.
    */
    /* 
       We need 6 vectors for ox,oy,oz, fx, fy,fz, and
       12 vectors for imass,xmass,pcharge,eolj,rolj,eox4, 
       ro6lj, ro12lj, dro6,dro12, rhs, and rhs2.

       But we will add another 11 vectors for the dljpot_input_stream,
       ordered as fx,fy,fz,pc_scale_pcharge,pcharge,mthree_pcharte,
       eox4,ro6lj,ro12lj,dro6,dro12
    ndoubles = ((18 * vec_len) + (8 * param_len) + ranlist_size + (2*itn*inp));
   */
    ndoubles = ((30 * vec_len) + (8 * param_len) + ranlist_size + (2*itn*inp));
    ask_for = ndoubles * (int64_t)sizeof(double);
    ox = (double *)aligned_calloc(one_l,ask_for,l_64);
    if (ox == NULL) {
      success = 0;
      fprintf(stderr,
	      "mobcal_alloc1: unable to allocate %ld bytes for data vectors\n",
	      ask_for);
      fflush(stderr);
    }
    if ((((int64_t)ox) & ((int64_t)63)) != 0) {
      success = 0;
      fprintf(stderr,"mobcal_alloc1: ox not 64 byte aligned ox = %px\n",ox);
      fflush(stderr);
    }
  }
  if (success) {
    state->usage += ask_for;
    itn_inp        = itn * inp;
    state->itn_inp = itn_inp;
    oy 	     	   = (double *)&ox[vec_len];
    oz 	     	   = (double *)&oy[vec_len];
    fx 	     	   = (double *)&oz[vec_len];
    fy 	     	   = (double *)&fx[vec_len];
    fz 	     	   = (double *)&fy[vec_len];
    imass          = (double *)&fz[vec_len];
    xmass          = (double *)&imass[vec_len];
    pcharge        = (double *)&xmass[vec_len];
    eolj           = (double *)&pcharge[vec_len];
    rolj     	   = (double *)&eolj[vec_len];
    eox4     	   = (double *)&rolj[vec_len];
    ro6lj    	   = (double *)&eox4[vec_len];
    ro12lj   	   = (double *)&ro6lj[vec_len];
    dro6     	   = (double *)&ro12lj[vec_len];
    dro12    	   = (double *)&dro6[vec_len];
    rhs      	   = (double *)&dro12[vec_len];
    rhs2     	   = (double *)&rhs[vec_len];

    tmc    	   = (double *)&rhs2[vec_len];
    tmm    	   = (double *)&tmc[param_len];
    asympp 	   = (double *)&tmm[param_len];
    ehsc   	   = (double *)&asympp[param_len];
    ehsm   	   = (double *)&ehsc[param_len];
    pac    	   = (double *)&ehsm[param_len];
    pam    	   = (double *)&pac[param_len];
    ranlist        = (double *)&pam[param_len+param_len];
    temp1  	   = (double *)&ranlist[ranlist_size];
    temp2  	   = (double *)&temp1[itn_inp];
    asym_work      = (double*)&temp2[itn_inp];
    dljpot_input_stream = (double*)&asym_work[vec_len];
    state->ox      = ox;
    state->oy      = oy;
    state->oz      = oz;
    state->fx      = fx;
    state->fy      = fy;
    state->fz      = fz;
    state->imass   = imass;
    state->xmass   = xmass;
    state->pcharge = pcharge;
    state->eolj    = eolj;
    state->rolj    = rolj;
    state->eox4    = eox4;
    state->ro6lj   = ro6lj;
    state->ro12lj  = ro12lj;
    state->dro6    = dro6;
    state->dro12   = dro12;
    state->rhs     = rhs;
    state->rhs2    = rhs2;
    state->tmc     = tmc;
    state->tmm     = tmm;
    state->asympp  = asympp;
    state->ehsc    = ehsc;
    state->ehsm    = ehsm;
    state->pac     = pac;
    state->pam     = pam;
    state->ranlist = ranlist;
    state->temp1   = temp1;
    state->temp2   = temp2;
    state->asym_work = asym_work;
    state->dljpot_input_stream = dljpot_input_stream;
  }
  /*
    But we also need some work vectors for mobil2:
    pgst, wgst, b2max, q1st, q2st, all of length inp,
    om11st, om12st, om13st, om22st of length itn,
    and cosx of length ibstmax
  */
  if (success) {
    inp = state->inp;
    inp_align = (align_size - (inp & align_mask)) & align_mask;
    inp += inp_align;
    itn = state->itn;
    itn_align = (align_size - (itn & align_mask)) & align_mask;
    itn += itn_align;
    ibstmax = state->ibstmax + 1; 
    /* we want to be able to refer to cos[ibstmax] */
    num_inp_vecs = (int64_t)5;

    num_itn_vecs = (int64_t)4;

    num_ibstmax_vecs = (int64_t)1;
    num_doubles = (int64_t)(num_inp_vecs * inp) + (num_itn_vecs * itn) + (int64_t)(num_ibstmax_vecs * ibstmax);
    ask_for = ndoubles * (int64_t)sizeof(double);
    pgst = (double *)aligned_calloc(one_l,ask_for,l_64);
    if (pgst == NULL) {
      success = 0;
      fprintf(stderr,
	      "mobcal_alloc1: unable to allocate %ld bytes for mobil2 scratch space\n",
	      ask_for);
      fflush(stderr);
    }
  }
  if (success) {
    state->usage += ask_for;
    wgst   = (double*)&pgst[inp];
    b2max  = (double*)&wgst[inp];
    q1st   = (double*)&b2max[inp];
    q2st   = (double*)&q1st[inp];
    om11st = (double*)&q2st[inp];
    om12st = (double*)&om11st[itn];
    om13st = (double*)&om12st[itn];
    om22st = (double*)&om13st[itn];
    cosx   = (double*)&om22st[itn];
    state->pgst   = pgst;
    state->wgst   = wgst;
    state->b2max  = b2max;
    state->q1st   = q1st;
    state->q2st   = q2st;
    state->q2st   = q2st;
    state->om11st = om11st;
    state->om12st = om12st;
    state->om13st = om13st;
    state->om22st = om22st;
    state->cosx   = cosx;
  }    
  if (success) {
    num_doubles = 2*vec_len;
    ask_for = num_doubles * (int64_t)sizeof(double);
    state->usage += ask_for;
    mthree_pcharge = (double *)aligned_calloc(one_l,ask_for,l_64);
    if (mthree_pcharge == NULL) {
      success = 0;
      fprintf(stderr,
	      "mobcal_alloc1: unable to allocate %ld bytes for mobil2 auxilliary vectors\n",
	      ask_for);
      fflush(stderr);
    } else {
      pc_scale_pcharge 	 = (double*)&mthree_pcharge[vec_len];
      state->mthree_pcharge = mthree_pcharge;
      state->pc_scale_pcharge = pc_scale_pcharge;
    }
  }
  if (success) {
    /*
      Here we shall allocate space for the dljpot_workspace, this would 
      replace all of the vectors in the above block except the mthree_pcharge,
      and pc_scale_pcharge vectors.
    num_doubles = 6*vec_len + 64 * align_size;
    for dljpot_workspace we need 41 * avec_len doubles.
    Here we round up to 64 and add some space to play with alignment
    issues as well.
    */
    avec_len = state->avec_len;
    num_doubles = 64*(avec_len + align_size);
    /*
    dmax_offset = 64*align_size;
    state->dmax_offset = dmax_offset;
    */
    ask_for = num_doubles * (int64_t)sizeof(double);
    state->usage += ask_for;
    dljpot_workspace = (double *)aligned_calloc(one_l,ask_for,l_64);
    if (dljpot_workspace == NULL) {
      success = 0;
      fprintf(stderr,
	      "mobcal_alloc1: unable to allocate %ld bytes for "
	      "dljpot workspace",
	      ask_for);
      fflush(stderr);
    } else {
      state->dljpot_workspace = dljpot_workspace;
    }
  }
  return(success);
}

