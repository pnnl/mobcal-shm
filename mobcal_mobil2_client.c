#include "system_includes.h"
#include "mobcal_state_struct.h"
/*
#include "mobcal_xrand.h"
#include "mobcal_rantate.h"
*/
#include "mobcal_rantate2.h"
#include "mobcal_gsang.h"
#include "mobcal_vec_set.h"
#include "mobcal_acc_temps.h"
#include "mobcal_mobil2_client.h"
void mobcal_mobil2_client(struct mobcal_state_struct *state) {
  /*
    Schedule the main mobcal loops, over ic[1:itn]  ig[1:inp], and im[1:imp]
    sets the om11st, om12st, om13st and om22st, q1st, w2st, fields of state.
    Assumes ip != 1 and igs != 1
    Called by: mobcal_shm_mobil2
    Calls:     fprintf, fflsuh, mobcal_rantate, mobcal_gsang, sin
  */
  void *server_area_base;
  void *client_area_base;
  volatile void *my_server_area;
  volatile void *my_client_area;
  volatile int64_t *jclient_area;
  volatile int64_t *jserver_area;
  volatile int64_t *server_flag1;
  volatile double *dclient_area;
  volatile double *dserver_area;
  volatile double *temp1;
  volatile double *temp2;
  
  double *b2max;
  double *ranlist;
  double v;
  double eo;
  double b;
  double ro;
  double erat;
  double ang;
  double d1;
  double rnb;
  double bst2;
  double hold1;
  double sin_ang;
  double hold2;
  double recip_imp;
  double recip_3_tst;
  double recip_12_tst_2;
  double recip_half_mu;
  double zero;
  double temp1d;
  double temp2d;

  int64_t zero_l;
  int64_t one_l;
  int64_t client_spin_count;
  int im2;
  int ip;

  int inp;
  int itn;

  int imp;
  int igs;

  int ig;
  int ic;


  int istep;
  int iic;

  int iran;
  int itn_inp_imp_4;
  int inp_imp_4;
  int imp_4;

  int iranp1;
  int itn_inp;

  int shm_block_size;
  int my_thread_id;

  int num_threads;
  int temp_pos;

  int my_flag1;
  int not_done;

  int im;
  int iran_off;

  int use_iu_dljpot;
  int padi;

  FILE *ofp;
  FILE *lfp;
  /*
#define DBG 1
  */
  zero_l         = (int64_t)0;
  one_l          = (int64_t)1;
  zero           = 0.0;
  im2 		 = state->im2;
  ip  		 = state->ip;
  inp 		 = state->inp;
  itn 		 = state->itn;
  imp 		 = state->imp;
  igs 		 = state->igs;
  eo             = state->eo;
  ro             = state->ro;
  ofp 		 = state->ofp;
  b2max  	 = state->b2max;
  recip_3_tst    = state->recip_3_tst;
  recip_12_tst_2 = state->recip_12_tst_2;
  recip_half_mu  = state->recip_half_mu;
  itn_inp        = state->itn_inp;
  itn_inp_imp_4  = state->itn_inp_imp_4;
  inp_imp_4      = state->inp_imp_4;
  imp_4          = state->imp_4;
  ranlist        = state->ranlist;
  lfp            = state->lfp;
  use_iu_dljpot  = state->use_iu_dljpot;
  recip_imp = 1.0/imp;
  iic = state->iic;
  server_area_base = state->server_area_base;
  client_area_base = state->client_area_base;
  shm_block_size   = state->shm_block_size;
  my_thread_id     = state->thread_id;
  num_threads      = state->num_threads;
  my_server_area   = server_area_base + my_thread_id * shm_block_size;
  my_client_area   = client_area_base + my_thread_id * shm_block_size;
  jclient_area     = (int64_t *)my_client_area;
  dclient_area     = (double *)my_client_area;
  jserver_area     = (int64_t *)my_server_area;
  dserver_area     = (double *)my_server_area;
  server_flag1     = &jserver_area[_FLAG1_];
  temp1            = &dclient_area[8];
  temp2            = &temp1[itn_inp];
  mobcal_vec_set(itn_inp,temp1,zero);
  mobcal_vec_set(itn_inp,temp2,zero);
  /*
    Start by setting ready to communicate message.
  */
  jclient_area[_IIC_] = -1;
  jclient_area[_FLAG1_] = 1;
  my_flag1 = 1;
  not_done = 1;
  iran_off  = 4;
  if (use_iu_dljpot == 1) {
    /*
      If using the iu dljpot function, match their use
      of random numbers.
    */
    iran_off = 1;
    /*
      When we add mobil4 this will becom 
      iran_off = 1+5*inum to match iu code.
    */
  }
  client_spin_count = zero_l;
  while (not_done) {
    /*
      Wait for server message.
    */
    while (*server_flag1 != my_flag1) {
      client_spin_count += one_l;
    }
    iic = jserver_area[_IIC_];
    if (iic < zero_l) {
      /*
	All tasks have been handed out exit loop.
      */
      not_done = 0;
    } else {
      ic = jserver_area[_IC_];
      ig = jserver_area[_IG_];
      im = jserver_area[_IM_];
      v  = dserver_area[_V_];
      jclient_area[_IIC_] = iic;
      jclient_area[_IC_]  = ic;
      jclient_area[_IG_]  = ig;
      jclient_area[_IM_]  = im;
      dclient_area[_V_]   = v;
      state->v = v;
      temp_pos = (ig * itn) + ic;
#ifdef DBG
      if (lfp) {
	fprintf(lfp,"mobcal_mobil2_client, iic = %d, ic = %d, ig = %d, im = %d, v = %le, temp_pos = %d\n",iic,ic,ig,im,v,temp_pos);
	fflush(lfp);
      }
#endif

      /*
	rnb = mobcal_xrand(state);
	mobcal_rantate(state);
      */
      /*
      iran = (iic * itn_inp_imp_4) + (ic * inp_imp_4) + (ig *  imp_4) +
	(im * 4) + 1;
      */
      iran = (iic * itn_inp_imp_4) + (ic * inp_imp_4) + (ig *  imp_4) +
	(im * 4) + iran_off;
      rnb = ranlist[iran];
#ifdef DBG
      if (lfp) {
	fprintf(lfp,"iran = %d, rnb = %le, ",iran,rnb);
	fflush(lfp);
      }
#endif
      iranp1 = iran + 1;
      mobcal_rantate2(state,iranp1);
      bst2 = rnb * b2max[ig];
      b    = ro *sqrt(bst2);
      state->b = b;
      mobcal_gsang(state,v,b,&erat,&ang,&d1,&istep);
	/*
	  Right here I want to look at istep, ic, ig, and im, v and b
	  fprintf(stdout,"istep = %d, ic = %d, ig = %d, im = %d, "
	          "v = %le, b = %le\n",istep,ic,ig,im,v,b);
          fflush(stdout);		  
	*/
      hold1 = 1.0 - cos(ang);
      sin_ang = sin(ang);
      hold2 = sin_ang * sin_ang;
      /*
	temp1 += ((hold1 * b2max[ig])*recip_imp);
	temp2 += (((1.5 * hold2) * b2max[ig])*recip_imp);
      */
      /*
	Linear Position in the temp* array for element [ic,ig] will
	be ig * itn + ic (fortran column oriented order).
      */
      temp1d = ((hold1 * b2max[ig])*recip_imp);
      temp1[temp_pos] += temp1d;
      temp2d = (((1.5 * hold2) * b2max[ig])*recip_imp);
      temp2[temp_pos] += temp2d;
#ifdef DBG
      if (lfp) {
	fprintf(lfp,"mobcal_mobil2_client: ic = %d, ig = %d, im = %d, ang = %le,temp1d = %le, temp2d = %le\n",ic,ig,im,ang,temp1d,temp2d);
	fflush(ofp);
      }
#endif
      /*
	Signal task finish by incrementing client flag1 location.
      */
      /*
      jclient_area[_FLAG1_] = my_flag1;
      */
      my_flag1 += 1;
      jclient_area[_FLAG1_] += one_l;
    } /* end else we had a task */
  } /* end while not_done */
  /*
    Now we will need to accumulate results for temp1 and temp2 
    arrays accross processors.
    Here 
  */
#ifdef DBG
  if (lfp) {
    fprintf(lfp,"mobcal_mobil2_client before mobcal_acc_temps\n");
    fflush(lfp);
  }
#endif

  mobcal_acc_temps(state);
  state->client_spin_count = client_spin_count;
}
