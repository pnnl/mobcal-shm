#include "system_includes.h"
#include "mobcal_state_struct.h"
/*
#include "mobcal_xrand.h"
#include "mobcal_rantate.h"
#include "mobcal_rantate2.h"
#include "mobcal_gsang.h"
*/
#include "mobcal_vec_set.h"
#include "mobcal_acc_temps.h"
#include "mobcal_mobil2_server.h"
void mobcal_mobil2_server(struct mobcal_state_struct *state) {
  /*
    Schedule the main mobcal loops, over ic[1:itn]  ig[1:inp], and im[1:imp]
    sets the om11st, om12st, om13st and om22st, q1st, w2st, fields of state.
    Assumes ip != 1 and igs != 1
    Called by: mobcal_mobil2
    Calls:     fprintf, fflsuh, mobcal_rantate, mobcal_gsang, sin
  */
  void *server_area_base;
  void *client_area_base;
  volatile void *server_area;
  volatile void *client_area;
  volatile int64_t *client_flag1;
  volatile int64_t *jclient_area;
  volatile int64_t *jserver_area;
  volatile double  *dserver_area;
  volatile double  *dmaster_area;
  double *om11st;
  double *om12st;
  double *om13st;
  double *om22st;
  double *q1st;
  double *q2st;
  double *b2max;
  double *pgst;
  double *wgst;
  double *ranlist;
  double *temp1;
  double *temp2;
  double pgst_ig;
  double wgst_ig;
  double pgst_ig2;
  double pgst_ig4;
  double v;
  double eo;
  double ro;
  double recip_imp;
  double recip_3_tst;
  double recip_12_tst_2;
  double recip_half_mu;
  double temp1v;
  double temp2v;
  double igval1;
  double igval2;
  double dzero;

  int64_t zero_l;
  int64_t one_l;
  int64_t server_flag1_p1;
  int64_t dist_work_poll_miss;
  int64_t cleanup_poll_miss;

  int     next_client;
  int     client_ready;
  int     client_finished;
  int     completed_tasks;     


  int im2;
  int ip;

  int inp;
  int itn;

  int imp;
  int igs;

  int ig;
  int ic;

  int im;
  int tot_num_points;

  int iic;
  int padi;

  int inp_imp_4;
  int imp_4;

  int itn_inp_imp_4;
  int itn_inp;

  int shm_block_size;
  int my_thread_id;

  int num_threads;
  int itn_inp_imp;

  int i;
  int ipos;

  FILE *ofp;
  FILE *lfp;
  /*
#define DBG 1
  */
  zero_l         = (int64_t)0;
  one_l          = (int64_t)1;
  dzero          = 0.0;
  im2 		 = state->im2;
  ip  		 = state->ip;
  inp 		 = state->inp;
  itn 		 = state->itn;
  imp 		 = state->imp;
  igs 		 = state->igs;
  eo             = state->eo;
  ro             = state->ro;
  q1st 		 = state->q1st;
  q2st 		 = state->q2st;
  om11st 	 = state->om11st;
  om12st 	 = state->om12st;
  om13st 	 = state->om13st;
  om22st 	 = state->om22st;
  b2max  	 = state->b2max;
  pgst   	 = state->pgst;
  wgst   	 = state->wgst;
  recip_3_tst    = state->recip_3_tst;
  recip_12_tst_2 = state->recip_12_tst_2;
  recip_half_mu  = state->recip_half_mu;
  itn_inp        = state->itn_inp;
  itn_inp_imp    = state->itn_inp_imp;
  itn_inp_imp_4  = state->itn_inp_imp_4;
  inp_imp_4      = state->inp_imp_4;
  imp_4          = state->imp_4;
  ranlist        = state->ranlist;
  ofp 		 = state->ofp;
  lfp            = state->lfp;
  recip_imp = 1.0/imp;
  iic = state->iic;
#ifdef DBG
  if (lfp) {
    fprintf(lfp,"mobcal_mobil2_server top\n");
    fflush(lfp);
  }
#endif  
  if (im2 == 0) {
    tot_num_points = itn * inp * imp;
    if (ofp) {
      fprintf(ofp,"number of complete cycles (itn) = %d\n",itn);
      fprintf(ofp,"number of velocity points (inp) = %d\n",inp);
      fprintf(ofp,"number of random points (imp) = %d\n",imp);
      fprintf(ofp,"total number of random points = %d\n",tot_num_points);
      fflush(ofp);
    }
  }
  for (ig=0;ig<inp;ig++) {
    q1st[ig] = 0;
    q2st[ig] = 0;
  }
  for (ic =0;ic <itn;ic++) {
    om11st[ic] = 0.0;
    om12st[ic] = 0.0;
    om13st[ic] = 0.0;
    om22st[ic] = 0.0;
  }    
  server_area_base = state->server_area_base;
  client_area_base = state->client_area_base;
  shm_block_size   = state->shm_block_size;
  my_thread_id     = state->thread_id;
  num_threads      = state->num_threads;
  dmaster_area = (double*) client_area_base;
  temp1        = (double *)&dmaster_area[8];
  temp2        = (double *)&temp1[itn_inp];
  mobcal_vec_set(itn_inp,temp1,dzero);
  mobcal_vec_set(itn_inp,temp2,dzero);
  /*
    Initialize communication with clients.
  */
  server_area = server_area_base + shm_block_size;
  for (i=1;i<num_threads;i++) {
    jserver_area = (int64_t*)server_area;
    jserver_area[_FLAG1_] = 0;
    jserver_area[_FLAG2_] = 0;
    server_area += shm_block_size; /* Caution address arithmetic here */
  }

  /*
    We will switch the order of the loops.
  */
  server_area   = server_area_base + shm_block_size;
  client_area   = client_area_base + shm_block_size;
  jserver_area  = (int64_t *) server_area;
  dserver_area  = (double  *) server_area;
  jclient_area  = (int64_t *) client_area;
  next_client   = 1;
  completed_tasks  = 0;
  dist_work_poll_miss = zero_l;
#ifdef DBG
  if (lfp) {
    fprintf(lfp,"mobcal_mobil2_server before main_loop, server_area = %p, client_area = %p, next_client = %d, completed_tasks = %d\n",server_area,client_area,next_client,completed_tasks);
    fflush(lfp);
  }
#endif
  /*
    main loop.
  */
  for (ig = 0; ig < inp; ig++) {
    pgst_ig = pgst[ig];
    v    = pgst_ig * sqrt(eo*recip_half_mu);
    for (ic = 0;ic < itn;ic++) {
      for (im = 0; im < imp; im++) {
	/*
	  Polling loop. We might want to count the number of poll 
	  misses.
	*/
	client_ready = 0;
	while (client_ready == 0) {
	  /*
	    Look for a ready client
	  */
	  server_flag1_p1 = jserver_area[_FLAG1_] + one_l;
	  client_flag1 = &jclient_area[_FLAG1_];
	  if (server_flag1_p1 != *client_flag1) {
	    /*
	      Client is not ready.
	    */
	    dist_work_poll_miss += 1;

	  } else {
	    if (jclient_area[_IIC_] >= 0) {
	      completed_tasks += 1;
	    }
	    jserver_area[_IIC_] = iic;
	    jserver_area[_IC_]  = ic;
	    jserver_area[_IG_]  = ig;
	    jserver_area[_IM_]  = im;
	    dserver_area[_V_]   = v;
#ifdef DBG
  if (lfp) {
    fprintf(lfp,"mobcal_mobil2_server sending iic = %d, ic = %d, ig = %d, im = %d, v = %le to %d\n",iic,ic,ig,im,v,next_client);
    fflush(lfp);
  }
#endif
	    /*
	      Might need some kine of msync in heare to force these fields to
	      be flushed out to memory where client can see.
	    jserver_area[_FLAG1_] = server_flag1_p1;
	    */
	    jserver_area[_FLAG1_] += one_l;
	    client_ready = 1;
	  }
	  next_client += 1;
	  if (next_client == num_threads) {
	    server_area = server_area_base;
	    client_area = client_area_base;
	    next_client = 1;
	  }
	  server_area += shm_block_size; /* Caution address arithmetic */
	  client_area += shm_block_size;
	  jserver_area = (int64_t *)server_area;
	  dserver_area = (double  *)server_area;
	  jclient_area = (int64_t *)client_area;
	} /* end while (client_ready == 0) */
      } /* end for im */
    } /* end for ic */
  } /* end for ig */
#ifdef DBG
  if (lfp) {
    fprintf(lfp,"mobcal_mobil2_server After distribute tasks loop\n");
    fflush(lfp);
  }
#endif
  cleanup_poll_miss = (int64_t)0;
  for (i = completed_tasks; i < itn_inp_imp; i ++) {
    client_finished = 0;
    while (client_finished == 0) {
      /*
	Look for a ready client
      */
      server_flag1_p1 = jserver_area[_FLAG1_] + one_l;
      client_flag1 = &jclient_area[_FLAG1_];
      if (server_flag1_p1 != *client_flag1) {
	/*
	  Client is not finished.
	*/
	cleanup_poll_miss += 1;
      } else {
	if (jclient_area[_IIC_] >= 0) {
	  client_finished = 1;
	}
	jserver_area[_IIC_] = -1;
	/*
	jserver_area[_FLAG1_] = server_flag1_p1;
	*/
	jserver_area[_FLAG1_] += one_l;
      }
      next_client += 1;
      if (next_client == num_threads) {
	server_area = server_area_base;
	client_area = client_area_base;
	next_client = 1;
      }
      server_area += shm_block_size; /* Caution address arithmetic here. */
      client_area += shm_block_size; /* Caution address arithmetic here. */
      jserver_area = (int64_t *)server_area;
      jclient_area = (int64_t *)client_area;
    } /* end while (client_finished == 0) */
  } /* end for (i=completed_tasks ... ) */
#ifdef DBG
  if (lfp) {
    fprintf(lfp,"mobcal_mobil2_server after cleanup loop.\n");
    fflush(lfp);
  }
#endif
  /*
    Now we will need to accumulate results for temp1 and temp2 
    arrays accross processors.
    Here 
  */
  mobcal_acc_temps(state);
  ipos = 0;
  for (ig = 0; ig < inp; ig++) {
    pgst_ig = pgst[ig];
    wgst_ig = wgst[ig];
    pgst_ig2 = pgst_ig * pgst_ig;
    pgst_ig4 = pgst_ig2 * pgst_ig2;
    igval1 = pgst_ig2 * wgst_ig * recip_3_tst;
    igval2 = pgst_ig4 * wgst_ig * recip_12_tst_2;
    for (ic = 0;ic < itn;ic++) {
#ifdef DBG
      if (lfp) { 
	fprintf(lfp,"mobcal_mobil2_server: ig = %d, ic = %d, ipos = %d,temp1[ipos] = %le, temp2[ipos] = %le\n",ig,ic,ipos,temp1[ipos],temp2[ipos]);
	fflush(lfp);
      }
#endif      
      temp1v = temp1[ipos];
      temp2v = temp2[ipos];
      om11st[ic] += (temp1v * wgst_ig);
      om12st[ic] += (temp1v * igval1);
      om13st[ic] += (temp1v * igval2);
      om22st[ic] += (temp2v * igval1);
      q1st[ig]   += temp1v;
      q2st[ig]   += temp2v;
      ipos += 1;
    } /* end for (ic...) */
  } /* end for (ig...) */
  state->dist_work_poll_miss = dist_work_poll_miss;
  state->cleanup_poll_miss   = cleanup_poll_miss;
}
