#ifndef _MOBCAL_STATE_STRUCT_H_
#define _MOBCAL_STATE_STRUCT_H_ 1
#include "ranlux_state_struct.h"
#include "mt_state_struct.h"
#include "mobcal_diffeq_state_struct.h"
/*
  Unrolling Parameter, for dljpot _IVEC_LEN_
*/
#define _IVEC_LEN_ 4
#define _SIX_IVEC_LEN_ 24
#define _AVEC_LEN_ 8
/*
  Coordinate offsets.
*/
#define _X_ 0
#define _Y_ 1
#define _Z_ 2
/*
  Properties fields offsets.
*/
#define _IMASS_   0
#define _XMASS_   1
#define _PCHARGE_ 2
#define _EOLJ_    3
#define _EOX4_    4
#define _ROLJ_    5
#define _RO6LJ_   6
#define _RO12LJ_  7
#define _DRO6_    8
#define _DRO12_   9
#define _RHS_     10
#define _RHS2_    11 
/*
  gprops field offsets.
*/
#define _GMASS_ 0
#define _GEOLJ_ 1
#define _GROLJ_ 2
#define _GRHS_  3

/*
  Thread communication fields 
*/
#define _IIC_ 0
#define _IC_  1
#define _IG_  2
#define _IM_  3
#define _DATA_OFFSET_ 4
#define _V_   5
#define _FLAG1_ 6
#define _FLAG2_ 7


struct mobcal_state_struct {
  void *server_area_base;
  void *client_area_base;
  void *my_server_area;
  void *my_client_area;
  struct ranlux_state_struct *ranlux_state;
  struct mt_state_struct *mt_state;
  struct mobcal_diffeq_state_struct *diffeq_state;
  /*
    Vector to store all random numbers generated.
    allocated in alloc2.
  */
  double *ranlist;
  /* 
    coordinates, lengths = number_data_points * coord_lda
  */
  double *ox;
  double *oy;
  double *oz;
  /*
  double *orig_coords;
  */
  double *fx;
  double *fy;
  double *fz;
  /*
  double *rota_coords;
  */
  /*
    Scratch space needed by mobcal_mobil4,
    it rotates the coordinates and 
  double *rota2_coords;
  */
  /*
    ljparameters, lengths = number_data_points
  */
  /*
    Properties length = number_data_points * prop_lda
  */
  double *imass; 
  double *xmass;
  double *pcharge;
  double *eolj;
  double *eox4;
  double *rolj;
  double *ro6lj;
  double *ro12lj;
  double *dro6;
  double *dro12;
  double *rhs;
  double *rhs2;
  /*
    hsparameters, lengths = number_data_points.
  double *properties;
  */
  /*
    atom type Parameters for xmass, eolj, rolj, and rhs 
  */
  double *gprops;
  /*
    output doubles, lengths = number_parameters
  */
  double *tmc;
  double *tmm;
  double *asympp;
  double *ehsc;
  double *ehsm;
  double *pac;
  double *pam;
  /*
    vectors needed by mobil2
  */
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
  /*
    pcharge scalings of length inatoms.
    used in dljpot
  */
  double *mthree_pcharge;
  double *pc_scale_pcharge;

  double *asym_work;
  double *dljpot_input_stream;
  double *dljpot_workspace;
  /*
    Parallel accumulation vectors for mobil2.
  */
  double *temp1;
  double *temp2;
  /*
    position doubles (center of mass coordinates).
  */
  double fxo;
  double fyo;
  double fzo;
  /*
    tracjectory doubles
  */
  double sw1;
  double sw2;
  double dtsf1;
  double dtsf2;
  double cmin;
  /*
    angles 
  */
  double theta;
  double phi;
  double gamma;
  /*
    constant doubles.
  */
  double mu;
  double recip_half_mu;
  double ro;
  double eo;
  double pi;
  double two_pi;
  double four_pi;
  double eight_pi;
  double half_pi;
  double quarter_pi;
  double cang;
  double ro2;
  double pi_ro2;
  double dipol;
  double emax;
  double xe;
  double xeo;
  double xe2_o_4pi_xeo;
  double pc;
  double pc_scale;
  double xmv;
  double xk;
  double xn;
  double correct;
  double romax;
  double m1;
  double m2;
  double mconst;
  double recip_mconst;
  double temp;
  double tst;
  double tst3;
  double recip_3_tst;
  double recip_12_tst_2;
  double romax_adj;
  double bg_dipole_mult;
  /*
    TRAJONE doubles;
  */
  double v;
  double b;
  double ntheta;
  double nphi;
  double ngamma;
  /*
    space usage data.
  */
  int64_t vec_len;
  int64_t param_len;
  int64_t usage;
  int64_t dmax_offset;
  /*
    Spin wait and polling stats.
  */
  /*
    These 2 set by mobcal_mobil2_server
  */
  int64_t dist_work_poll_miss;
  int64_t cleanup_poll_miss;
  /*
    These 3 set mob mobcal_acc_temp2
  */
  int64_t client_spin_count;
  int64_t left_spin_count;
  int64_t right_spin_count;
  /*
    This one set by mobcal_mobil2_client or mobcal_mobil2_client_dbg.
  */
  int64_t finish_spin_count;

  int num_data_points;
  int num_parameters;

  int thread_id;
  int num_threads;

  int shmid;
  int shm_block_size;

  int shm_size;
  int shm_prjid;

  int mpi_rank;
  int num_mpi_ranks;
  
  int ipr; /* number of random rotations per potential calculation.*/
  int itn; /* Nomber of complete cycles for average mobility calc in mobil2 */

  int inp; /* Number of points in velocity integration in mobil2.*/
  int imp; /* Number of points in Monte Carlo integrations of impact
	      parameter and orienation in mobil2. */

  // random number generator parameters.
  int i1;
  int i2; /* should change to random_seed. */

  int i3; /* kount */
  int i4; /* mkount */

  int i5; /* rn_choice */
  int i6; /* */
  /*
    trajectory ints;
  */
  int ifail;
  int ifailc;

  int inwr;
  /*
    constant integers.
  */
  int inatom;

  int icoord;
  int iic;

  /* 
    Monte Carlo integers
  */
  int inum;
  int inor;

  // print switches
  int ip;
  int it;

  int iu1;
  int iu2;

  int iu3;
  int iv;

  int im2;
  int im4;

  int igs;
  int iunit;
  
  int charge_dist;
  int imm;

  int immmin;
  int immmax;

  int max_atom_imass;
  int input_buffer_len;

  int ibstmax;
  int irn; /* see mobcal_rmax_emax_r00 */

  int ranlist_size;
  int itn_inp;

  int itn_inp_imp;
  int itn_inp_imp_4;

  int inp_imp_4;  
  int imp_4;

  int mobility_type;
  /* Use Dennis Thomas' eolj,rolj scheme, if 1, else use iu scheme */
  int use_dgt; 

  /* Atomic number of the buffer gas, used to set m1 and dipol in
     mobcal_init_constants 1. 14 for nitrogen, 4 for  helium 
     else nitrogen is assumed. */
  int buffer_gas;
  /*
    Use the mersenne-twister random number generator if 1, ranlux otherwise.
  */
  int use_mt;

  int filename_limit;
  /*
    Use the Indianna University LJ schema if set to 1, otherwise use
    Dennis Thomas' view.
  */
  int use_iu_dljpot;

  int avec_len;
  int padi;

// input characters.
  char *param_file;
  char *at_param_file;
  char *mfj_file; // formerly filen1.
  char *output_file; // formerly filen2.
  /*
  char *filen1;
  char *filen2;
  */
  char *unit;
  char *dchar;
  char *xlabel;
  char *input_buffer; /* buffer for reading atomtype_parameters and 
			 coordinates files. */

  FILE *ifp;
  FILE *ofp;
  FILE *lfp;
  FILE *efp;
}
;
#endif
