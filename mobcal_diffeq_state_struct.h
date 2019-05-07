#ifndef _MOBCAL_DIFFEQ_STATE_STRUCT_H_
#define _MOBCAL_DIFFEQ_STATE_STRUCT_H_ 1
struct mobcal_diffeq_state_struct {
  double *a;        /* length 4 */
  double *b;        /* length 4 */
  double *c;        /* length 4 */
  double *three_a;  /* length 4 */
  double *ampc;     /* length 5 */
  double *amcc;     /* length 4 */
  double *savw;     /* length 6 */
  double *savdw;    /* length 6 */
  double *q;        /* length 6 */
  double *r;        /* length 6 */
  double *array;    /* length 36 */
  double *predictor; /* points to column 5 of array (starting from 0). */
  double *corrector; /* points to column 5 of array (starting from 0). */
  double *col_4;     /* points to column 4 of array (starting from 0). */
  double *col_1;     /* points to column 1 of array (starting from 0). */
  double var;
  double cvar;
  double acst;
  double hvar;
  double hcvar;
  int    ld_array;
  int    four_ld_array;
}
;
#endif
