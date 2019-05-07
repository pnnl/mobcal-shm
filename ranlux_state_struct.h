#ifndef _RANLUX_STATE_STRUCT_H_
#define _RANLUX_STATE_STRUCT_H_ 1
struct ranlux_state_struct {
  double *seeds; /* length = 24 */

  double carry;
  double twom48;

  double twom24;
  double twom12;

  int *inext;   /* length = 24 */
  int *ndskip;  /* length = 8 although we only use ther first 5 */
  int *isdext;  /* length = 32 although only needs to be 25 */

  int uninitialized;
  int inseed;

  int maxlev;
  int lxdflt;

  int igiga;
  int jsdflt;

  int icons;

  int in24;
  int nskip;

  int i24;
  int j24;

  int kount;
  int mkount;

  int luxlev;
  int padi;
  
  int padj;
  int padk;
}
;
#endif
