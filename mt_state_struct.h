#ifndef _MT_STATE_STRUCT_H_
#define _MT_STATE_STRUCT_H_ 1
#define NN 312
#define MM 156
#define MATRIX_A 0xB5026F5AA96619E9ULL
#define UM 0xFFFFFFFF80000000ULL /* Most significant 33 bits */
#define LM 0x7FFFFFFFULL /* Least significant 31 bits */

struct mt_state_struct {
  unsigned long long mt[NN];
  unsigned long long mag01[2];
  int mti;
  int padi;
}
;
#endif
