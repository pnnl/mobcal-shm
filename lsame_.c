#include "system_includes.h"
#include "lsame.h"
int lsame_(char *ca_p, char *cb_p) {
/*
  lsame.c adapted from lapack lsame.f

   Arguments: ca_p pointer to character ca
              cb_p pointer to character cb
   Returns:
   LSAME returns 1 if ca is the same letter as cb regardless of
   case.
*/
  int  inta;
  int  intb;
  int  result;
  int  izcode;
  char ca;
  char cb;
  char zcode;
  char cpad;
  ca = *ca_p;
  cb = *cb_p;
  zcode = 'Z';
  result = (ca == cb);
  if (result == 0) {
    inta = (int)ca;
    intb = (int)cb;
    izcode = (int)zcode;
    if ((izcode == 90) || (izcode == 122)) {
      /*
        ASCII is assumed - ZCODE is the ASCII code of either lower or
        upper case 'Z'.
      */
      if ((inta >= 97) && (inta <=122)) inta = inta - 32;
      if ((intb >= 97) && (intb <=122)) intb = intb - 32;
    } else {
      if ((izcode == 233) || (izcode == 169)) {
/*
        EBCDIC is assumed - ZCODE is the EBCDIC code of either lower or
        upper case 'Z'.
*/

	  if ( ((inta >= 129) && (inta<=137)) ||
	       ((inta >= 145) && (inta<=153)) ||
	       ((inta >= 162) && (inta<=169))) inta = inta + 64;
	  if ( ((intb >= 129) && (intb<=137)) ||
	       ((intb >= 145) && (intb<=153)) ||
	       ((intb >= 162) && (intb<=169))) intb = intb + 64;
      } else {
        if ((izcode == 218) || (izcode == 250 )) {
/*
           ASCII is assumed, on Prime machines - ZCODE is the ASCII code
           plus 128 of either lower or upper case 'Z'.
*/
	  if ( (inta >= 225) && (inta <= 250) ) inta = inta - 32;
	  if ( (intb >= 225) && (intb <= 250) ) intb = intb - 32;
	} /* end else Prime machine */
      } /* end else not EBCDIC */
    } /* end else not asci */
    result = (inta == intb);
  } /* end if (result == 0) */
  return(result);
}
