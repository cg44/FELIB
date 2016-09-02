#include "f2c.h"


int asful_(doublereal *sysk, integer *isysk, integer *jsysk, doublereal *elk, 
		   integer *ielk, integer *jelk, integer *steer, integer *isteer, 
		   integer *dofel, integer *itest)
{

	integer errmes_(integer *itest, integer *ierror, char *srname, ftnlen srname_len);

    /* Initialized data */

    static char srname[6+1] = "ASFUL ";

    /* System generated locals */
    integer elk_dim1, elk_offset, sysk_dim1, sysk_offset, i__1, i__2;

    /* Local variables */
    static integer i, j, jtest, steeri, steerj;
    extern integer errmes_();
    static integer ierror;

/* 
-----------------------------------------------------------------------
 PURPOSE

      ASFUL assembles full real system matrix 

 HISTORY 

      Copyright (C) 2000 : CCLRC, Rutherford Appleton Laboratory
                           Chilton, Didcot, Oxfordshire OX11 0QX 

      Release 1.1    29 Oct 1978 (CG)
      Commented       6 Feb 1980 (KR)
      Release 4.0     2 Oct 1999 (CG)
      Release 1.0 (C) 1 May 2000 (CG)

 ARGUMENTS in 
      SYSK    contains system matrix prior to addition of 
              current elemnt matrix contribution 
      ISYSK   first dimension of SYSK (.GE. total number of 
              unconstrained degrees of freedom) 
      JSYSK   second dimension of SYSK (.GE. total number of 
              unconstrained degrees of freedom) 
      ELK     element matrix 
      IELK    first dimension of ELK (.GE. DOFEL) 
      JELK    second dimension of ELK (.GE. DOFEL) 
      STEER   contains freedom numbers associated with element 
              matrix contributions to system matrix 
      ISTEER  dimension of STEER (.GE. DOFEL) 
      DOFEL   maximum number of degrees of freedom associated 
              with element type 
      ITEST   error checking option 

 ARGUMENTS out 
      SYSK    system matrix 

 ROUTINES called 
      ERRMES 

***********************************************************************
*/
/* Parameter adjustments */
    --steer;
    elk_dim1 = *ielk;
    elk_offset = elk_dim1 + 1;
    elk -= elk_offset;
    sysk_dim1 = *isysk;
    sysk_offset = sysk_dim1 + 1;
    sysk -= sysk_offset;

    /* Function Body */

/*     Parameter checking */

    jtest = *itest;
    if (jtest != -1) {
	ierror = 0;
	if (*isteer < *dofel) {
	    ierror = 3;
	}
	if (*ielk < *dofel || *jelk < *dofel) {
	    ierror = 2;
	}
	if (*dofel <= 0) {
	    ierror = 1;
	}
	*itest = errmes_(&jtest, &ierror, srname, 6L);
	if (*itest != 0) {
	    return 0;
	}
    }

/*     Main loops */

    i__1 = *dofel;
    for (i = 1; i <= i__1; ++i) {
	steeri = steer[i];
	if (steeri > 0) {
	    i__2 = *dofel;
	    for (j = 1; j <= i__2; ++j) {
		steerj = steer[j];
		if (steerj > 0) {

/*     Range checking on STEERI and STEERJ */

		    if (jtest != -1) {
			ierror = 0;
			if (*isysk < steeri || *jsysk < steerj) {
			    ierror = 4;
			}
			*itest = errmes_(&jtest, &ierror, srname, 6L);
			if (*itest != 0) {
			    return 0;
			}
		    }
		    sysk[steeri + steerj * sysk_dim1] += elk[i + j * elk_dim1];
		}
	    }
	}
    }
return 0;
} 
/* asful_ */

