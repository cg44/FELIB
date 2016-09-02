/* asusmg.f -- translated by f2c (version of 28 August 1991  0:07:02).
   You must link the resulting object file with the libraries:
	-lF77 -lI77 -lm -lc   (in that order)
*/

#include "f2c.h"


/* Subroutine */ int asusmg_(sysk, isysk, jsysk, elk, ielk, jelk, str1, istr1,
	 str2, istr2, hband, dofel1, dofel2, itest)
doublereal *sysk;
integer *isysk, *jsysk;
doublereal *elk;
integer *ielk, *jelk, *str1, *istr1, *str2, *istr2, *hband, *dofel1, *dofel2, 
	*itest;
{
    /* Initialized data */

    static char srname[6+1] = "ASUSMG";

    /* System generated locals */
    integer elk_dim1, elk_offset, sysk_dim1, sysk_offset, i__1, i__2;

    /* Local variables */
    static integer i, j, jtest, offset, steeri, steerj;
    extern integer errmes_();
    static integer ierror;

/* -----------------------------------------------------------------------
 */
/* PURPOSE */
/*      ASUSMG adds the contribution from an general element matrix into 
*/
/*      a real unsymmetric system matrix, */

/* HISTORY */

/*      Copyright (C) 1996 : CCLRC, Rutherford Appleton Laboratory */
/*                           Chilton, Didcot, Oxfordshire OX11 0QX */

/*      Release 3.0  27 Mar 1983 (CG) */
/*      Release 4.0   2 Oct 1996 (CG) */

/* ARGUMENTS in */
/*      SYSK    contains system matrix prior to addition of */
/*              current element matrix contribution */
/*      ISYSK   first dimension of SYSK (.GE. total number of */
/*              unconstrained degrees of freedom) */
/*      JSYSK   second dimension of SYSK (.GE. HBAND) */
/*      ELK     element matrix */
/*      IELK    first dimension of ELK (.GE. dofel1) */
/*      JELK    second dimension of ELK (.GE. dofel2) */
/*      STR1    contains freedom numbers associated with element */
/*              matrix contributions to system matrix */
/*      ISTR1   first dimension of STR1 (.GE. dofel1) */
/*      STR2    contains freedom numbers associated with element */
/*              matrix contributions to system matrix */
/*      ISTR2   first dimension of STR2 (.GE. dofel2) */
/*      HBAND   semi-bandwidth of system matrix, including */
/*              diagonal */
/*      DOFEL1  maximum degrees of freedom 1 */
/*      DOFEL2  maximum degrees of freedom 2 */
/*      ITEST   error checking option */

/* ARGUMENTS out */
/*      SYSK    system matrix */

/* ROUTINES called */
/*      ERRMES */


/*      SUBROUTINE ASUSMG(SYSK,ISYSK,JSYSK,ELK,IELK,JELK,STR1,ISTR1,STR2, 
*/
/*     *                  ISTR2,HBAND,DOFEL1,DOFEL2,ITEST) */
/* ***********************************************************************
 */



    /* Parameter adjustments */
    --str2;
    --str1;
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
	if (*istr1 < *dofel1 || *istr2 < *dofel2) {
	    ierror = 4;
	}
	if (*ielk < *dofel1 || *jelk < *dofel2) {
	    ierror = 3;
	}
	if (*jsysk < (*hband << 1) - 1) {
	    ierror = 2;
	}
	if (*hband <= 0 || *dofel1 <= 0 || *dofel2 <= 0) {
	    ierror = 1;
	}
	*itest = errmes_(&jtest, &ierror, srname, 6L);
	if (*itest != 0) {
	    return 0;
	}
    }

/*     Assembly loops */

    i__1 = *dofel1;
    for (i = 1; i <= i__1; ++i) {
	steeri = str1[i];
	if (steeri > 0) {
	    i__2 = *dofel2;
	    for (j = 1; j <= i__2; ++j) {
		steerj = str2[j];
		if (steerj > 0) {
		    offset = steerj - steeri + *hband;

/*     Range checking */

		    if (jtest != -1) {
			ierror = 0;
			if (*isysk < steeri || *jsysk < offset) {
			    ierror = 5;
			}
			*itest = errmes_(itest, &ierror, srname, 6L);
			if (*itest != 0) {
			    return 0;
			}
		    }

		    sysk[steeri + offset * sysk_dim1] += elk[i + j * elk_dim1]
			    ;
		}
/* L1000: */
	    }

	}
/* L1010: */
    }

} /* asusmg_ */

