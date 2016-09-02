/* iassym.f -- translated by f2c (version of 28 August 1991  0:07:02).
   You must link the resulting object file with the libraries:
	-lF77 -lI77 -lm -lc   (in that order)
*/

#include "f2c.h"


/* Subroutine */ int iassym_(sysk, isysk, jsysk, elk, ielk, jelk, steer, 
	isteer, hband, dofel, itest)
doublereal *sysk;
integer *isysk, *jsysk;
doublereal *elk;
integer *ielk, *jelk, *steer, *isteer, *hband, *dofel, *itest;
{
    /* Initialized data */

    static char srname[6+1] = "IASSYM";

    /* System generated locals */
    integer elk_dim1, elk_offset, sysk_dim2, sysk_offset, i__1, i__2;

    /* Local variables */
    static integer i, j, jtest, cd, steeri;
    extern integer errmes_();
    static integer ierror;

/* -----------------------------------------------------------------------
 */
/* PURPOSE */
/*      IASSYM for real symmetric system matrix, adds the contribution */
/*      from an element matrix into real part of complex SYSK */

/* HISTORY */

/*      Copyright (C) 1996 : CCLRC, Rutherford Appleton Laboratory */
/*                           Chilton, Didcot, Oxfordshire OX11 0QX */

/*      Release 3.0  29 Oct 1985 (CRIE) */
/*      Commented     7 Nov 1985 (CG) */
/*      Release 4.0   2 Oct 1996 (CG) */

/* ARGUMENTS in */
/*      SYSK    contains system matrix prior to addition of */
/*              current element matrix contribution */
/*      ISYSK   first dimension of SYSK (.GE. total number of */
/*              unconstrained degrees of freedom) */
/*      JSYSK   second dimension of SYSK (.GE. HBAND) */
/*      ELK     element matrix */
/*      IELK    first dimension of ELK (.GE. DOFEL) */
/*      JELK    second dimension of ELK (.GE. DOFEL) */
/*      STEER   contains freedom numbers associated with element */
/*              matrix contributions to system matrix */
/*      ISTEER  first dimension of STEER (.GE. DOFEL) */
/*      HBAND   semi-bandwidth of system matrix, including */
/*              diagonal */
/*      DOFEL   maximum degrees of freedom associated with */
/*              element type */
/*      ITEST   error checking option */

/* ARGUMENTS out */
/*      SYSK    system matrix -   ordered pairs */

/* ROUTINES called */
/*      ERRMES */

/*      SUBROUTINE IASSYM(SYSK,ISYSK,JSYSK,ELK,IELK,JELK,STEER,ISTEER, */
/*     *                  HBAND,DOFEL,ITEST) */
/* ***********************************************************************
 */



    /* Parameter adjustments */
    --steer;
    elk_dim1 = *ielk;
    elk_offset = elk_dim1 + 1;
    elk -= elk_offset;
    sysk_dim2 = *isysk;
    sysk_offset = (sysk_dim2 + 1 << 1) + 1;
    sysk -= sysk_offset;

    /* Function Body */

/*     Parameter checking */

    jtest = *itest;
    if (*itest != -1) {
	ierror = 0;
	if (*isteer < *dofel) {
	    ierror = 4;
	}
	if (*ielk < *dofel || *jelk < *dofel) {
	    ierror = 3;
	}
	if (*jsysk < *hband) {
	    ierror = 2;
	}
	if (*hband <= 0 || *dofel <= 0) {
	    ierror = 1;
	}
	*itest = errmes_(itest, &ierror, srname, 6L);
	if (*itest != 0) {
	    return 0;
	}
    }
    i__1 = *dofel;
    for (i = 1; i <= i__1; ++i) {
	if (steer[i] != 0) {
	    i__2 = *dofel;
	    for (j = 1; j <= i__2; ++j) {
		if (steer[j] != 0) {
		    cd = steer[j] - steer[i] + *hband;
		    if (cd <= *hband) {
			steeri = steer[i];
			if (jtest != -1) {
			    ierror = 0;
			    if (*isysk < steeri) {
				ierror = 5;
			    }
			    *itest = errmes_(itest, &ierror, srname, 6L);
			    if (*itest != 0) {
				return 0;
			    }
			}
			sysk[(steeri + cd * sysk_dim2 << 1) + 2] += elk[i + j 
				* elk_dim1];
		    }
		}
/* L1000: */
	    }
	}
/* L1010: */
    }

} /* iassym_ */

