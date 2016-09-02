/* asusm.f -- translated by f2c (version of 28 August 1991  0:07:02).
   You must link the resulting object file with the libraries:
	-lF77 -lI77 -lm -lc   (in that order)
*/

#include "f2c.h"


/* Subroutine */ int asusm_(sysk, isysk, jsysk, elk, ielk, jelk, steer, 
	isteer, hband, dofel, itest)
doublereal *sysk;
integer *isysk, *jsysk;
doublereal *elk;
integer *ielk, *jelk, *steer, *isteer, *hband, *dofel, *itest;
{
    /* Initialized data */

    static char srname[6+1] = "ASUSM ";

    /* System generated locals */
    integer elk_dim1, elk_offset, sysk_dim1, sysk_offset, i__1, i__2;

    /* Local variables */
    static integer i, j, jtest, offset, steeri, steerj;
    extern integer errmes_();
    static integer ierror;

/* -----------------------------------------------------------------------
 */
/* PURPOSE */
/*      ASUSM assembles the contribution from an element matrix into */
/*      into a real unsymmetric system matrix */

/* HISTORY */

/*      Copyright (C) 1999 : CCLRC, Rutherford Appleton Laboratory */
/*                           Chilton, Didcot, Oxfordshire OX11 0QX */

/*      Release 1.1  29 Oct 1979 (IMS) */
/*      Commented     7 Feb 1980 (KR) */
/*      Release 4.0   2 Oct 1999 (CG) */

/* ARGUMENTS in */
/*      SYSK    contains system matrix into which contributions */
/*              from current element matrix are to be assembled */
/*      ISYSK   first dimension of SYSK (.GE. total number of */
/*              unconstrained degrees of freedom) */
/*      JSYSK   second dimension of SYSK (.GE. 2*HBAND-1) */
/*      ELK     element matrix */
/*      IELK    first dimension of ELK (.GE. DOFEL) */
/*      JELK    second dimension of ELK (.GE. DOFEL) */
/*      STEER   contains freedom numbers associated with element */
/*              matrix contributions to system matrix */
/*      ISTEER  first dimension of STEER (.GE. DOFEL) */
/*      HBAND   semi-bandwidth of system matrix, including */
/*              diagonal */
/*      DOFEL   maximum number of degrees of freedom associated */
/*              with this element type */
/*      ITEST   error checking option */

/* ARGUMENTS out */
/*      SYSK    system matrix */

/* ROUTINES called */
/*      ERRMES */


/*      SUBROUTINE ASUSM(SYSK,ISYSK,JSYSK,ELK,IELK,JELK,STEER,ISTEER, */
/*     *                 HBAND,DOFEL,ITEST) */
/* ***********************************************************************
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
	    ierror = 4;
	}
	if (*ielk < *dofel || *jelk < *dofel) {
	    ierror = 3;
	}
	if (*jsysk < (*hband << 1) - 1) {
	    ierror = 2;
	}
	if (*hband <= 0 || *dofel <= 0) {
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
		    offset = steerj - steeri + *hband;

/*     Range checking on STEERI and STEERJ */

		    if (*itest != -1) {
			ierror = 0;
			if (*isysk < steeri || *jsysk < offset) {
			    ierror = 5;
			}
			*itest = errmes_(&jtest, &ierror, srname, 6L);
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

} /* asusm_ */

