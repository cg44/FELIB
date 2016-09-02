/* scaprd.f -- translated by f2c (version of 28 August 1991  0:07:02).
   You must link the resulting object file with the libraries:
	-lF77 -lI77 -lm -lc   (in that order)
*/

#include "f2c.h"


/* Subroutine */ int scaprd_(v, iv, w, iw, n, prodct, itest)
doublereal *v;
integer *iv;
doublereal *w;
integer *iw, *n;
doublereal *prodct;
integer *itest;
{
    /* Initialized data */

    static char srname[6+1] = "SCAPRD";

    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i;
    extern integer errmes_();
    static integer ierror;

/* -----------------------------------------------------------------------
 */
/* PURPOSE */
/*      SCAPR forms the scalar product of two vectors V and W, storing */
/*      the result in PRODCT */

/* HISTORY */

/*      Copyright (C) 1996 : CCLRC, Rutherford Appleton Laboratory */
/*                           Chilton, Didcot, Oxfordshire OX11 0QX */

/*      Release 1.1  29 Oct 1979 (CG) */
/*      Commented    22 Oct 1980 (KR) */
/*      Release 4.0   2 Oct 1996 (CG) */

/* ARGUMENTS in */
/*      V       vector of length IV */
/*      IV      dimension of V (.GE. N) */
/*      W       vector of length IW */
/*      IW      dimension of W (.GE. N) */
/*      N       number of elements of V (and W) to be used in */
/*              forming scalar product */
/*      ITEST   error checking option */

/* ARGUMENTS out */
/*      PRODCT  contains: sigma I = 1 to N of (V(I)*W(I)) */

/* ROUTINES called */
/*      ERRMES */

/*     SUBROUTINE SCAPRD(V,IV,W,IW,N,PRODCT,ITEST) */
/* ***********************************************************************
 */

    /* Parameter adjustments */
    --w;
    --v;

    /* Function Body */

/*     Parameter checking */

    if (*itest != -1) {
	ierror = 0;
	if (*n > *iw) {
	    ierror = 3;
	}
	if (*n > *iv) {
	    ierror = 2;
	}
	if (*n <= 0) {
	    ierror = 1;
	}
	*itest = errmes_(itest, &ierror, srname, 6L);
	if (ierror != 0) {
	    return 0;
	}
    }

/*     Main body */

    *prodct = 0.;
    i__1 = *n;
    for (i = 1; i <= i__1; ++i) {
	*prodct += v[i] * w[i];
/* L1000: */
    }

} /* scaprd_ */

