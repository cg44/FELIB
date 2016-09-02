/* vecnul.f -- translated by f2c (version of 28 August 1991  0:07:02).
   You must link the resulting object file with the libraries:
	-lF77 -lI77 -lm -lc   (in that order)
*/

#include "f2c.h"


/* Subroutine */ int vecnul_(vec, ivec, n, itest)
doublereal *vec;
integer *ivec, *n, *itest;
{
    /* Initialized data */

    static char srname[6+1] = "VECNUL";

    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i;
    extern integer errmes_();
    static integer ierror;

/* -----------------------------------------------------------------------
 */
/* PURPOSE */
/*      VECNUL sets the first N elements of a vector to zero. */

/* HISTORY */

/*      Copyright (C) 1996 : CCLRC, Rutherford Appleton Laboratory */
/*                           Chilton, Didcot, Oxfordshire OX11 0QX */

/*      Release 1.1  29 Oct 1979 (IMS) */
/*      Commented    23 Oct 1980 (KR) */
/*      Release 4.0   2 Oct 1996 (CG) */

/* ARGUMENTS in */
/*      IVEC    length of vector VEC (.GE. N) */
/*      N       number of elements to be set to zero */
/*      ITEST   error checking option */

/* ARGUMENTS out */
/*      VEC     vector of length IVEC. VEC(I)=0.0d0 for I=1(1)N */

/* ROUTINES called */
/*      ERRMES */

/*      SUBROUTINE VECNUL(VEC,IVEC,N,ITEST) */
/* ***********************************************************************
 */


    /* Parameter adjustments */
    --vec;

    /* Function Body */

/*     Parameter checking */

    if (*itest != -1) {
	ierror = 0;
	if (*n > *ivec) {
	    ierror = 2;
	}
	if (*n <= 0) {
	    ierror = 1;
	}
	*itest = errmes_(itest, &ierror, srname, 6L);
	if (*itest != 0) {
	    return 0;
	}
    }

/*     Main body */

    i__1 = *n;
    for (i = 1; i <= i__1; ++i) {
	vec[i] = 0.;
/* L1000: */
    }

} /* vecnul_ */

