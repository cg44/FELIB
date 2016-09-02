/* vecadd.f -- translated by f2c (version of 28 August 1991  0:07:02).
   You must link the resulting object file with the libraries:
	-lF77 -lI77 -lm -lc   (in that order)
*/

#include "f2c.h"


/* Subroutine */ int vecadd_(v, iv, w, iw, n, itest)
doublereal *v;
integer *iv;
doublereal *w;
integer *iw, *n, *itest;
{
    /* Initialized data */

    static char srname[6+1] = "VECADD";

    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i;
    extern integer errmes_();
    static integer ierror;

/* -----------------------------------------------------------------------
 */
/* PURPOSE */
/*      VECADD adds the vector V to the vector W, storing the result */
/*      in V */

/* HISTORY */

/*      Copyright (C) 1996 : CCLRC, Rutherford Appleton Laboratory */
/*                           Chilton, Didcot, Oxfordshire OX11 0QX */

/*      Release 1.1  29 Oct 1979 (IMS) */
/*      Commented    23 Oct 1980 (KR) */
/*      Release 4.0   2 Oct 1996 (CG) */

/* ARGUMENTS in */
/*      V       vector of length IV. On entry, contains the */
/*              values of one of the vectors to be added */
/*              in elements 1 to N */
/*      IV      length of vector V (.GE. N) */
/*      W       vector of length IW. Contains the values of the */
/*              second vector to be added in elements 1 to N */
/*      IW      length of vector W (.GE. N) */
/*      N       number of elements of V and W to be added */
/*      ITEST   error checking option */

/* ARGUMENTS out */
/*      V       vector of length IV. On exit, V(I) contains */
/*              the sum of V(I) and W(I) for I=1(1)N */

/* ROUTINES called */
/*      ERRMES */

/*     SUBROUTINE VECADD(V,IV,W,IW,N,ITEST) */
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
	if (*itest != 0) {
	    return 0;
	}
    }

/*     Main body */

    i__1 = *n;
    for (i = 1; i <= i__1; ++i) {
	v[i] += w[i];
/* L1000: */
    }

} /* vecadd_ */

