/* dyad.f -- translated by f2c (version of 28 August 1991  0:07:02).
   You must link the resulting object file with the libraries:
	-lF77 -lI77 -lm -lc   (in that order)
*/

#include "f2c.h"


/* Subroutine */ int dyad_(v, iv, w, iw, a, ia, ja, n, itest)
doublereal *v;
integer *iv;
doublereal *w;
integer *iw;
doublereal *a;
integer *ia, *ja, *n, *itest;
{
    /* Initialized data */

    static char srname[6+1] = "DYAD  ";

    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2;

    /* Local variables */
    static integer i, j;
    extern integer errmes_();
    static integer ierror;

/* -----------------------------------------------------------------------
 */
/* PURPOSE */
/*      DYAD forms the dyad matrix A from two vectors V and W */

/* HISTORY */

/*      Copyright (C) 1996 : CCLRC, Rutherford Appleton Laboratory */
/*                           Chilton, Didcot, Oxfordshire OX11 0QX */

/*      Release 1.1  29 Oct 1979 (CG) */
/*      Commented    14 Feb 1980 (KR) */
/*      Release 4.0   2 Oct 1996 (CG) */

/* ARGUMENTS in */
/*      V       first vector in multipcation */
/*      IV      dimension of vector V (.GE. N) */
/*      W       second vector in multiplication */
/*      IW      dimension of vector W (.GE. N) */
/*      IA      first dimension of dyadic array A (.GE. N) */
/*      JA      second dimension of A (.GE. N) */
/*      N       number of elements of vectors V and W to be used */
/*              in forming the dyad */
/*      ITEST   error checking option */

/* ARGUMENTS out */
/*      A       dyadic array formed by product of V and W */

/* ROUTINES called */
/*      ERRMES */

/*      SUBROUTINE DYAD(V,IV,W,IW,A,IA,JA,N,ITEST) */
/* ***********************************************************************
 */

    /* Parameter adjustments */
    a_dim1 = *ia;
    a_offset = a_dim1 + 1;
    a -= a_offset;
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
	if (*n > *ia || *n > *ja) {
	    ierror = 4;
	}
	*itest = errmes_(itest, &ierror, srname, 6L);
	if (*itest != 0) {
	    return 0;
	}
    }

/*     Main body */

    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	i__2 = *n;
	for (i = 1; i <= i__2; ++i) {
	    a[i + j * a_dim1] = v[i] * w[j];
/* L1000: */
	}
/* L1010: */
    }

} /* dyad_ */

