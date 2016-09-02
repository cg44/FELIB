/* vecmat.f -- translated by f2c (version of 28 August 1991  0:07:02).
   You must link the resulting object file with the libraries:
	-lF77 -lI77 -lm -lc   (in that order)
*/

#include "f2c.h"


/* Subroutine */ int vecmat_(v, iv, a, ia, ja, m, n, w, iw, itest)
doublereal *v;
integer *iv;
doublereal *a;
integer *ia, *ja, *m, *n;
doublereal *w;
integer *iw, *itest;
{
    /* Initialized data */

    static char srname[6+1] = "VECMAT";

    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2;

    /* Local variables */
    static integer i, j;
    static doublereal x;
    extern integer errmes_();
    static integer ierror;

/* -----------------------------------------------------------------------
 */
/* PURPOSE */
/*      VECMAT pre-multiplies the matrix A by the vector V, storing */
/*      the result in the vector W */

/* HISTORY */

/*      Copyright (C) 1996 : CCLRC, Rutherford Appleton Laboratory */
/*                           Chilton, Didcot, Oxfordshire OX11 0QX */

/*      Release 2.0   1 Feb 1981 (CG) */
/*      Commented     1 Feb 1981 (CG) */
/*      Release 4.0   2 Oct 1996 (CG) */

/* ARGUMENTS in */
/*      V       vector of dimension IV */
/*      IV      dimension of V (.GE. M) */
/*      A       array of dimension (IA, JA) */
/*      IA      first dimension of A (.GE. M) */
/*      JA      second dimension of A (.GE. N) */
/*      M       number of rows of A to be used in the */
/*              multiplication */
/*      N       number of columns of A and the number of */
/*              elemenets of V to be used in the multiplication */
/*      IW      dimension of vector W (.GE. N) */
/*      ITEST   error checking option */

/* ARGUMENTS out */
/*      W       vector of dimension IW; contains the result of */
/*              the operation W=V*A */

/* ROUTINES called */
/*      ERRMES */

/*     SUBROUTINE VECMAT(V,IV,A,IA,JA,M,N,W,IW,ITEST) */
/* ***********************************************************************
 */


    /* Parameter adjustments */
    --w;
    a_dim1 = *ia;
    a_offset = a_dim1 + 1;
    a -= a_offset;
    --v;

    /* Function Body */

/*     Parameter checking */

    if (*itest != -1) {
	ierror = 0;
	if (*n > *iw) {
	    ierror = 4;
	}
	if (*m > *iv) {
	    ierror = 3;
	}
	if (*m > *ia || *n > *ja) {
	    ierror = 2;
	}
	if (*m <= 0 || *n <= 0) {
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
	x = 0.;
	i__2 = *m;
	for (j = 1; j <= i__2; ++j) {
	    x += a[j + i * a_dim1] * v[j];
/* L1000: */
	}
	w[i] = x;
/* L1010: */
    }

} /* vecmat_ */

