/* matvec.f -- translated by f2c (version of 28 August 1991  0:07:02).
   You must link the resulting object file with the libraries:
	-lF77 -lI77 -lm -lc   (in that order)
*/

#include "f2c.h"


/* Subroutine */ int matvec_(a, ia, ja, v, iv, m, n, w, iw, itest)
doublereal *a;
integer *ia, *ja;
doublereal *v;
integer *iv, *m, *n;
doublereal *w;
integer *iw, *itest;
{
    /* Initialized data */

    static char srname[6+1] = "MATVEC";

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
/*      MATVEC post-multiplies the matrix A by the vector V, storing */
/*      the result in the vector W */

/* HISTORY */

/*      Copyright (C) 1996 : CCLRC, Rutherford Appleton Laboratory */
/*                           Chilton, Didcot, Oxfordshire OX11 0QX */

/*      Release 1.1  29 Oct 1979 (IMS) */
/*      Commented    14 Oct 1980 (KR) */
/*      Release 4.0   2 Oct 1996 (CG) */

/* ARGUMENTS in */
/*      A       array of dimension (IA, JA) */
/*      IA      first dimension of A (.GE. M) */
/*      JA      second dimension of A (.GE. N) */
/*      V       vector of dimension IV */
/*      IV      dimension of V (.GE. N) */
/*      M       number of rows of A to be used in the */
/*              multiplication */
/*      N       number of columns of A and the number of */
/*              elemenets of V to be used in the multiplication */
/*      IW      dimension of vector W (.GE. M) */
/*      ITEST   error checking option */

/* ARGUMENTS out */
/*      W       vector of dimension IW; contains the result of */
/*              the operation W=A*V */

/* ROUTINES called */
/*      ERRMES */


/*      SUBROUTINE MATVEC(A,IA,JA,V,IV,M,N,W,IW,ITEST) */
/* ***********************************************************************
 */

    /* Parameter adjustments */
    --w;
    --v;
    a_dim1 = *ia;
    a_offset = a_dim1 + 1;
    a -= a_offset;

    /* Function Body */

/*     Parameter checking */

    if (*itest != -1) {
	ierror = 0;
	if (*m > *iw) {
	    ierror = 4;
	}
	if (*n > *iv) {
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

    i__1 = *m;
    for (i = 1; i <= i__1; ++i) {
	x = 0.;
	i__2 = *n;
	for (j = 1; j <= i__2; ++j) {
	    x += a[i + j * a_dim1] * v[j];
/* L1000: */
	}
	w[i] = x;
/* L1010: */
    }

} /* matvec_ */

