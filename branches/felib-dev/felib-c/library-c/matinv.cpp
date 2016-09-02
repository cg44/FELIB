/* matinv.f -- translated by f2c (version of 28 August 1991  0:07:02).
   You must link the resulting object file with the libraries:
	-lF77 -lI77 -lm -lc   (in that order)
*/

#include "f2c.h"


/* Subroutine */ int matinv_(a, ia, ja, b, ib, jb, n, det, itest)
doublereal *a;
integer *ia, *ja;
doublereal *b;
integer *ib, *jb, *n;
doublereal *det;
integer *itest;
{
    /* Initialized data */

    static char srname[6+1] = "MATINV";

    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset;

    /* Local variables */
    static integer k, l, m;
    static doublereal x;
    extern doublereal unflo_();
    static integer jtest;
    extern integer errmes_();
    static integer ierror;

/* -----------------------------------------------------------------------
 */
/* PURPOSE */
/*      MATINV forms the inverse of matrix A in B */

/* HISTORY */

/*      Copyright (C) 1996 : CCLRC, Rutherford Appleton Laboratory */
/*                           Chilton, Didcot, Oxfordshire OX11 0QX */

/*      Release 1.1  29 Oct 1979 (IMS) */
/*      Commented    12 Oct 1980 (KR) */
/*      Release 4.0   2 Oct 1996 (CG) */

/* ARGUMENTS in */
/*      A       array of dimension (IA, JA) to be inverted */
/*      IA      first dimension of A (.GE. N) */
/*      JA      second dimension of A (.GE. N) */
/*      IB      first dimension of array B (.GE. N) */
/*      JB      second dimension of array B (.GE. N) */
/*      N       order of matrix A (.GT. 1 .AND. .LT. 4) */
/*      ITEST   error checking option */

/* ARGUMENTS out */
/*      B       array of dimension (IB, JB) containing inverse of */
/*              A */
/*      DET     contains value of determinant of matrix A */

/* ROUTINES called */
/*      ERRMES */

/*      SUBROUTINE MATINV(A,IA,JA,B,IB,JB,N,DET,ITEST) */
/* ***********************************************************************
 */

    /* Parameter adjustments */
    b_dim1 = *ib;
    b_offset = b_dim1 + 1;
    b -= b_offset;
    a_dim1 = *ia;
    a_offset = a_dim1 + 1;
    a -= a_offset;

    /* Function Body */

/*     Parameter checking */

    jtest = *itest;
    if (jtest != -1) {
	ierror = 0;
	if (*n > *ib || *n > *jb) {
	    ierror = 3;
	}
	if (*n > *ia || *n > *ja) {
	    ierror = 2;
	}
	if (*n <= 1 || *n >= 4) {
	    ierror = 1;
	}
	*itest = errmes_(&jtest, &ierror, srname, 6L);
	if (*itest != 0) {
	    return 0;
	}
    }

/*     Main body */

    m = *n - 1;
    switch ((int)m) {
	case 1:  goto L1000;
	case 2:  goto L1030;
    }

/*     Code for 2x2 matrix */

L1000:
    *det = a[a_dim1 + 1] * a[(a_dim1 << 1) + 2] - a[(a_dim1 << 1) + 1] * a[
	    a_dim1 + 2];

/*     Check that determinant is not near zero */

    if (jtest != -1) {
	ierror = 0;
	if (abs(*det) < unflo_(&x)) {
	    ierror = 4;
	}
	*itest = errmes_(&jtest, &ierror, srname, 6L);
	if (*itest != 0) {
	    return 0;
	}
    }

    b[b_dim1 + 1] = a[(a_dim1 << 1) + 2];
    b[(b_dim1 << 1) + 1] = -a[(a_dim1 << 1) + 1];
    b[b_dim1 + 2] = -a[a_dim1 + 2];
    b[(b_dim1 << 1) + 2] = a[a_dim1 + 1];
    for (k = 1; k <= 2; ++k) {
	for (l = 1; l <= 2; ++l) {
	    b[k + l * b_dim1] /= *det;
/* L1010: */
	}
/* L1020: */
    }
    return 0;

/*     Code for 3x3 matrix */

L1030:
    *det = a[a_dim1 + 1] * (a[(a_dim1 << 1) + 2] * a[a_dim1 * 3 + 3] - a[(
	    a_dim1 << 1) + 3] * a[a_dim1 * 3 + 2]);
    *det -= a[(a_dim1 << 1) + 1] * (a[a_dim1 + 2] * a[a_dim1 * 3 + 3] - a[
	    a_dim1 + 3] * a[a_dim1 * 3 + 2]);
    *det += a[a_dim1 * 3 + 1] * (a[a_dim1 + 2] * a[(a_dim1 << 1) + 3] - a[
	    a_dim1 + 3] * a[(a_dim1 << 1) + 2]);

/*     Check on determinant not near zero */

    if (jtest != -1) {
	ierror = 0;
	if (abs(*det) < unflo_(&x)) {
	    ierror = 4;
	}
	*itest = errmes_(&jtest, &ierror, srname, 6L);
	if (*itest != 0) {
	    return 0;
	}
    }

    b[b_dim1 + 1] = a[(a_dim1 << 1) + 2] * a[a_dim1 * 3 + 3] - a[(a_dim1 << 1)
	     + 3] * a[a_dim1 * 3 + 2];
    b[b_dim1 + 2] = -a[a_dim1 + 2] * a[a_dim1 * 3 + 3] + a[a_dim1 + 3] * a[
	    a_dim1 * 3 + 2];
    b[b_dim1 + 3] = a[a_dim1 + 2] * a[(a_dim1 << 1) + 3] - a[a_dim1 + 3] * a[(
	    a_dim1 << 1) + 2];
    b[(b_dim1 << 1) + 1] = -a[(a_dim1 << 1) + 1] * a[a_dim1 * 3 + 3] + a[(
	    a_dim1 << 1) + 3] * a[a_dim1 * 3 + 1];
    b[(b_dim1 << 1) + 2] = a[a_dim1 + 1] * a[a_dim1 * 3 + 3] - a[a_dim1 + 3] *
	     a[a_dim1 * 3 + 1];
    b[(b_dim1 << 1) + 3] = -a[a_dim1 + 1] * a[(a_dim1 << 1) + 3] + a[a_dim1 + 
	    3] * a[(a_dim1 << 1) + 1];
    b[b_dim1 * 3 + 1] = a[(a_dim1 << 1) + 1] * a[a_dim1 * 3 + 2] - a[(a_dim1 
	    << 1) + 2] * a[a_dim1 * 3 + 1];
    b[b_dim1 * 3 + 2] = -a[a_dim1 + 1] * a[a_dim1 * 3 + 2] + a[a_dim1 + 2] * 
	    a[a_dim1 * 3 + 1];
    b[b_dim1 * 3 + 3] = a[a_dim1 + 1] * a[(a_dim1 << 1) + 2] - a[a_dim1 + 2] *
	     a[(a_dim1 << 1) + 1];

    for (k = 1; k <= 3; ++k) {
	for (l = 1; l <= 3; ++l) {
	    b[k + l * b_dim1] /= *det;
/* L1040: */
	}
/* L1050: */
    }

} /* matinv_ */

