/* chosub.f -- translated by f2c (version of 28 August 1991  0:07:02).
   You must link the resulting object file with the libraries:
	-lF77 -lI77 -lm -lc   (in that order)
*/

#include "f2c.h"


/* Subroutine */ int chosub_(a, ia, ja, r, ir, n, hband, itest)
doublereal *a;
integer *ia, *ja;
doublereal *r;
integer *ir, *n, *hband, *itest;
{
    /* Initialized data */

    static char srname[6+1] = "CHOSUB";

    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2;

    /* Local variables */
    static integer i, j, k, l, m, w;
    static doublereal x;
    static integer ij;
    extern integer errmes_();
    static integer ierror;

/* -----------------------------------------------------------------------
 */
/* PURPOSE */
/*      CHOSUB performs forward and backward substitution on a matrix */
/*      reduced by CHORDN */

/* HISTORY */

/*      Copyright (C) 1996 : CCLRC, Rutherford Appleton Laboratory */
/*                           Chilton, Didcot, Oxfordshire OX11 0QX */

/*      Release 1.1  29 Oct 1979 (CG) */
/*      Commented    10 Oct 1980 (KR) */
/*      Release 4.0   2 Oct 1996 (CG) */

/* ARGUMENTS in */
/*      A       array of dimension (IA,JA).  contains the */
/*              elements of the lower half of the positive */
/*              definite symmetric band matrix of order N and */
/*              semi-bandwidth HBAND.  A should previously have */
/*              been reduced using CHORDN or CHOSOL */
/*      IA      first dimension of A (.GE.N) */
/*      JA      second dimension of A (.GE.HBAND) */
/*      R       on entry, contains the vector of rhs's */
/*      IR      dimension of R (.GE.N) */
/*      N       order of matrix A */
/*      HBAND   semi-bandwidth of A */
/*      ITEST   error checking option */

/* ARGUMENTS out */
/*      R       on exit, contains the solution vector */

/* ROUTINES called */
/*      ERRMES */


/*      SUBROUTINE CHOSUB(A,IA,JA,R,IR,N,HBAND,ITEST) */
/* ***********************************************************************
 */

    /* Parameter adjustments */
    --r;
    a_dim1 = *ia;
    a_offset = a_dim1 + 1;
    a -= a_offset;

    /* Function Body */

/*     Parameter checking */

    if (*itest != -1) {
	ierror = 0;
	if (*ir < *n) {
	    ierror = 3;
	}
	if (*ia < *n || *ja < *hband) {
	    ierror = 2;
	}
	if (*n <= 0 || *hband <= 0) {
	    ierror = 1;
	}
	*itest = errmes_(itest, &ierror, srname, 6L);
	if (*itest != 0) {
	    return 0;
	}
    }

/*     Main body */

    w = *hband - 1;
    r[1] /= a[(w + 1) * a_dim1 + 1];
    i__1 = *n;
    for (i = 2; i <= i__1; ++i) {
	x = 0.;
	k = 1;
	if (i <= w + 1) {
	    k = w - i + 2;
	}
	i__2 = w;
	for (j = k; j <= i__2; ++j) {
	    ij = i + j - w - 1;
	    x += a[i + j * a_dim1] * r[ij];
/* L1000: */
	}
	r[i] = (r[i] - x) / a[i + (w + 1) * a_dim1];
/* L1010: */
    }
    r[*n] /= a[*n + (w + 1) * a_dim1];
    i = *n - 1;
L1020:
    x = 0.;
    l = i + w;
    if (i > *n - w) {
	l = *n;
    }
    m = i + 1;
    i__1 = l;
    for (j = m; j <= i__1; ++j) {
	ij = w + i - j + 1;
	x += a[j + ij * a_dim1] * r[j];
/* L1030: */
    }
    r[i] = (r[i] - x) / a[i + (w + 1) * a_dim1];
    --i;
    if (i != 0) {
	goto L1020;
    }

} /* chosub_ */

