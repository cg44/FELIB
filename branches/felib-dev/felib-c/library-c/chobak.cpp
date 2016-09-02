/* chobak.f -- translated by f2c (version of 28 August 1991  0:07:02).
   You must link the resulting object file with the libraries:
	-lF77 -lI77 -lm -lc   (in that order)
*/

#include "f2c.h"


/* Subroutine */ int chobak_(a, ia, ja, r, ir, n, hband, itest)
doublereal *a;
integer *ia, *ja;
doublereal *r;
integer *ir, *n, *hband, *itest;
{
    /* Initialized data */

    static char srname[6+1] = "CHOBAK";

    /* System generated locals */
    integer a_dim1, a_offset, i__1;

    /* Local variables */
    static integer i, j, l, m, w;
    static doublereal x;
    static integer ij;
    extern integer errmes_();
    static integer ierror;

/* -----------------------------------------------------------------------
 */
/* PURPOSE */
/*      CHOBAK performs bakward substitution on a matrix processed by */
/*      CHORDN and CHOFWD */

/* HISTORY */

/*      Copyright (C) 1996 : CCLRC, Rutherford Appleton Laboratory */
/*                           Chilton, Didcot, Oxfordshire OX11 0QX */

/*      Release 1.1  29 Oct 1979 (CG) */
/*      Commented    10 Oct 1980 (KR) */
/*      Release 4.0   2 Oct 1996 (CG) */

/* ARGUMENTS in */
/*      A       array of dimension (IA, JA).  contains the */
/*              elements of the lower half of the positive */
/*              definite band matrix of order N and with semi- */
/*              bandwidth HBAND, reduced by CHORDN or CHOSOL */
/*      IA      first dimension of A (.GE. N) */
/*      JA      second dimension of A (.GE. HBAND) */
/*      R       on entry, contains elements of N rhs's after */
/*              processing by CHOFWD */
/*      IR      dimension of vector R (.GE. N) */
/*      N       order of matrix A */
/*      HBAND   semi-bandwidth of matrix A */
/*      ITEST   error checking option */

/* ARGUMENTS out */
/*      R       on exit, contains solution vector */

/* ROUTINES called */
/*      ERRMES */


/*      SUBROUTINE CHOBAK(A,IA,JA,R,IR,N,HBAND,ITEST) */
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
    r[*n] /= a[*n + (w + 1) * a_dim1];
    i = *n - 1;
L1000:
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
/* L1010: */
    }
    r[i] = (r[i] - x) / a[i + (w + 1) * a_dim1];
    --i;
    if (i != 0) {
	goto L1000;
    }

} /* chobak_ */

