/* chosol.f -- translated by f2c (version of 28 August 1991  0:07:02).
   You must link the resulting object file with the libraries:
	-lF77 -lI77 -lm -lc   (in that order)
*/

#include "f2c.h"


/* Subroutine */ int chosol_(a, ia, ja, r, ir, n, hband, itest)
doublereal *a;
integer *ia, *ja;
doublereal *r;
integer *ir, *n, *hband, *itest;
{
    /* Initialized data */

    static char srname[6+1] = "CHOSOL";

    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2;

    /* Builtin functions */
    double sqrt();

    /* Local variables */
    static integer i, j, k, l, m, w;
    static doublereal x;
    static integer jtest, la, lb, ij, ik, lk;
    extern integer errmes_();
    static integer ierror;

/* -----------------------------------------------------------------------
 */
/* PURPOSE */
/*      CHOSOL solves a set of real symmetric positive definite banded */
/*      equations with a single right hand side by choleski decomposition.
 */
/*      Only the lower band and diagonal are stored in a rectangular */
/*      array A. */

/* HISTORY */

/*      Copyright (C) 1996 : CCLRC, Rutherford Appleton Laboratory */
/*                           Chilton, Didcot, Oxfordshire OX11 0QX */

/*      Release 1.1  29 Oct 1979 (IMS) */
/*      Commented    12 Feb 1980 (KR) */
/*      Release 4.0   2 Oct 1996 (CG) */

/* ARGUMENTS in */
/*      A       on entry contains lower half of pd symmetric */
/*              band matrix stored as a rectangular array */
/*      IA      first dimension of A (.GE. N) */
/*      JA      second dimension of A (.GE. HBAND) */
/*      R       contains elements of right hand side */
/*      IR      dimension of R (.GE. N) */
/*      N       order of matrix A */
/*      HBAND   semi-bandwidth of A (includes diagonal) */
/*      ITEST   error checking option */

/* ARGUMENTS out */
/*      A       on exit, contains lower triangular reduced */
/*      R       matrix solution vector */

/* ROUTINES called */
/*      ERRMES */


/*      SUBROUTINE CHOSOL(A,IA,JA,R,IR,N,HBAND,ITEST) */
/* ***********************************************************************
 */

    /* Parameter adjustments */
    --r;
    a_dim1 = *ia;
    a_offset = a_dim1 + 1;
    a -= a_offset;

    /* Function Body */

/*     Parameter checking */

    jtest = *itest;
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
	*itest = errmes_(&jtest, &ierror, srname, 6L);
	if (*itest != 0) {
	    return 0;
	}
    }

/*     Main body */

    w = *hband - 1;
    i__1 = *n;
    for (i = 1; i <= i__1; ++i) {
	x = 0.;
	i__2 = w;
	for (j = 1; j <= i__2; ++j) {
	    x += a[i + j * a_dim1] * a[i + j * a_dim1];
/* L1000: */
	}

/*     Range checking on A(I,W+1) */

	if (jtest != -1) {
	    ierror = 0;
	    if (a[i + (w + 1) * a_dim1] - x <= 0.) {
		ierror = 4;
	    }
	    *itest = errmes_(&jtest, &ierror, srname, 6L);
	    if (*itest != 0) {
		return 0;
	    }
	}

	a[i + (w + 1) * a_dim1] = sqrt(a[i + (w + 1) * a_dim1] - x);

	i__2 = w;
	for (k = 1; k <= i__2; ++k) {
	    x = 0.;
	    if (i + k <= *n) {
		if (k != w) {
		    l = w - k;
L1010:
		    ik = i + k;
		    lk = l + k;
		    x += a[ik + l * a_dim1] * a[i + lk * a_dim1];
		    --l;
		    if (l != 0) {
			goto L1010;
		    }
		}
		la = i + k;
		lb = w - k + 1;
		a[la + lb * a_dim1] = (a[la + lb * a_dim1] - x) / a[i + (w + 
			1) * a_dim1];
	    }
/* L1020: */
	}
/* L1030: */
    }
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
/* L1040: */
	}
	r[i] = (r[i] - x) / a[i + (w + 1) * a_dim1];
/* L1050: */
    }
    r[*n] /= a[*n + (w + 1) * a_dim1];
    i = *n - 1;
L1060:
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
/* L1070: */
    }
    r[i] = (r[i] - x) / a[i + (w + 1) * a_dim1];
    --i;
    if (i != 0) {
	goto L1060;
    }

} /* chosol_ */

