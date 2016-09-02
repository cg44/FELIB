/* chordn.f -- translated by f2c (version of 28 August 1991  0:07:02).
   You must link the resulting object file with the libraries:
	-lF77 -lI77 -lm -lc   (in that order)
*/

#include "f2c.h"


/* Subroutine */ int chordn_(a, ia, ja, n, hband, itest)
doublereal *a;
integer *ia, *ja, *n, *hband, *itest;
{
    /* Initialized data */

    static char srname[6+1] = "CHORDN";

    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2;

    /* Builtin functions */
    double sqrt();

    /* Local variables */
    static integer i, j, k, l, w;
    static doublereal x;
    static integer jtest, la, lb, ik, lk;
    extern integer errmes_();
    static integer ierror;

/* -----------------------------------------------------------------------
 */
/* PURPOSE */
/*      CHORDN performs choleski reduction on a real symmetric positive */

/*      definite banded matrix.  Only the lower band and */
/*      diagonal are stored in a rectangular array A. */

/* HISTORY */

/*      Copyright (C) 1996 : CCLRC, Rutherford Appleton Laboratory */
/*                           Chilton, Didcot, Oxfordshire OX11 0QX */

/*      Release 1.1  29 Oct 1979 (IMS) */
/*      Commented    12 Feb 1980 (CG) */
/*      Release 4.0   2 Oct 1996 (CG) */

/* ARGUMENTS in */
/*      A       on entry, contains lower band and diagonal of */
/*              pd real symmetric matrix */
/*      IA      first dimension of A (.GE. N) */
/*      JA      second dimension of A (.GE. HBAND) */
/*      N       order of matrix A */
/*      HBAND   semi-bandwidth of A (including diagonal) */
/*      ITEST   error checking option */

/* ARGUMENTS out */
/*      A       reduced matrix L, where the input matrix A has */
/*              been reduced to triangular matrices L and lt */
/*              where A=L lt */

/* ROUTINES called */
/*      ERRMES */

/*     SUBROUTINE CHORDN(A,IA,JA,N,HBAND,ITEST) */
/* ***********************************************************************
 */

    /* Parameter adjustments */
    a_dim1 = *ia;
    a_offset = a_dim1 + 1;
    a -= a_offset;

    /* Function Body */

/*     Parameter checking */

    jtest = *itest;
    if (jtest != -1) {
	ierror = 0;
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

/*     Check on value of A(I,J) */

	if (jtest != -1) {
	    ierror = 0;
	    if (a[i + (w + 1) * a_dim1] - x <= 0.) {
		ierror = 3;
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

} /* chordn_ */

