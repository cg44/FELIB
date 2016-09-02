/* qlval.f -- translated by f2c (version of 28 August 1991  0:07:02).
   You must link the resulting object file with the libraries:
	-lF77 -lI77 -lm -lc   (in that order)
*/

#include "f2c.h"


/* Subroutine */ int qlval_(diag, idiag, sub, isub, n, eps, itest)
doublereal *diag;
integer *idiag;
doublereal *sub;
integer *isub, *n;
doublereal *eps;
integer *itest;
{
    /* Initialized data */

    static char srname[6+1] = "QLVAL ";

    /* System generated locals */
    integer i__1, i__2;
    doublereal d__1, d__2;

    /* Builtin functions */
    double sqrt();

    /* Local variables */
    static doublereal b, c, f, g, h;
    static integer i, j, l, m;
    static doublereal p, r, s;
    static integer jtest, i1, m1, ii;
    extern integer errmes_();
    static integer ierror;

/* -----------------------------------------------------------------------
 */
/* PURPOSE */
/*      QLVAL calculates all the eigenvalues of a real symmetric */
/*      tridiagonal matrix or of a full real symmetric matrix */
/*      that has been reduced to tridiagonal form using HOUSE */

/* HISTORY */

/*      Copyright (C) 1996 : CCLRC, Rutherford Appleton Laboratory */
/*                           Chilton, Didcot, Oxfordshire OX11 0QX */

/*      Release 1.1  29 Oct 1979 (IMS) */
/*      Commented    15 Oct 1980 (KR) */
/*      Release 4.0   2 Oct 1996 (CG) */

/* ARGUMENTS in */
/*      DIAG    vector of length IDIAG.  On entry, contains the */
/*              diagonal elements of the tridiagonal matrix */
/*      IDIAG   dimension of vector DIAG (.GE. N) */
/*      SUB     vector of dimension ISUB.  On entry, contains */
/*              SUB-DIAGONAL elements of tridiagonal matrix in */
/*              elements SUB(2) to SUB(N). SUB(1) May be */
/*              arbitrary. Contents destroyed during execution */
/*              of HOUSE */
/*      ISUB    dimension of SUB (.GE. N) */
/*      N       order of tridiagonal matrix */
/*      EPS     smallest positive number such that 1.+EPS>1. */
/*      ITEST   error checking option */

/* ARGUMENTS out */
/*      DIAG    vector of dimension IDIAG.  on exit, contains */
/*              eigenvalues in ascending order */

/* ROUTINES called */
/*      ERRMES */

/*      SUBROUTINE QLVAL(DIAG,IDIAG,SUB,ISUB,N,EPS,ITEST) */
/* ***********************************************************************
 */

    /* Parameter adjustments */
    --sub;
    --diag;

    /* Function Body */

/*     Parameter checking */

    jtest = *itest;
    if (jtest != -1) {
	ierror = 0;
	if (*idiag < *n || *isub < *n) {
	    ierror = 2;
	}
	if (*n <= 0) {
	    ierror = 1;
	}
	*itest = errmes_(&jtest, &ierror, srname, 6L);
	if (*itest != 0) {
	    return 0;
	}
    }

/*     Main body */

    if (*n != 1) {
	i__1 = *n;
	for (i = 2; i <= i__1; ++i) {
	    sub[i - 1] = sub[i];
/* L1000: */
	}
    }

    sub[*n] = 0.;
    b = 0.;
    f = 0.;

    i__1 = *n;
    for (l = 1; l <= i__1; ++l) {
	j = 0;
	h = *eps * ((d__1 = diag[l], abs(d__1)) + (d__2 = sub[l], abs(d__2)));


/*     Look for small SUB-DIAGONAL element */

	if (b < h) {
	    b = h;
	}
	i__2 = *n;
	for (m = l; m <= i__2; ++m) {
	    if ((d__1 = sub[m], abs(d__1)) <= b) {
		goto L1020;
	    }
/* L1010: */
	}

L1020:
	if (m != l) {
L1030:
	    if (j == 30) {
		goto L1090;
	    } else {

/*     Form shift */

		++j;
		g = diag[l];
		h = diag[l + 1] - g;
		if (abs(h) >= (d__1 = sub[l], abs(d__1))) {

		    p = sub[l] * 2. / h;
		    r = sqrt(p * p + 1.);
		    diag[l] = sub[l] * p / (r + 1.);
		} else {
		    p = h * .5 / sub[l];
		    r = sqrt(p * p + 1.);
		    h = p + r;
		    if (p < 0.) {
			h = p - r;
		    }
		    diag[l] = sub[l] / h;
		}
		h = g - diag[l];
		i1 = l + 1;
		if (i1 <= *n) {

		    i__2 = *n;
		    for (i = i1; i <= i__2; ++i) {
			diag[i] -= h;
/* L1040: */
		    }
		}

/*     Ql transformation */

		f += h;
		p = diag[m];
		c = 1.;
		s = 0.;
		m1 = m - 1;
		i__2 = m1;
		for (ii = l; ii <= i__2; ++ii) {
		    i = m1 - ii + l;
		    g = c * sub[i];
		    h = c * p;
		    if (abs(p) < (d__1 = sub[i], abs(d__1))) {
			c = p / sub[i];
			r = sqrt(c * c + 1.);
			sub[i + 1] = s * sub[i] * r;
			s = 1. / r;
			c /= r;
		    } else {
			c = sub[i] / p;
			r = sqrt(c * c + 1.);
			sub[i + 1] = s * p * r;
			s = c / r;
			c = 1. / r;
		    }
		    p = c * diag[i] - s * g;
		    diag[i + 1] = h + s * (c * g + s * diag[i]);
/* L1050: */
		}

		sub[l] = s * p;
		diag[l] = c * p;
		if ((d__1 = sub[l], abs(d__1)) > b) {
		    goto L1030;
		}
	    }
	}

/*     Order eigenvalue */

	p = diag[l] + f;
	if (l != 1) {
	    i__2 = l;
	    for (ii = 2; ii <= i__2; ++ii) {
		i = l - ii + 2;
		if (p >= diag[i - 1]) {
		    goto L1070;
		} else {
		    diag[i] = diag[i - 1];
		}
/* L1060: */
	    }
	}
	i = 1;
L1070:
	diag[i] = p;
/* L1080: */
    }

    *itest = 0;
    return 0;

L1090:
    if (jtest != -1) {
	ierror = 3;
	*itest = errmes_(&jtest, &ierror, srname, 6L);

    }
} /* qlval_ */

