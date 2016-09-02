/* qlvec.f -- translated by f2c (version of 28 August 1991  0:07:02).
   You must link the resulting object file with the libraries:
	-lF77 -lI77 -lm -lc   (in that order)
*/

#include "f2c.h"


/* Subroutine */ int qlvec_(diag, idiag, sub, isub, t, it, jt, n, eps, itest)
doublereal *diag;
integer *idiag;
doublereal *sub;
integer *isub;
doublereal *t;
integer *it, *jt, *n;
doublereal *eps;
integer *itest;
{
    /* Initialized data */

    static char srname[6+1] = "QLVEC ";

    /* System generated locals */
    integer t_dim1, t_offset, i__1, i__2, i__3;
    doublereal d__1, d__2;

    /* Builtin functions */
    double sqrt();

    /* Local variables */
    static doublereal b, c, f, g, h;
    static integer i, j, k, l, m;
    static doublereal p, r, s;
    static integer jtest, i1, m1, ii;
    extern integer errmes_();
    static integer ierror;

/* -----------------------------------------------------------------------
 */
/* PURPOSE */
/*      QLVEC calculates all the eigenvalues and eigenvectors of a */
/*      real symmetric tridiagonal matrix, or of a full real */
/*      symmetric matrix reduced to tridiagonal form by HOUSE */

/* HISTORY */

/*      Copyright (C) 1996 : CCLRC, Rutherford Appleton Laboratory */
/*                           Chilton, Didcot, Oxfordshire OX11 0QX */

/*      Release 1.1  29 Oct 1979 (ims) */
/*      Commented    16 Oct 1980 (KR) */
/*      Release 4.0   2 Oct 1996 (CG) */

/* ARGUMENTS in */
/*      DIAG    vector of length IDIAG. On entry, contains */
/*              diagonal elements of tridiagonal matrix */
/*      IDIAG   dimension of vector DIAG (.GE.N) */
/*      SUB     vector of length ISUB. On entry, contains */
/*              SUB-DIAGONAL elements of tridiagonal matrix in */
/*              SUB(2) to SUB(N). On exit, contents of SUB are */
/*              destroyed */
/*      ISUB    dimension of vector SUB (.GE. N) */
/*      T       array of dimension (IT, JT).  There are two modes */
/*              of operation: */
/*                (I)  eigen vectors of tridiagonal matrix */
/*                     on entry T should contain the identity */
/*                     matrix of order N */
/*                (II) eigenvectors of full symmetric matrix */
/*                     on entry T should contain the orthogonal */
/*                     matrix obtained from HOUSE */
/*      IT      first dimension of array T (.GE. N) */
/*      JT      second dimension of T (.GE. N) */
/*      N       order of tridiagonal matrix */
/*      EPS     smallest positive number such that 1.+EPS>1. */
/*      ITEST   error checking option */

/* ARGUMENTS out */
/*      DIAG    vector of dimension IDIAG.  on exit, contains */
/*              the eigenvalues in ascending order */
/*      T       array of dimension (IT, JT).  on exit, T contains */
/*              the normalised eigenvectors , with T(I, J), */
/*              I=1(1)N containing the eigenvector corresponding */
/*              to the eigenvalue in DIAG(J) */

/* ROUTINES called */
/*      ERRMES */

/*      SUBROUTINE QLVEC(DIAG,IDIAG,SUB,ISUB,T,IT,JT,N,EPS,ITEST) */
/* ***********************************************************************
 */

    /* Parameter adjustments */
    t_dim1 = *it;
    t_offset = t_dim1 + 1;
    t -= t_offset;
    --sub;
    --diag;

    /* Function Body */

/*     Parameter checking */

    jtest = *itest;
    if (jtest != -1) {
	ierror = 0;
	if (*it < *n || *jt < *n) {
	    ierror = 3;
	}
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
		goto L1110;
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

/*     Form vector */

		    diag[i + 1] = h + s * (c * g + s * diag[i]);
		    i__3 = *n;
		    for (k = 1; k <= i__3; ++k) {
			h = t[k + (i + 1) * t_dim1];
			t[k + (i + 1) * t_dim1] = s * t[k + i * t_dim1] + c * 
				h;
			t[k + i * t_dim1] = c * t[k + i * t_dim1] - s * h;
/* L1050: */
		    }
/* L1060: */
		}

		sub[l] = s * p;
		diag[l] = c * p;
		if ((d__1 = sub[l], abs(d__1)) > b) {
		    goto L1030;
		}
	    }
	}
	diag[l] += f;
/* L1070: */
    }

/*     Order eigenvalues and eigenvectors */

    i__1 = *n;
    for (i = 1; i <= i__1; ++i) {
	k = i;
	p = diag[i];
	i1 = i + 1;
	if (i1 <= *n) {
	    i__2 = *n;
	    for (j = i1; j <= i__2; ++j) {
		if (diag[j] < p) {
		    k = j;
		    p = diag[j];
		}
/* L1080: */
	    }
	}

	if (k != i) {
	    diag[k] = diag[i];
	    diag[i] = p;
	    i__2 = *n;
	    for (j = 1; j <= i__2; ++j) {
		p = t[j + i * t_dim1];
		t[j + i * t_dim1] = t[j + k * t_dim1];
		t[j + k * t_dim1] = p;
/* L1090: */
	    }
	}
/* L1100: */
    }

    *itest = 0;
    return 0;

L1110:
    if (jtest != -1) {
	ierror = 1;
	*itest = errmes_(&jtest, &ierror, srname, 6L);
    }

} /* qlvec_ */

