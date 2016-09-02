/* house.f -- translated by f2c (version of 28 August 1991  0:07:02).
   You must link the resulting object file with the libraries:
	-lF77 -lI77 -lm -lc   (in that order)
*/

#include "f2c.h"


/* Subroutine */ int house_(a, ia, ja, t, it, jt, diag, idiag, sub, isub, n, 
	tol, itest)
doublereal *a;
integer *ia, *ja;
doublereal *t;
integer *it, *jt;
doublereal *diag;
integer *idiag;
doublereal *sub;
integer *isub, *n;
doublereal *tol;
integer *itest;
{
    /* Initialized data */

    static char srname[6+1] = "HOUSE ";

    /* System generated locals */
    integer a_dim1, a_offset, t_dim1, t_offset, i__1, i__2, i__3;
    doublereal d__1;

    /* Builtin functions */
    double sqrt();

    /* Local variables */
    extern doublereal veps_();
    static doublereal f, g, h;
    static integer i, j, k, l;
    static doublereal x;
    static integer j1;
    static doublereal hh;
    static integer ii;
    extern integer errmes_();
    static integer ierror;

/* -----------------------------------------------------------------------
 */
/* PURPOSE */
/*      HOUSE uses householder's method to reduce a real symmetric */
/*      matrix to tridiagonal form for use with QLVAL, QLVEC */
/*      A is stored as A full matrix */

/* HISTORY */

/*      Copyright (C) 1996 : CCLRC, Rutherford Appleton Laboratory */
/*                           Chilton, Didcot, Oxfordshire OX11 0QX */

/*      Release 1.1  29 Oct 1979 (IMS) */
/*      Commented    19 Feb 1980 (KR) */
/*      Release 4.0   2 Oct 1996 (CG) */

/* ARGUMENTS in */
/*      A       array containing the elements of the symmetric */
/*              matrix */
/*      IA      first dimension of A (.GE. N) */
/*      JA      second dimension of A (.GE. N) */
/*      IT      first dimension of array T (.GE. N) */
/*      JT      second dimension of T (.GE. N) */
/*      IDIAG   dimension of vector DIAG (.GE. N) */
/*      ISUB    dimension of vector SUB (.GE. N) */
/*      N       order of matrix A */
/*      TOL     value of RMIN/EPS, where RMIN is the smallest */
/*              positive number exactly representable on the */
/*              computer, and EPS is the smallest positive */
/*              number such that 1.+EPS>1. */
/*      ITEST   error checking option */

/* ARGUMENTS out */
/*      T       contains the orthogonal matrix 'q', the product */
/*              of the householder transformation matrices */
/*      DIAG    contains diagonal elements of tridiagonal matrix */
/*      SUB     contains the N-1 off-diagonal elements of the */
/*              tridiagonal matrix stored in SUB(2) to SUB(N). */
/*              SUB(1)=0. */

/* ROUTINES called */
/*      ERRMES */

/*      SUBROUTINE HOUSE(A,IA,JA,T,IT,JT,DIAG,IDIAG,SUB,ISUB,N,TOL,ITEST) 
*/
/* ***********************************************************************
 */



    /* Parameter adjustments */
    --sub;
    --diag;
    t_dim1 = *it;
    t_offset = t_dim1 + 1;
    t -= t_offset;
    a_dim1 = *ia;
    a_offset = a_dim1 + 1;
    a -= a_offset;

    /* Function Body */

/*     Parameter checking */

    if (*itest != -1) {
	ierror = 0;
	if (*idiag < *n || *isub < *n) {
	    ierror = 2;
	}
	if (*ia < *n || *ja < *n) {
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
	i__2 = i;
	for (j = 1; j <= i__2; ++j) {
	    t[i + j * t_dim1] = a[i + j * a_dim1];
/* L1000: */
	}
/* L1010: */
    }

    if (*n != 1) {
	i__1 = *n;
	for (ii = 2; ii <= i__1; ++ii) {
	    i = *n - ii + 2;
	    l = i - 2;
	    f = t[i + (i - 1) * t_dim1];
	    g = 0.;
	    if (l != 0) {
		i__2 = l;
		for (k = 1; k <= i__2; ++k) {
		    g += t[i + k * t_dim1] * t[i + k * t_dim1];
/* L1020: */
		}

/*     If G is too small for orthogonality to be guaranteed th
e */
/*     transformation is skipped */

	    }
	    h = g + f * f;
	    if (g > *tol) {
		++l;
		g = sqrt(h);
		if (f >= 0.) {
		    g = -g;
		}
		sub[i] = g;
		h -= f * g;
		t[i + (i - 1) * t_dim1] = f - g;
		f = 0.;
		i__2 = l;
		for (j = 1; j <= i__2; ++j) {

/*     Form element of A*U */

		    t[j + i * t_dim1] = t[i + j * t_dim1] / h;
		    g = 0.;
		    i__3 = j;
		    for (k = 1; k <= i__3; ++k) {
			g += t[j + k * t_dim1] * t[i + k * t_dim1];
/* L1030: */
		    }
		    j1 = j + 1;
		    if (j1 <= l) {
			i__3 = l;
			for (k = j1; k <= i__3; ++k) {

/*     Form element of p */

			    g += t[k + j * t_dim1] * t[i + k * t_dim1];
/* L1040: */
			}
		    }
		    sub[j] = g / h;

/*     Form K */

		    f += g * t[j + i * t_dim1];
/* L1050: */
		}

/*     Form reduced A */

		hh = f / (h + h);
		i__2 = l;
		for (j = 1; j <= i__2; ++j) {
		    f = t[i + j * t_dim1];
		    g = sub[j] - hh * f;
		    sub[j] = g;
		    i__3 = j;
		    for (k = 1; k <= i__3; ++k) {
			t[j + k * t_dim1] = t[j + k * t_dim1] - f * sub[k] - 
				g * t[i + k * t_dim1];
/* L1060: */
		    }
/* L1070: */
		}
	    } else {
		sub[i] = f;
		h = 0.;
	    }
	    diag[i] = h;
/* L1080: */
	}
    }

/*     Accumulation of transformation matrices */

    sub[1] = 0.;
    diag[1] = 0.;
    i__1 = *n;
    for (i = 1; i <= i__1; ++i) {
	l = i - 1;
	if ((d__1 = diag[i], abs(d__1)) > veps_(&x)) {
	    i__2 = l;
	    for (j = 1; j <= i__2; ++j) {
		g = 0.;
		i__3 = l;
		for (k = 1; k <= i__3; ++k) {
		    g += t[i + k * t_dim1] * t[k + j * t_dim1];
/* L1090: */
		}
		i__3 = l;
		for (k = 1; k <= i__3; ++k) {
		    t[k + j * t_dim1] -= g * t[k + i * t_dim1];
/* L1100: */
		}
/* L1110: */
	    }
	}
	diag[i] = t[i + i * t_dim1];
	t[i + i * t_dim1] = 1.;
	if (l != 0) {
	    i__2 = l;
	    for (j = 1; j <= i__2; ++j) {
		t[i + j * t_dim1] = 0.;
		t[j + i * t_dim1] = 0.;
/* L1120: */
	    }
	}
/* L1130: */
    }

} /* house_ */

