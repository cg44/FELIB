/* jaco.f -- translated by f2c (version of 28 August 1991  0:07:02).
   You must link the resulting object file with the libraries:
	-lF77 -lI77 -lm -lc   (in that order)
*/

#include "f2c.h"


/* Subroutine */ int jaco_(a, ia, ja, diag, idiag, sub, isub, n, hband, itest)

doublereal *a;
integer *ia, *ja;
doublereal *diag;
integer *idiag;
doublereal *sub;
integer *isub, *n, *hband, *itest;
{
    /* Initialized data */

    static char srname[6+1] = "JACO  ";

    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3, i__4, i__5;
    doublereal d__1;

    /* Builtin functions */
    double sqrt();

    /* Local variables */
    static integer iugl, maxl, maxr;
    extern doublereal veps_();
    static doublereal b, c, g;
    static integer i, j, k, l, m;
    static doublereal s, u, x, c2;
    static integer j2, n2;
    static doublereal s2, u1;
    static integer jl, jm;
    static doublereal cs;
    static integer ir, kr;
    extern integer errmes_();
    static integer ierror, irr;

/* -----------------------------------------------------------------------
 */
/* PURPOSE */
/*      JACO reduces a real symmetric band matrix to tridiagonal form */
/*      using jacobi rotations, for use with QLVAL or QLVEC. */
/*      the lower triangle of the matrix is stored in a */
/*      rectangular array. */

/* HISTORY */

/*      Copyright (C) 1996 : CCLRC, Rutherford Appleton Laboratory */
/*                           Chilton, Didcot, Oxfordshire OX11 0QX */

/*      Release 1.1  29 Oct 1979 (IMS) */
/*      Commented    26 Feb 1980 (KR) */
/*      Release 4.0   2 Oct 1996 (CG) */

/* ARGUMENTS in */
/*      A       contains the elements of the lower triangle of */
/*              the positive definite band matrix */
/*      IA      first dimension of A (.GE. N) */
/*      JA      second dimension of A (.GE. HBAND) */
/*      IDIAG   dimension of vector DIAG (.GE. N) */
/*      ISUB    dimension of vector SUB (.GE. N) */
/*      N       order of matrix A */
/*      HBAND   semi-bandwidth of matrix A (includes diagonal) */
/*      ITEST   error checking option */

/* ARGUMENTS out */
/*      A       destroyed on succesful exit */
/*      DIAG    contains the diagonal elements of the */
/*              tridiagonal matrix */
/*      SUB     contains the N-1 off-diagonal elements of the */
/*              tridiagonal matrix stored in SUB(2) to SUB(N), */
/*              with SUB(1)=0 */

/* ROUTINES called */
/*      ERRMES */

/*      SUBROUTINE JACO(A,IA,JA,DIAG,IDIAG,SUB,ISUB,N,HBAND,ITEST) */
/* ***********************************************************************
 */



    /* Parameter adjustments */
    --sub;
    --diag;
    a_dim1 = *ia;
    a_offset = a_dim1 + 1;
    a -= a_offset;

    /* Function Body */

/*     Parameter checking */

    if (*itest != -1) {
	ierror = 0;
	if (*isub < *n) {
	    ierror = 4;
	}
	if (*idiag < *n) {
	    ierror = 3;
	}
	if (*ia < *n || *ja < *hband) {
	    ierror = 2;
	}
	if (*n <= 0) {
	    ierror = 1;
	}
	*itest = errmes_(itest, &ierror, srname, 6L);
	if (*itest != 0) {
	    return 0;
	}
    }

/*     Main body */

    g = 0.;
    i__1 = *n;
    for (i = 1; i <= i__1; ++i) {
	i__2 = *hband;
	for (j = 1; j <= i__2; ++j) {
	    l = *hband + 1 - j;
	    k = i + 1 - l;
	    if (k > 0) {
		a[k + l * a_dim1] = a[i + j * a_dim1];
	    }
/* L1000: */
	}
/* L1010: */
    }

    kr = *hband - 1;
    i__1 = kr;
    for (i = 1; i <= i__1; ++i) {
	m = *hband - i;
	i__2 = m;
	for (j = 1; j <= i__2; ++j) {
	    k = *n + 1 - i;
	    l = *hband + 1 - j;
	    a[k + l * a_dim1] = 0.;
/* L1020: */
	}
/* L1030: */
    }

    m = *hband - 1;
    n2 = *n - 2;
    if (n2 >= 1) {
	i__1 = n2;
	for (k = 1; k <= i__1; ++k) {
	    maxr = m;
	    if (*n - k < m) {
		maxr = *n - k;
	    }
	    i__2 = maxr;
	    for (irr = 2; irr <= i__2; ++irr) {
		ir = maxr + 2 - irr;
		kr = k + ir;
		i__3 = *n;
		i__4 = m;
		for (j = kr; i__4 < 0 ? j >= i__3 : j <= i__3; j += i__4) {
		    if (j == kr) {
			if ((d__1 = a[k + (ir + 1) * a_dim1], abs(d__1)) < 
				veps_(&x)) {
			    goto L1070;
			} else {
			    b = -a[k + ir * a_dim1] / a[k + (ir + 1) * a_dim1]
				    ;
			    iugl = k;
			}
		    } else if (abs(g) < veps_(&x)) {
			goto L1070;
		    } else {
			jm = j - m;
			b = -a[jm - 1 + (m + 1) * a_dim1] / g;
			iugl = j - m;
		    }
		    s = 1. / sqrt(b * b + 1.);
		    c = b * s;
		    c2 = c * c;
		    s2 = s * s;
		    cs = c * s;
		    u = c2 * a[j - 1 + a_dim1] - cs * 2. * a[j - 1 + (a_dim1 
			    << 1)] + s2 * a[j + a_dim1];
		    u1 = s2 * a[j - 1 + a_dim1] + cs * 2. * a[j - 1 + (a_dim1 
			    << 1)] + c2 * a[j + a_dim1];
		    a[j - 1 + (a_dim1 << 1)] = cs * (a[j - 1 + a_dim1] - a[j 
			    + a_dim1]) + (c2 - s2) * a[j - 1 + (a_dim1 << 1)];

		    a[j - 1 + a_dim1] = u;
		    a[j + a_dim1] = u1;
		    j2 = j - 2;
		    i__5 = j2;
		    for (l = iugl; l <= i__5; ++l) {
			jl = j - l;
			u = c * a[l + jl * a_dim1] - s * a[l + (jl + 1) * 
				a_dim1];
			a[l + (jl + 1) * a_dim1] = s * a[l + jl * a_dim1] + c 
				* a[l + (jl + 1) * a_dim1];
			a[l + jl * a_dim1] = u;
/* L1040: */
		    }
		    jm = j - m;
		    if (j != kr) {
			a[jm - 1 + (m + 1) * a_dim1] = c * a[jm - 1 + (m + 1) 
				* a_dim1] - s * g;
		    }
		    maxl = m - 1;
		    if (*n - j < m - 1) {
			maxl = *n - j;
		    }
		    if (maxl > 0) {
			i__5 = maxl;
			for (l = 1; l <= i__5; ++l) {
			    u = c * a[j - 1 + (l + 2) * a_dim1] - s * a[j + (
				    l + 1) * a_dim1];
			    a[j + (l + 1) * a_dim1] = s * a[j - 1 + (l + 2) * 
				    a_dim1] + c * a[j + (l + 1) * a_dim1];
			    a[j - 1 + (l + 2) * a_dim1] = u;
/* L1050: */
			}
		    }
		    if (j + m <= *n) {
			g = -s * a[j + (m + 1) * a_dim1];
			a[j + (m + 1) * a_dim1] = c * a[j + (m + 1) * a_dim1];

		    }
/* L1060: */
		}
L1070:
		;
	    }
/* L1080: */
	}
    }

    sub[1] = 0.;
    i__1 = *n;
    for (i = 1; i <= i__1; ++i) {
	diag[i] = a[i + a_dim1];
/* L1090: */
    }

    if (2 <= *n) {
	i__1 = *n;
	for (i = 2; i <= i__1; ++i) {
	    sub[i] = a[i - 1 + (a_dim1 << 1)];
/* L1100: */
	}
    }
/*     $st$ unreachable comments ... */


} /* jaco_ */

