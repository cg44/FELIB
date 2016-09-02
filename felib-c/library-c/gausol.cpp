/* gausol.f -- translated by f2c (version of 28 August 1991  0:07:02).
   You must link the resulting object file with the libraries:
	-lF77 -lI77 -lm -lc   (in that order)
*/

#include "f2c.h"


/* Subroutine */ int gausol_(a, ia, ja, al, ial, jal, n, hband, ropiv, iropiv,
	 r, ir, itest)
doublereal *a;
integer *ia, *ja;
doublereal *al;
integer *ial, *jal, *n, *hband, *ropiv, *iropiv;
doublereal *r;
integer *ir, *itest;
{
    /* Initialized data */

    static doublereal one = 1.;
    static char srname[6+1] = "GAUSOL";
    static doublereal zero = 0.;

    /* System generated locals */
    integer a_dim1, a_offset, al_dim1, al_offset, i__1, i__2, i__3;
    doublereal d__1;

    /* Local variables */
    extern doublereal veps_();
    static integer i, j, k, l, m;
    static doublereal x, y;
    static integer jtest, ii, ik, jj, kk, jr, iw;
    extern integer errmes_();
    static integer ierror;
    static doublereal eps;
    static integer iro;

/* -----------------------------------------------------------------------
 */
/* PURPOSE */
/*      GAUSOLc alculates the solution of a set of unsymmetric real */
/*      banded equations with A single rhs using gaussian */
/*      elimination with partial pivoting */

/* HISTORY */

/*      Copyright (C) 1996 : CCLRC, Rutherford Appleton Laboratory */
/*                           Chilton, Didcot, Oxfordshire OX11 0QX */

/*      Release 1.1  29 Oct 1979 (CG) */
/*      Commented    11 Oct 1980 (KR) */
/*      Release 4.0   2 Oct 1996 (CG) */

/* ARGUMENTS in */
/*      A       array of dimension (IA, JA).  On entry, contains */
/*              the elements of the band matrix */
/*      IA      first dimension of A (.GE. N) */
/*      JA      second dimension of A (.GE. MIN(2*HBAND-1, N)) */
/*      IAL     first dimension of array AL (.GE. N) */
/*      JAL     second dimension of array AL (.GE. HBAND-1) */
/*      N       order of band matrix A */
/*      HBAND   semi-bandwidth of matrix A */
/*      IROPIV  dimension of vector ROPIV (.GE. N) */
/*      R       on entry, contains the vector of rhs's */
/*      IR      dimension of R (.GE. N) */
/*      ITEST   error checking option */

/* ARGUMENTS out */
/*      A       on exit, contains the upper triangular matrix U, */
/*              with the diagonal elements of U stored as their */
/*              reciprocals.  the i'th row of U is stored in the */
/*              i'th row of A, with the diagonal element of U in */
/*              the location A(I, 1) */
/*      AL      contains the sub-diagonal elements of L, the */
/*              lower triangular matrix.  The multipliers L(I, R) */
/*              obtained at the r'th major step of the */
/*              elimination are stored in A(R, I-R) */
/*      ROPIV   contains details of the row interchanges. */
/*              ROPIV(R)=R if no interchange occurs at the r'th */
/*              major step; if rows R and J are interchanged */
/*              then ROPIV(R)=J */
/*      R       on exit, contains the solution vector */

/* ROUTINES called */
/*      VEPS    ERRMES */


/*      SUBROUTINE GAUSOL(A,IA,JA,AL,IAL,JAL,N,HBAND,ROPIV,IROPIV,R,IR, */

/*     *                  ITEST) */
/* ***********************************************************************
 */



    /* Parameter adjustments */
    --r;
    --ropiv;
    al_dim1 = *ial;
    al_offset = al_dim1 + 1;
    al -= al_offset;
    a_dim1 = *ia;
    a_offset = a_dim1 + 1;
    a -= a_offset;

    /* Function Body */

/*     Parameter checking */

    jtest = *itest;
    if (jtest != -1) {
	ierror = 0;
	if (*ia < *n || *ja < (*hband << 1) - 1) {
	    ierror = 5;
	}
	if (*ial < *n || *jal < *hband) {
	    ierror = 4;
	}
	if (*iropiv < *n) {
	    ierror = 3;
	}
	if (*ir < *n) {
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

    ierror = 6;
    eps = veps_(&x);
/* Computing MIN */
    i__1 = *n, i__2 = (*hband << 1) - 1;
    iw = min(i__1,i__2);
    m = *hband;
    k = iw - *hband;
    if (k > 0) {
	i__1 = k;
	for (i = 1; i <= i__1; ++i) {
	    l = iw - m;
	    i__2 = m;
	    for (j = 1; j <= i__2; ++j) {
		jj = j + l;
		a[i + j * a_dim1] = a[i + jj * a_dim1];
/* L1000: */
	    }
	    ++m;
	    i__2 = iw;
	    for (j = m; j <= i__2; ++j) {
		a[i + j * a_dim1] = zero;
/* L1010: */
	    }
/* L1020: */
	}
    }
    m = *n - iw + *hband + 1;
    j = iw + 1;
    if (m <= *n) {
	i__1 = *n;
	for (i = m; i <= i__1; ++i) {
	    --j;
	    i__2 = iw;
	    for (k = j; k <= i__2; ++k) {
		a[i + k * a_dim1] = zero;
/* L1030: */
	    }
/* L1040: */
	}

/*     Insert zeros */

    }

    i__1 = *n;
    for (i = 1; i <= i__1; ++i) {
	x = zero;
	i__2 = iw;
	for (j = 1; j <= i__2; ++j) {
	    x += (d__1 = a[i + j * a_dim1], abs(d__1));
/* L1050: */
	}
	if (x > zero) {

/*     ROPIV norms of A calculated and theiro reciprocals stored i
n */
/*     firost column of AL */

	    al[i + al_dim1] = one / x;
	} else {
	    goto L1170;
	}
/* L1060: */
    }

    ierror = 7;
    i__1 = *n;
    for (iro = 1; iro <= i__1; ++iro) {
	x = zero;
/* Computing MIN */
	i__2 = iro + *hband - 1;
	m = min(i__2,*n);
	i__2 = m;
	for (i = iro; i <= i__2; ++i) {
	    y = (d__1 = a[i + a_dim1], abs(d__1)) * al[i + al_dim1];
	    if (y > x) {
		x = y;
		j = i;
	    }
/* L1070: */
	}

/*     IRO'TH pivot element selected */

	if (x < eps) {
	    goto L1180;
	} else {
	    ropiv[iro] = j;
	    if (j != iro) {
		i__2 = iw;
		for (i = 1; i <= i__2; ++i) {
		    x = a[iro + i * a_dim1];
		    a[iro + i * a_dim1] = a[j + i * a_dim1];
		    a[j + i * a_dim1] = x;
/* L1080: */
		}

/*     Row pivots IRO and J interchanged */

		al[j + al_dim1] = al[iro + al_dim1];
	    }
	    jr = iro + 1;
	    y = one / a[iro + a_dim1];
	    if (jr <= m) {
		i__2 = m;
		for (i = jr; i <= i__2; ++i) {
		    x = a[i + a_dim1] * y;
		    if (iw >= 2) {
			i__3 = iw;
			for (j = 2; j <= i__3; ++j) {
			    a[i + (j - 1) * a_dim1] = a[i + j * a_dim1] - x * 
				    a[iro + j * a_dim1];
/* L1090: */
			}
		    }
		    ik = i - iro;
		    al[iro + ik * al_dim1] = x;
		    a[i + iw * a_dim1] = zero;
/* L1100: */
		}
	    }

/*     Elimination complete */

	    a[iro + a_dim1] = y;
	}
/* L1110: */
    }

    m = *hband - 1;
    i__1 = *n;
    for (k = 1; k <= i__1; ++k) {
/* Computing MIN */
	i__2 = m + 1;
	m = min(i__2,*n);
	j = ropiv[k];
	if (j != k) {
	    x = r[k];

/*     Row pivots K and J interchanged */

	    r[k] = r[j];
	    r[j] = x;
	}
	ik = k + 1;
	if (ik > m) {
	    goto L1140;
	} else {
	    x = r[k];
	    i__2 = m;
	    for (i = ik; i <= i__2; ++i) {
		ii = i - k;
		r[i] -= x * al[k + ii * al_dim1];
/* L1120: */
	    }

/*     Forward substitution complete */

	}
/* L1130: */
    }

L1140:
    i__1 = *n;
    for (k = 1; k <= i__1; ++k) {
	m = min(k,iw);
	i = *n + 1 - k;
	ii = i - 1;
	y = a[i + a_dim1];
	x = r[i];
	if (m != 1) {
	    i__2 = m;
	    for (j = 2; j <= i__2; ++j) {
		kk = j + ii;
		x -= a[i + j * a_dim1] * r[kk];
/* L1150: */
	    }
	}

/*     Backward substitution complete */

	r[i] = x * y;
/* L1160: */
    }

    return 0;
L1170:
    iro = i;

L1180:
    a[iro + a_dim1] = zero;
    *itest = errmes_(&jtest, &ierror, srname, 6L);

} /* gausol_ */

