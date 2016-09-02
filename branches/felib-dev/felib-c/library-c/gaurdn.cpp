/* gaurdn.f -- translated by f2c (version of 28 August 1991  0:07:02).
   You must link the resulting object file with the libraries:
	-lF77 -lI77 -lm -lc   (in that order)
*/

#include "f2c.h"


/* Subroutine */ int gaurdn_(a, ia, ja, al, ial, jal, n, hband, ropiv, iropiv,
	 itest)
doublereal *a;
integer *ia, *ja;
doublereal *al;
integer *ial, *jal, *n, *hband, *ropiv, *iropiv, *itest;
{
    /* Initialized data */

    static doublereal one = 1.;
    static char srname[6+1] = "GAURDN";
    static doublereal zero = 0.;

    /* System generated locals */
    integer a_dim1, a_offset, al_dim1, al_offset, i__1, i__2, i__3;
    doublereal d__1;

    /* Local variables */
    extern doublereal veps_();
    static integer i, j, k, l, m;
    static doublereal x, y;
    static integer jtest, ik, jj, ir, jr, iw;
    extern integer errmes_();
    static integer ierror;
    static doublereal eps;

/* -----------------------------------------------------------------------
 */
/* PURPOSE */
/*      GAURDN decomposes A real unsymmetric matrix of order N into */
/*      triangular matrices using gaussian elimination with */
/*      partial pivoting */

/* HISTORY */

/*      Copyright (C) 1996 : CCLRC, Rutherford Appleton Laboratory */
/*                           Chilton, Didcot, Oxfordshire OX11 0QX */

/*      Release 1.1  29 Oct 1979 (CG) */
/*      Commented    10 Oct 1980 (KR) */
/*      Release 4.0   2 Oct 1996 (CG) */

/* ARGUMENTS in */
/*      A       array of dimensions (IA, JA).  on entry, contains */
/*              the elements of the band matrix */
/*      IA      first dimension of A (.GE. N) */
/*      JA      second dimension of A (.GE. MIN(2*HBAND-1,N)) */
/*      IAL     first dimension of array AL (.GE. N) */
/*      JAL     second dimension of array AL (.GE. HBAND-1) */
/*      N       order of matrix A */
/*      HBAND   semi-bandwidth of matrix A */
/*      IROPIV  dimension of vector ROPIV */
/*      ITEST   error checking option */

/* ARGUMENTS out */
/*      A       array of dimensions (IA, JA).  on exit, contains */
/*              the upper triangular matrix U, with the */
/*              diagonal elements of U stored as their */
/*              reciprocals.  the i'th row of U is stored in the */
/*              i'th row of A, starting with the diagonal */
/*              element of U in A(I, 1) */
/*      ROPIV   vector of length IROPIV, containing details of */
/*              row interchanges.  if no interchange occurs at */
/*              the r'th major step then ROPIV(R)=R; if the R */
/*              and J rows are interchanged then ROPIV(R)=J */

/* ROUTINES called */
/*      VEPS    ERRMES */

/*      SUBROUTINE GAURDN(A,IA,JA,AL,IAL,JAL,N,HBAND,ROPIV,IROPIV,ITEST) 
*/
/* ***********************************************************************
 */


/*      INTRINSIC MIN */
/*      EXTERNAL ERRMES, VEPS */

    /* Parameter adjustments */
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
	    ierror = 4;
	}
	if (*ial < *n || *jal < *hband) {
	    ierror = 3;
	}
	if (*iropiv < *n) {
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

    ierror = 5;
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

/*     Zeros inserted */

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

/*     Row norms of A calculated and their reciprocals stored in 
*/
/*     first column of AL */

	    al[i + al_dim1] = one / x;
	} else {
	    goto L1120;
	}
/* L1060: */
    }

    ierror = 6;
    i__1 = *n;
    for (ir = 1; ir <= i__1; ++ir) {
	x = zero;
/* Computing MIN */
	i__2 = ir + *hband - 1;
	m = min(i__2,*n);
	i__2 = m;
	for (i = ir; i <= i__2; ++i) {
	    y = (d__1 = a[i + a_dim1], abs(d__1)) * al[i + al_dim1];
	    if (y > x) {
		x = y;
		j = i;
	    }
/* L1070: */
	}

/*     IR'TH pivot element selected. */

	if (x < eps) {
	    goto L1130;
	} else {
	    ropiv[ir] = j;
	    if (j != ir) {
		i__2 = iw;
		for (i = 1; i <= i__2; ++i) {
		    x = a[ir + i * a_dim1];
		    a[ir + i * a_dim1] = a[j + i * a_dim1];
		    a[j + i * a_dim1] = x;
/* L1080: */
		}

/*     Row IR and J interchanged. */


		al[j + al_dim1] = al[ir + al_dim1];
	    }
	    jr = ir + 1;
	    y = one / a[ir + a_dim1];
	    if (jr <= m) {
		i__2 = m;
		for (i = jr; i <= i__2; ++i) {
		    x = a[i + a_dim1] * y;
		    if (iw >= 2) {
			i__3 = iw;
			for (j = 2; j <= i__3; ++j) {
			    a[i + (j - 1) * a_dim1] = a[i + j * a_dim1] - x * 
				    a[ir + j * a_dim1];
/* L1090: */
			}
		    }
		    ik = i - ir;
		    al[ir + ik * al_dim1] = x;
		    a[i + iw * a_dim1] = zero;
/* L1100: */
		}
	    }

/*     Elimination completed */

	    a[ir + a_dim1] = y;
	}
/* L1110: */
    }

    return 0;
L1120:
    ir = i;

L1130:
    a[ir + a_dim1] = zero;
    *itest = errmes_(&jtest, &ierror, srname, 6L);

} /* gaurdn_ */

