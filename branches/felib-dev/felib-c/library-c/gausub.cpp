/* gausub.f -- translated by f2c (version of 28 August 1991  0:07:02).
   You must link the resulting object file with the libraries:
	-lF77 -lI77 -lm -lc   (in that order)
*/

#include "f2c.h"


/* Subroutine */ int gausub_(a, ia, ja, al, ial, jal, n, hband, ropiv, iropiv,
	 r, ir, itest)
doublereal *a;
integer *ia, *ja;
doublereal *al;
integer *ial, *jal, *n, *hband, *ropiv, *iropiv;
doublereal *r;
integer *ir, *itest;
{
    /* Initialized data */

    static char srname[6+1] = "GAUSUB";

    /* System generated locals */
    integer a_dim1, a_offset, al_dim1, al_offset, i__1, i__2;

    /* Local variables */
    static integer i, j, k, m;
    static doublereal x, y;
    static integer ii, ik, kk, iw;
    extern integer errmes_();
    static integer ierror;

/* -----------------------------------------------------------------------
 */
/* PURPOSE */
/*      GAUSUB calculates the solution of a set of unsymmetric real */
/*      banded linear equations with a single rhs.  The banded */
/*      matrix has previously been decomposed into triangular */
/*      matrices using GAURDN */

/* HISTORY */

/*      Copyright (C) 1996 : CCLRC, Rutherford Appleton Laboratory */
/*                           Chilton, Didcot, Oxfordshire OX11 0QX */

/*      Release 1.1  29 Oct 1979 (CG) */
/*      Commented    11 Oct 1980 (KR) */
/*      Release 4.0   2 Oct 1996 (CG) */

/* ARGUMENTS in */
/*      A       array of dimensions (IA, JA).  On entry, contains */
/*              the elements of the band matrix in LU form, */
/*              after processing by GAURDN */
/*      IA      first dimension of A (.GE. N) */
/*      JA      second dimension of A (.GE. MIN(2*HBAND-1, N)) */
/*      IAL     first dimension of AL (.GE. N) */
/*      JAL     second dimension of AL (.GE. HBAND-1) */
/*      N       order of banded matrix A */
/*      HBAND   semi-bandwidth of A */
/*      ROPIV   vector of dimension IROPIV.  Contains details */
/*              of row interchanges performed by GAURDN */
/*      IROPIV  dimension of ROPIV (.GE. N) */
/*      R       on entry, contains the vector of the rhs, */
/*              length IR */
/*      IR      dimension of R (.GE. N) */
/*      ITEST   error checking option */

/* ARGUMENTS out */
/*      R       on exit, contains solution vector */

/* ROUTINES called */
/*      ERRMES */


/*      SUBROUTINE GAUSUB(A,IA,JA,AL,IAL,JAL,N,HBAND,ROPIV,IROPIV,R,IR, */

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

    if (*itest != -1) {
	ierror = 0;
	if (*ir < *n) {
	    ierror = 5;
	}
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
	*itest = errmes_(itest, &ierror, srname, 6L);
	if (*itest != 0) {
	    return 0;
	}
    }

/*     Main body */

/* Computing MIN */
    i__1 = *n, i__2 = (*hband << 1) - 1;
    iw = min(i__1,i__2);
    m = *hband - 1;
    i__1 = *n;
    for (k = 1; k <= i__1; ++k) {
/* Computing MIN */
	i__2 = m + 1;
	m = min(i__2,*n);
	j = ropiv[k];
	if (j != k) {
	    x = r[k];

/*     Rows K and J interchanged */

	    r[k] = r[j];
	    r[j] = x;
	}
	ik = k + 1;
	if (ik > m) {
	    goto L1020;
	} else {
	    x = r[k];
	    i__2 = m;
	    for (i = ik; i <= i__2; ++i) {
		ii = i - k;
		r[i] -= x * al[k + ii * al_dim1];
/* L1000: */
	    }

/*     Forward substitution complete */

	}
/* L1010: */
    }

L1020:
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
/* L1030: */
	    }
	}

/*     Backward substitution complete */

	r[i] = x * y;
/* L1040: */
    }

} /* gausub_ */

