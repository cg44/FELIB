/* matidn.f -- translated by f2c (version of 28 August 1991  0:07:02).
   You must link the resulting object file with the libraries:
	-lF77 -lI77 -lm -lc   (in that order)
*/

#include "f2c.h"


/* Subroutine */ int matidn_(a, ia, ja, m, n, itest)
doublereal *a;
integer *ia, *ja, *m, *n, *itest;
{
    /* Initialized data */

    static char srname[6+1] = "MATIDN";

    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2;

    /* Local variables */
    static integer i, j, l;
    extern integer errmes_();
    static integer ierror;

/* -----------------------------------------------------------------------
 */
/* PURPOSE */
/*      MATIDN sets the matrix A to the identity matrix */

/* HISTORY */

/*      Copyright (C) 1996 : CCLRC, Rutherford Appleton Laboratory */
/*                           Chilton, Didcot, Oxfordshire OX11 0QX */

/*      Release 1.1  29 Oct 1979 (CG) */
/*      Commented    12 Oct 1980 (KR) */
/*      Release 4.0   2 Oct 1996 (CG) */

/* ARGUMENTS in */
/*      IA      first dimension of array A (.GE. M) */
/*      JA      second dimension of array A (.GE. N) */
/*      M       number of rows of A to be assigned values */
/*      N       number of columns of A to be assigned values */
/*      ITEST   error checking option */

/* ARGUMENTS out */
/*      A       array of dimension (IA, JA).  A(I, J) is set to 1 */
/*              if I=J, 0 otherwise */

/* ROUTINES called */
/*      ERRMES */


/*      SUBROUTINE MATIDN(A,IA,JA,M,N,ITEST) */
/* ***********************************************************************
 */

    /* Parameter adjustments */
    a_dim1 = *ia;
    a_offset = a_dim1 + 1;
    a -= a_offset;

    /* Function Body */

/*     Parameter checking */

    if (*itest != -1) {
	ierror = 0;
	if (*m > *ia || *n > *ja) {
	    ierror = 2;
	}
	if (*n <= 0 || *m <= 0) {
	    ierror = 1;
	}
	*itest = errmes_(itest, &ierror, srname, 6L);
	if (ierror != 0) {
	    return 0;
	}
    }

/*     Main body */

    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	i__2 = *m;
	for (i = 1; i <= i__2; ++i) {
	    a[i + j * a_dim1] = 0.;
/* L1000: */
	}
/* L1010: */
    }

/*     Set diagonal entry */

    l = min(*m,*n);
    i__1 = l;
    for (i = 1; i <= i__1; ++i) {
	a[i + i * a_dim1] = 1.;
/* L1020: */
    }

} /* matidn_ */

