/* matnul.f -- translated by f2c (version of 28 August 1991  0:07:02).
   You must link the resulting object file with the libraries:
	-lF77 -lI77 -lm -lc   (in that order)
*/

#include "f2c.h"


/* Subroutine */ int matnul_(a, ia, ja, m, n, itest)
doublereal *a;
integer *ia, *ja, *m, *n, *itest;
{
    /* Initialized data */

    static char srname[6+1] = "MATNUL";

    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2;

    /* Local variables */
    static integer i, j;
    extern integer errmes_();
    static integer ierror;

/* -----------------------------------------------------------------------
 */
/* PURPOSE */
/*      MATNUL sets matrix A to the null matrix */

/* HISTORY */

/*      Copyright (C) 1996 : CCLRC, Rutherford Appleton Laboratory */
/*                           Chilton, Didcot, Oxfordshire OX11 0QX */

/*      Release 1.1  29 Oct 1979 (CG) */
/*      Commented    12 Oct 1980 (KR) */
/*      Release 4.0   2 Oct 1996 (CG) */

/* ARGUMENTS in */
/*      IA      first dimension of array A (.GE. M) */
/*      JA      second dimension of array A (.GE. N) */
/*      M       number of rows of A to be set to zero */
/*      N       number of columns of A to be set to zero */
/*      ITEST   error checking option */

/* ARGUMENTS out */
/*      A       array of dimension (IA, JA). A(I, J)=0 for */
/*              I=1(1)M and J=1(1)N */

/* ROUTINES called */
/*      ERRMES */


/*      SUBROUTINE MATNUL(A,IA,JA,M,N,ITEST) */
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
	if (*m <= 0 || *n <= 0) {
	    ierror = 1;
	}
	*itest = errmes_(itest, &ierror, srname, 6L);
	if (*itest != 0) {
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

} /* matnul_ */

