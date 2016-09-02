/* matran.f -- translated by f2c (version of 28 August 1991  0:07:02).
   You must link the resulting object file with the libraries:
	-lF77 -lI77 -lm -lc   (in that order)
*/

#include "f2c.h"


/* Subroutine */ int matran_(a, ia, ja, b, ib, jb, m, n, itest)
doublereal *a;
integer *ia, *ja;
doublereal *b;
integer *ib, *jb, *m, *n, *itest;
{
    /* Initialized data */

    static char srname[6+1] = "MATRAN";

    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, i__1, i__2;

    /* Local variables */
    static integer i, j;
    extern integer errmes_();
    static integer ierror;

/* -----------------------------------------------------------------------
 */
/* PURPOSE */
/*      MATRAN forms the transpose of the matrix A */

/* HISTORY */

/*      Copyright (C) 1996 : CCLRC, Rutherford Appleton Laboratory */
/*                           Chilton, Didcot, Oxfordshire OX11 0QX */

/*      Release 1.1  20 Oct 1979 (IMS) */
/*      Commented    12 Oct 1980 (KR) */
/*      Release 4.0   2 Oct 1996 (CG) */

/* ARGUMENTS in */
/*      A       array of dimension (IA, JA) to be transposed */
/*      IA      first dimension of A (.GE. M) */
/*      JA      second dimension of A (.GE. N) */
/*      IB      first dimension of array B (.GE. N) */
/*      JB      second dimension of array B (.GE. M) */
/*      M       number of rows of A to be transposed */
/*      N       number of columns of A to be transposed */
/*      ITEST   error checking option */

/* ARGUMENTS out */
/*      B       array of dimension (IB, JB). B(J, I) = A(I, J) for */
/*              I=1(1)M and J=1(1)N */

/* ROUTINES called */
/*      ERRMES */


/*      SUBROUTINE MATRAN(A,IA,JA,B,IB,JB,M,N,ITEST) */
/* ***********************************************************************
 */

    /* Parameter adjustments */
    b_dim1 = *ib;
    b_offset = b_dim1 + 1;
    b -= b_offset;
    a_dim1 = *ia;
    a_offset = a_dim1 + 1;
    a -= a_offset;

    /* Function Body */

/*     Parameter checking */

    if (*itest != -1) {
	ierror = 0;
	if (*n > *ib || *m > *jb) {
	    ierror = 3;
	}
	if (*m > *ia || *n > *ja) {
	    ierror = 2;
	}
	if (*m <= 0 || *n <= 0) {
	    ierror = 1;
	}
	*itest = errmes_(itest, &ierror, srname, 6L);
	if (ierror != 0) {
	    return 0;
	}
    }

/*     Main body */

    i__1 = *m;
    for (i = 1; i <= i__1; ++i) {
	i__2 = *n;
	for (j = 1; j <= i__2; ++j) {
	    b[j + i * b_dim1] = a[i + j * a_dim1];
/* L1000: */
	}
/* L1010: */
    }

} /* matran_ */

