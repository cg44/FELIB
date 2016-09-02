/* matmul.f -- translated by f2c (version of 28 August 1991  0:07:02).
   You must link the resulting object file with the libraries:
	-lF77 -lI77 -lm -lc   (in that order)
*/

#include "f2c.h"


/* Subroutine */ int matmul_(a, ia, ja, b, ib, jb, c, ic, jc, l, m, n, itest)
doublereal *a;
integer *ia, *ja;
doublereal *b;
integer *ib, *jb;
doublereal *c;
integer *ic, *jc, *l, *m, *n, *itest;
{
    /* Initialized data */

    static char srname[6+1] = "MATMUL";

    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, c_dim1, c_offset, i__1, i__2, 
	    i__3;

    /* Local variables */
    static integer i, j, k;
    static doublereal x;
    extern integer errmes_();
    static integer ierror;

/* -----------------------------------------------------------------------
 */
/* PURPOSE */
/*      MATMUL pre-multiplies matrix B by A, storing the result in C */

/* HISTORY */

/*      Copyright (C) 1996 : CCLRC, Rutherford Appleton Laboratory */
/*                           Chilton, Didcot, Oxfordshire OX11 0QX */

/*      Release 1.1  29 Oct 1979 (IMS) */
/*      Commented    12 Oct 1980 (KR) */
/*      Release 4.0   2 Oct 1996 (CG) */

/* ARGUMENTS in */
/*      A       array of dimension (IA, JA) */
/*      IA      first dimension of A (.GE. L) */
/*      JA      second dimension of A (.GE. M) */
/*      B       array of dimension (IB, JB) */
/*      IB      first dimension of B (.GE. M) */
/*      JB      second dimension of B (.GE. N) */
/*      IC      first dimension of array C (.GE. L) */
/*      JC      second dimension of array C (.GE. N) */
/*      L       number of rows of A to be used in multiplication */
/*      M       number of columns of A and number of rows of B */
/*              to be used in multiplication */
/*      N       number of columns of B to be used in */
/*              multiplication */
/*      ITEST   error checking option */

/* ARGUMENTS out */
/*      C       contains result of matrix multiplication (C=A*B) */

/* ROUTINES called */
/*      ERRMES */


/*      SUBROUTINE MATMUL(A,IA,JA,B,IB,JB,C,IC,JC,L,M,N,ITEST) */
/* ***********************************************************************
 */

    /* Parameter adjustments */
    c_dim1 = *ic;
    c_offset = c_dim1 + 1;
    c -= c_offset;
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
	if (*l > *ic || *n > *jc) {
	    ierror = 4;
	}
	if (*m > *ib || *n > *jb) {
	    ierror = 3;
	}
	if (*l > *ia || *m > *ja) {
	    ierror = 2;
	}
	if (*l <= 0 || *m <= 0 || *n <= 0) {
	    ierror = 1;
	}
	*itest = errmes_(itest, &ierror, srname, 6L);
	if (*itest != 0) {
	    return 0;
	}
    }

/*     Main body */

    i__1 = *l;
    for (i = 1; i <= i__1; ++i) {
	i__2 = *n;
	for (j = 1; j <= i__2; ++j) {
	    x = 0.;
	    i__3 = *m;
	    for (k = 1; k <= i__3; ++k) {
		x += a[i + k * a_dim1] * b[k + j * b_dim1];
/* L1000: */
	    }
	    c[i + j * c_dim1] = x;
/* L1010: */
	}
/* L1020: */
    }

} /* matmul_ */

