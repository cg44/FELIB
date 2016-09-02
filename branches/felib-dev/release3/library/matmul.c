/* matmul.f -- translated by f2c (version 20000121).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"

/* *********************************************************************** */
/* Subroutine */ int matmul_(doublereal *a, integer *ia, integer *ja, 
	doublereal *b, integer *ib, integer *jb, doublereal *c__, integer *ic,
	 integer *jc, integer *l, integer *m, integer *n, integer *itest)
{
    /* Initialized data */

    static struct {
	char e_1[8];
	doublereal e_2;
	} equiv_7 = { " MATMUL ", 0. };

#define srname (*(doublereal *)&equiv_7)


    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, c_dim1, c_offset, i__1, i__2, 
	    i__3;

    /* Local variables */
    static integer i__, j, k;
    static doublereal x;
    extern integer errmes_(integer *, integer *, doublereal *);
    static integer ierror;

/* ----------------------------------------------------------------------- */
/* PURPOSE */
/*      PRE-MULTIPLIES MATRIX B BY A, STORING THE RESULT IN C */

/* HISTORY */

/*      COPYRIGHT (C) 1997 : CLRC, RUTHERFORD APPLETON LABORATORY */
/*                           CHILTON, DIDCOT, OXFORDSHIRE OX11 0QX */

/*      RELEASE 1.1  29 OCT 1979 (IMS) */
/*      COMMENTED    12 OCT 1980 (KR) */

/* ARGUMENTS IN */
/*      A       ARRAY OF DIMENSION (IA,JA) */
/*      IA      FIRST DIMENSION OF A (.GE.L) */
/*      JA      SECOND DIMENSION OF A (.GE.M) */
/*      B       ARRAY OF DIMENSION (IB,JB) */
/*      IB      FIRST DIMENSION OF B (.GE.M) */
/*      JB      SECOND DIMENSION OF B (.GE.N) */
/*      IC      FIRST DIMENSION OF ARRAY C (.GE.L) */
/*      JC      SECOND DIMENSION OF ARRAY C (.GE.N) */
/*      L       NUMBER OF ROWS OF A TO BE USED IN MULTIPLICATION */
/*      M       NUMBER OF COLUMNS OF A AND NUMBER OF ROWS OF B */
/*              TO BE USED IN MULTIPLICATION */
/*      N       NUMBER OF COLUMNS OF B TO BE USED IN */
/*              MULTIPLICATION */
/*      ITEST   ERROR CHECKING OPTION */

/* ARGUMENTS OUT */
/*      C       CONTAINS RESULT OF MATRIX MULTIPLICATION (C=A*B) */

/* ROUTINES CALLED */
/*      ERRMES */


/*     SUBROUTINE MATMUL(A, IA, JA, B, IB, JB, C, IC, JC, L, M, N, */
/*    *     ITEST) */
/* *********************************************************************** */

    /* Parameter adjustments */
    a_dim1 = *ia;
    a_offset = 1 + a_dim1 * 1;
    a -= a_offset;
    b_dim1 = *ib;
    b_offset = 1 + b_dim1 * 1;
    b -= b_offset;
    c_dim1 = *ic;
    c_offset = 1 + c_dim1 * 1;
    c__ -= c_offset;

    /* Function Body */

/*     PARAMETER CHECKING */

    if (*itest == -1) {
	goto L1010;
    }
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
    *itest = errmes_(itest, &ierror, &srname);
    if (*itest != 0) {
	return 0;
    }

/*     MAIN BODY */

L1010:
    i__1 = *l;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = *n;
	for (j = 1; j <= i__2; ++j) {
	    x = 0.;
	    i__3 = *m;
	    for (k = 1; k <= i__3; ++k) {
		x += a[i__ + k * a_dim1] * b[k + j * b_dim1];
/* L1020: */
	    }
	    c__[i__ + j * c_dim1] = x;
/* L1030: */
	}
/* L1040: */
    }

    return 0;
} /* matmul_ */

#undef srname


