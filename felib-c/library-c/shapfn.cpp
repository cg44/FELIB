/* shapfn.f -- translated by f2c (version of 28 August 1991  0:07:02).
   You must link the resulting object file with the libraries:
	-lF77 -lI77 -lm -lc   (in that order)
*/

#include "f2c.h"


/* Subroutine */ int shapfn_(n, in, jn, fun, ifun, nodel, nde, itest)
doublereal *n;
integer *in, *jn;
doublereal *fun;
integer *ifun, *nodel, *nde, *itest;
{
    /* Initialized data */

    static char srname[6+1] = "SHAPFN";

    /* System generated locals */
    integer n_dim1, n_offset, i__1, i__2;

    /* Local variables */
    static integer i, j, ij;
    extern integer errmes_();
    extern /* Subroutine */ int matnul_();
    static integer ierror;

/* -----------------------------------------------------------------------
 */
/* PURPOSE */
/*     SHAPFN constructs the shape function matrix N for a system of coupl
ed*/
/*      differential equations.  This routine is most frequently used */
/*      in forming the consistent element and system mass matrices */

/* HISTORY */

/*      Copyright (C) 1996 : CCLRC, Rutherford Appleton Laboratory */
/*                           Chilton, Didcot, Oxfordshire OX11 0QX */

/*      Release 1.1  29 Oct 1979 (CG) */
/*      Commented    22 Oct 1980 (KR) */

/* ARGUMENTS IN */
/*      IN      first dimension of array N (.GE. NDE) */
/*      JN      second dimension of N (.GE. NODEL*NDE) */
/*      FUN     vector of length IFUN. Contains values of the */
/*              shape functions at the point where shape */
/*              function matrix required */
/*      IFUN    length of vector FUN (.GE. NODEL) */
/*      NODEL   number of shape functions to be used in */
/*              constructing shape function matrix N (usually */
/*              number of nodes IN element under consideration) */
/*      NDE     number of coupled equations being considered */
/*      ITEST   error checking option */

/* ARGUMENTS out */
/*      N       array of dimension (IN, JN).  contains the */
/*              values of the shape function matrix */

/* ROUTINES called */
/*      ERRMES */

/*     SUBROUTINE SHAPFN(N,IN,JN,FUN,IFUN,NODEL,NDE,ITEST) */
/* ***********************************************************************
 */
/*      Release 4.0   2 Oct 1996 (CG) */

    /* Parameter adjustments */
    --fun;
    n_dim1 = *in;
    n_offset = n_dim1 + 1;
    n -= n_offset;

    /* Function Body */

/*     Parameter checking */

    if (*itest != -1) {
	ierror = 0;
	if (*nodel <= 0 || *nde <= 0) {
	    ierror = 1;
	}
	if (*in < *nde || *jn < *nodel * *nde) {
	    ierror = 2;
	}
	if (*ifun < *nodel) {
	    ierror = 3;
	}
	*itest = errmes_(itest, &ierror, srname, 6L);
	if (*itest != 0) {
	    return 0;
	}
    }

/*     Main body */

    i__1 = *nodel * *nde;
    matnul_(&n[n_offset], in, jn, nde, &i__1, itest);

    i__1 = *nodel;
    for (j = 1; j <= i__1; ++j) {
	i__2 = *nde;
	for (i = 1; i <= i__2; ++i) {
	    ij = (j - 1) * *nde + i;
	    n[i + ij * n_dim1] = fun[j];
/* L1000: */
	}
/* L1010: */
    }

} /* shapfn_ */

