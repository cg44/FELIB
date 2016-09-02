/* prtmat.f -- translated by f2c (version of 28 August 1991  0:07:02).
   You must link the resulting object file with the libraries:
	-lF77 -lI77 -lm -lc   (in that order)
*/

#include "f2c.h"

/* Table of constant values */

static integer c__1 = 1;


/* Subroutine */ int prtmat_(a, ia, ja, m, n, nout, itest)
doublereal *a;
integer *ia, *ja, *m, *n, *nout, *itest;
{
    /* Initialized data */

    static char srname[6+1] = "PRTMAT";

    /* Format strings */
    static char fmt_9990[] = "(\002 \002,6d12.4)";

    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2;

    /* Builtin functions */
    integer s_wsfe(), do_fio(), e_wsfe();

    /* Local variables */
    static integer i, j;
    extern integer errmes_();
    static integer ierror;

    /* Fortran I/O blocks */
    static cilist io___4 = { 0, 0, 0, fmt_9990, 0 };


/* -----------------------------------------------------------------------
 */
/* PURPOSE */
/*      PRTMAT prints a two-dimensional array in a standard format */

/* HISTORY */

/*      Copyright (C) 1996 : CCLRC, Rutherford Appleton Laboratory */
/*                           Chilton, Didcot, Oxfordshire OX11 0QX */

/*      Release 1.1  29 Oct 1979 (CG) */
/*      Commented    14 Oct 1980 (KR) */
/*      Release 4.0   2 Oct 1996 (CG) */

/* ARGUMENTS in */
/*      A       array of diimension (IA, JA) containing numbers */
/*              to be printed */
/*      IA      first dimension of A (.GE. M) */
/*      JA      second dimension of A (.GE. N) */
/*      M       number of rows of A to be printed */
/*      N       number of columns of A to be printed */
/*      NOUT    fortran unit number */
/*      ITEST   error checking option */

/* ROUTINES called */
/*      ERRMES */

/*      SUBROUTINE PRTMAT(A,IA,JA,M,N,NOUT,ITEST) */
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
	if (*ia < *m || *ja < *n) {
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

    i__1 = *m;
    for (i = 1; i <= i__1; ++i) {
	io___4.ciunit = *nout;
	s_wsfe(&io___4);
	i__2 = *n;
	for (j = 1; j <= i__2; ++j) {
	    do_fio(&c__1, (char *)&a[i + j * a_dim1], (ftnlen)sizeof(
		    doublereal));
	}
	e_wsfe();
/* L1000: */
    }

} /* prtmat_ */

