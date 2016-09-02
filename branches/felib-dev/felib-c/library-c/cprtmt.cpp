/* cprtmt.f -- translated by f2c (version of 28 August 1991  0:07:02).
   You must link the resulting object file with the libraries:
	-lF77 -lI77 -lm -lc   (in that order)
*/

#include "f2c.h"

/* Table of constant values */

static integer c__1 = 1;


/* Subroutine */ int cprtmt_(a, ia, ja, m, n, nout, itest)
doublereal *a;
integer *ia, *ja, *m, *n, *nout, *itest;
{
    /* Initialized data */

    static char srname[6+1] = "CPRTMT";

    /* Format strings */
    static char fmt_9990[] = "(\002 \002,3(2x,2d12.4))";

    /* System generated locals */
    integer a_dim2, a_offset, i__1, i__2;

    /* Builtin functions */
    integer s_wsfe(), do_fio(), e_wsfe();

    /* Local variables */
    static integer i, j, k;
    extern integer errmes_();
    static integer ierror;

    /* Fortran I/O blocks */
    static cilist io___3 = { 0, 0, 0, fmt_9990, 0 };


/* -----------------------------------------------------------------------
 */
/* PURPOSE */
/*      CPRTMT prints A complex matrix in a standard format */

/* HISTORY */

/*      Copyright (C) 1996 : CCLRC, Rutherford Appleton Laboratory */
/*                           Chilton, Didcot, Oxfordshire OX11 0QX */

/*      Release 3.0  29 Oct 1984 (CRIE) */
/*      Comments      1 Nov 1985 (CG) */
/*      Release 4.0   2 Oct 1996 (CG) */

/* ARGUMENTS in */
/*      A       matrix of dimension (2, IA, JA) containing numbers to be 
*/
/*              printed */
/*      IA      first dimension of A (.GE. M) */
/*      JA      second dimension of A (.GE. N) */
/*      M       number of rows of A to be printed (.LE. IA) */
/*      N       number of columes of A to be printed (.LE. JA) */
/*      NOUT    fortran unit number */
/*      ITEST   error checking option */

/* ROUTINES called */
/*      ERRMES */

/*      SUBROUTINE CPRTMT(A,IA,JA,M,N,NOUT,ITEST) */
/* ***********************************************************************
 */

    /* Parameter adjustments */
    a_dim2 = *ia;
    a_offset = (a_dim2 + 1 << 1) + 1;
    a -= a_offset;

    /* Function Body */

/*     Parameter checking */

    if (*itest != -1) {
	ierror = 0;
	if (*ia < *m || *ja < *n) {
	    ierror = 2;
	}
	if (*n <= 0 || *m <= 0) {
	    ierror = 1;
	}
	*itest = errmes_(itest, &ierror, srname, 6L);
	if (*itest != 0) {
	    return 0;
	}
    }

/*     Main loops */

    io___3.ciunit = *nout;
    s_wsfe(&io___3);
    i__1 = *n;
    for (k = 1; k <= i__1; ++k) {
	i__2 = *m;
	for (j = 1; j <= i__2; ++j) {
	    for (i = 1; i <= 2; ++i) {
		do_fio(&c__1, (char *)&a[i + (j + k * a_dim2 << 1)], (ftnlen)
			sizeof(doublereal));
	    }
	}
    }
    e_wsfe();

} /* cprtmt_ */

