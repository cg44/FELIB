/* prtvec.f -- translated by f2c (version of 28 August 1991  0:07:02).
   You must link the resulting object file with the libraries:
	-lF77 -lI77 -lm -lc   (in that order)
*/

#include "f2c.h"

/* Table of constant values */

static integer c__1 = 1;


/* Subroutine */ int prtvec_(v, iv, n, nout, itest)
doublereal *v;
integer *iv, *n, *nout, *itest;
{
    /* Initialized data */

    static char srname[6+1] = "PRTVEC";

    /* Format strings */
    static char fmt_9990[] = "(\002 \002,6d12.4)";

    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    integer s_wsfe(), do_fio(), e_wsfe();

    /* Local variables */
    static integer i;
    extern integer errmes_();
    static integer ierror;

    /* Fortran I/O blocks */
    static cilist io___3 = { 0, 0, 0, fmt_9990, 0 };


/* -----------------------------------------------------------------------
 */
/* PURPOSE */
/*      PRTVEC prints a vector in a standard format */

/* HISTORY */

/*      Copyright (C) 1996 : CCLRC, Rutherford Appleton Laboratory */
/*                           Chilton, Didcot, Oxfordshire OX11 0QX */

/*      Release 1.1  29 Oct 1979 (CG) */
/*      Commented    14 Oct 1980 (KR) */
/*      Release 4.0   2 Oct 1996 (CG) */

/* ARGUMENTS in */
/*      V       vector of dimension IV containing numbers to be */
/*              printed */
/*      IV      dimension of V (.GE. N) */
/*      N       number of elements of V to be printed */
/*      NOUT    fortran unit number */
/*      ITEST   error checking option */

/* ROUTINES called */
/*      ERRMES */

/*      SUBROUTINE PRTVEC(V,IV,N,NOUT,ITEST) */
/* ***********************************************************************
 */

    /* Parameter adjustments */
    --v;

    /* Function Body */

/*     Parameter checking */

    if (*itest != -1) {
	ierror = 0;
	if (*iv < *n) {
	    ierror = 2;
	}
	if (*n <= 0) {
	    ierror = 1;
	}
	*itest = errmes_(itest, &ierror, srname, 6L);
	if (*itest != 0) {
	    return 0;
	}
    }

/*     Main body */

    io___3.ciunit = *nout;
    s_wsfe(&io___3);
    i__1 = *n;
    for (i = 1; i <= i__1; ++i) {
	do_fio(&c__1, (char *)&v[i], (ftnlen)sizeof(doublereal));
    }
    e_wsfe();

} /* prtvec_ */

