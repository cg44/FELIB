/* cprtvc.f -- translated by f2c (version of 28 August 1991  0:07:02).
   You must link the resulting object file with the libraries:
	-lF77 -lI77 -lm -lc   (in that order)
*/

#include "f2c.h"

/* Table of constant values */

static integer c__1 = 1;


/* Subroutine */ int cprtvc_(v, iv, n, nout, itest)
doublereal *v;
integer *iv, *n, *nout, *itest;
{
    /* Initialized data */

    static char srname[6+1] = "CPRTVC";

    /* Format strings */
    static char fmt_9990[] = "(\002 \002,3(2x,2d12.4))";

    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    integer s_wsfe(), do_fio(), e_wsfe();

    /* Local variables */
    static integer i, j;
    extern integer errmes_();
    static integer ierror;

    /* Fortran I/O blocks */
    static cilist io___3 = { 0, 0, 0, fmt_9990, 0 };


/* -----------------------------------------------------------------------
 */
/* PURPOSE */
/*      CPRTVC prints a complex vector in a standard format */

/* HISTORY */

/*      Copyright (C) 1996 : CCLRC, Rutherford Appleton Laboratory */
/*                           Chilton, Didcot, Oxfordshire OX11 0QX */

/*      Release 3.0  29 Oct 1984 (CRIE) */
/*      comments      1 Nov 1985 (CG) */
/*      Release 4.0   2 Oct 1996 (CG) */

/* ARGUMENTS in */
/*      V       vector of dimension (2, IV) containing numbers to be */
/*              printed */
/*      IV      dimension of V (.GE. N) */
/*      N       number of elements of V to be printed */
/*      NOUT    fortran unit number */
/*      ITEST   error checking option */

/* ROUTINES called */
/*      ERRMES */


/*      SUBROUTINE CPRTVC(V,IV,N,NOUT,ITEST) */
/* ***********************************************************************
 */

    /* Parameter adjustments */
    v -= 3;

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

/*     Main loops */

    io___3.ciunit = *nout;
    s_wsfe(&io___3);
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	for (i = 1; i <= 2; ++i) {
	    do_fio(&c__1, (char *)&v[i + (j << 1)], (ftnlen)sizeof(doublereal)
		    );
	}
    }
    e_wsfe();

} /* cprtvc_ */

