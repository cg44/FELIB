/* prttop.f -- translated by f2c (version of 28 August 1991  0:07:02).
   You must link the resulting object file with the libraries:
	-lF77 -lI77 -lm -lc   (in that order)
*/

#include "f2c.h"

/* Table of constant values */

static integer c__1 = 1;


/* Subroutine */ int prttop_(totels, eltop, ieltop, jeltop, nout, itest)
integer *totels, *eltop, *ieltop, *jeltop, *nout, *itest;
{
    /* Initialized data */

    static char srname[6+1] = "PRTTOP";

    /* Format strings */
    static char fmt_9990[] = "(\002 \002,////\002 **** ELEMENT TOPOLOGY ***\
*\002,//\002 \002)";
    static char fmt_9980[] = "(\002 NUMBER OF ELEMENTS = \002,i3)";
    static char fmt_9970[] = "(/\002 \002,2x,\002ELEM\002,4x,\002ELTYP\002,4\
x,\002NODEL\002,4x,\002NODES\002,/\002 \002)";
    static char fmt_9960[] = "(\002 \002,2x,i4,4x,i3,7x,i3,4x,9(i4,2x),\
/\002 \002,27x,9(i4,2x))";

    /* System generated locals */
    integer eltop_dim1, eltop_offset, i__1, i__2;

    /* Builtin functions */
    integer s_wsfe(), e_wsfe(), do_fio();

    /* Local variables */
    static integer i, j, k, l, jtest;
    extern integer errmes_();
    static integer ierror;

    /* Fortran I/O blocks */
    static cilist io___4 = { 0, 0, 0, fmt_9990, 0 };
    static cilist io___5 = { 0, 0, 0, fmt_9980, 0 };
    static cilist io___6 = { 0, 0, 0, fmt_9970, 0 };
    static cilist io___10 = { 0, 0, 0, fmt_9960, 0 };


/* -----------------------------------------------------------------------
 */
/* PURPOSE */
/*      PRTTOP prints element topologies in a standard format */

/* HISTORY */

/*      Copyright (C) 1996 : CCLRC, Rutherford Appleton Laboratory */
/*                           Chilton, Didcot, Oxfordshire OX11 0QX */

/*      Release 1.1  29 Oct 1979 (CG) */
/*      Commented    14 Oct 1980 (KR) */
/*      Release 4.0   2 Oct 1996 (CG) */

/* ARGUMENTS in */
/*      TOTELS  total number of elements in the mesh */
/*      ELTOP   integer array of dimension (IELTOP, JELTOP) */
/*              containing element topologies, element type, and */
/*              number of nodes on the element */
/*      IELTOP  first dimension of ELTOP (.GE. TOTELS) */
/*      JELTOP  second dimension of ELTOP (.GE. NUMBER of nodes */
/*              on element + 2) */
/*      NOUT    fortran unit number */
/*      ITEST   error checking option */

/* ROUTINES called */
/*      ERRMES */


/*      SUBROUTINE PRTTOP(TOTELS,ELTOP,IELTOP,JELTOP,NOUT,ITEST) */
/* ***********************************************************************
 */

    /* Parameter adjustments */
    eltop_dim1 = *ieltop;
    eltop_offset = eltop_dim1 + 1;
    eltop -= eltop_offset;

    /* Function Body */

/*     Parameter checking */

    jtest = *itest;
    if (jtest != -1) {
	ierror = 0;
	if (*ieltop < *totels) {
	    ierror = 2;
	}
	if (*totels <= 0) {
	    ierror = 1;
	}
	*itest = errmes_(&jtest, &ierror, srname, 6L);
	if (*itest != 0) {
	    return 0;
	}
    }

/*     Main body */

    io___4.ciunit = *nout;
    s_wsfe(&io___4);
    e_wsfe();
    io___5.ciunit = *nout;
    s_wsfe(&io___5);
    do_fio(&c__1, (char *)&(*totels), (ftnlen)sizeof(integer));
    e_wsfe();
    io___6.ciunit = *nout;
    s_wsfe(&io___6);
    e_wsfe();
    i__1 = *totels;
    for (i = 1; i <= i__1; ++i) {
	l = eltop[i + (eltop_dim1 << 1)];
	k = l + 2;

/*     Range checking on K */

	if (jtest != -1) {
	    ierror = 0;
	    if (*jeltop < k) {
		ierror = 3;
	    }
	    *itest = errmes_(&jtest, &ierror, srname, 6L);
	    if (*itest != 0) {
		return 0;
	    }
	}

	io___10.ciunit = *nout;
	s_wsfe(&io___10);
	do_fio(&c__1, (char *)&i, (ftnlen)sizeof(integer));
	i__2 = k;
	for (j = 1; j <= i__2; ++j) {
	    do_fio(&c__1, (char *)&eltop[i + j * eltop_dim1], (ftnlen)sizeof(
		    integer));
	}
	e_wsfe();
/* L1000: */
    }

} /* prttop_ */

