/* prtgeo.f -- translated by f2c (version of 28 August 1991  0:07:02).
   You must link the resulting object file with the libraries:
	-lF77 -lI77 -lm -lc   (in that order)
*/

#include "f2c.h"

/* Table of constant values */

static integer c__1 = 1;


/* Subroutine */ int prtgeo_(totnod, dimen, coord, icoord, jcoord, nout, 
	itest)
integer *totnod, *dimen;
doublereal *coord;
integer *icoord, *jcoord, *nout, *itest;
{
    /* Initialized data */

    static char srname[6+1] = "PRTGEO";

    /* Format strings */
    static char fmt_9990[] = "(\002 \002,////\002 **** NODAL GEOMETRY ***\
*\002,//\002 \002)";
    static char fmt_9980[] = "(\002 NUMBER OF NODES = \002,i3)";
    static char fmt_9970[] = "(\002 DIMENSIONS      = \002,i3)";
    static char fmt_9960[] = "(/\002 \002,2x,\002NODE\002,9x,\002X1\002,/\
\002 \002)";
    static char fmt_9950[] = "(/\002 \002,2x,\002NODE\002,9x,\002X1\002,12x\
,\002X2\002,/\002 \002)";
    static char fmt_9940[] = "(/\002 \002,2x,\002NODE\002,9x,\002X1\002,12x\
,\002X2\002,12x,\002X3\002,/\002 \002)";
    static char fmt_9930[] = "(/\002 \002,2x,\002NODE\002,9x,\002X1\002,12x\
,\002X2\002,12x,\002X3\002,12x,\002X4\002,/\002 \002)";
    static char fmt_9920[] = "(\002 \002,2x,i3,5x,4(d12.4,2x))";

    /* System generated locals */
    integer coord_dim1, coord_offset, i__1, i__2;

    /* Builtin functions */
    integer s_wsfe(), e_wsfe(), do_fio();

    /* Local variables */
    static integer i, j;
    extern integer errmes_();
    static integer ierror;

    /* Fortran I/O blocks */
    static cilist io___3 = { 0, 0, 0, fmt_9990, 0 };
    static cilist io___4 = { 0, 0, 0, fmt_9980, 0 };
    static cilist io___5 = { 0, 0, 0, fmt_9970, 0 };
    static cilist io___6 = { 0, 0, 0, fmt_9960, 0 };
    static cilist io___7 = { 0, 0, 0, fmt_9950, 0 };
    static cilist io___8 = { 0, 0, 0, fmt_9940, 0 };
    static cilist io___9 = { 0, 0, 0, fmt_9930, 0 };
    static cilist io___11 = { 0, 0, 0, fmt_9920, 0 };


/* -----------------------------------------------------------------------
 */
/* PURPOSE */
/*      PRTGEO prints out element geometry in a standard format */

/* HISTORY */

/*      Copyright (C) 1996 : CCLRC, Rutherford Appleton Laboratory */
/*                           Chilton, Didcot, Oxfordshire OX11 0QX */

/*      Release 1.1  29 Oct 1979 (CG) */
/*      Commented    14 Oct 1980 (KR) */
/*      Release 4.0   2 Oct 1996 (CG) */

/* ARGUMENTS in */
/*      TOTNOD  total number of nodes in the mesh */
/*      DIMEN   dimensionality of the geometric data */
/*      COORD   array of dimension (ICOORD, JCOORD) containing */
/*              global coordinates of the nodes */
/*      ICOORD  first dimension of COORD (.GE. TOTNOD) */
/*      JCOORD  second dimension of COORD (.GE. DIMEN) */
/*      NOUT    fortran unit number */
/*      ITEST   error checking option */

/* ROUTINES called */
/*      ERRMES */

/*      SUBROUTINE PRTGEO(TOTNOD,DIMEN,COORD,ICOORD,JCOORD,NOUT,ITEST) */
/* ***********************************************************************
 */

    /* Parameter adjustments */
    coord_dim1 = *icoord;
    coord_offset = coord_dim1 + 1;
    coord -= coord_offset;

    /* Function Body */

/*     Parameter checking */

    if (*itest != -1) {
	ierror = 0;
	if (*icoord < *totnod || *jcoord < *dimen) {
	    ierror = 2;
	}
	if (*totnod <= 0 || *dimen <= 0) {
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
    e_wsfe();
    io___4.ciunit = *nout;
    s_wsfe(&io___4);
    do_fio(&c__1, (char *)&(*totnod), (ftnlen)sizeof(integer));
    e_wsfe();
    io___5.ciunit = *nout;
    s_wsfe(&io___5);
    do_fio(&c__1, (char *)&(*dimen), (ftnlen)sizeof(integer));
    e_wsfe();
    if (*dimen == 1) {
	io___6.ciunit = *nout;
	s_wsfe(&io___6);
	e_wsfe();
    }
    if (*dimen == 2) {
	io___7.ciunit = *nout;
	s_wsfe(&io___7);
	e_wsfe();
    }
    if (*dimen == 3) {
	io___8.ciunit = *nout;
	s_wsfe(&io___8);
	e_wsfe();
    }
    if (*dimen == 4) {
	io___9.ciunit = *nout;
	s_wsfe(&io___9);
	e_wsfe();
    }

    i__1 = *totnod;
    for (i = 1; i <= i__1; ++i) {
	io___11.ciunit = *nout;
	s_wsfe(&io___11);
	do_fio(&c__1, (char *)&i, (ftnlen)sizeof(integer));
	i__2 = *dimen;
	for (j = 1; j <= i__2; ++j) {
	    do_fio(&c__1, (char *)&coord[i + j * coord_dim1], (ftnlen)sizeof(
		    doublereal));
	}
	e_wsfe();
/* L1000: */
    }

} /* prtgeo_ */

