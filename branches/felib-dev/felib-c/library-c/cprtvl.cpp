/* cprtvl.f -- translated by f2c (version of 28 August 1991  0:07:02).
   You must link the resulting object file with the libraries:
	-lF77 -lI77 -lm -lc   (in that order)
*/

#include "f2c.h"

/* Table of constant values */

static integer c__1 = 1;


/* Subroutine */ int cprtvl_(val, ival, nf, inf, jnf, dofnod, totnod, nout, 
	itest)
doublereal *val;
integer *ival, *nf, *inf, *jnf, *dofnod, *totnod, *nout, *itest;
{
    /* Initialized data */

    static char srname[6+1] = "CPRTVL";

    /* Format strings */
    static char fmt_9920[] = "(/\002 \002,2(\002NODE\002,6x,\002REAL\002,7x\
,\002IMAGINARY\002,6x),/\002 \002)";
    static char fmt_9910[] = "(2(\002 \002,i4,2x,2(d12.5,2x),1x))";
    static char fmt_9990[] = "(/\002 \002,\002NODE\002,6x,2(\002REAL\002,7x\
,\002IMAGINARY\002,16x),/\002 \002)";
    static char fmt_9970[] = "(\002 \002,i4,2x,2(d12.5,2x),8x,2(d12.5,2x))";
    static char fmt_9960[] = "(\002 \002,i4,2x,2(d12.5,2x),8x,2(d12.5,2x),\
1(/\002 \002,6x,2(d12.5,2x),8x,2(d12.5,2x)))";
    static char fmt_9950[] = "(\002 \002,i4,2x,2(d12.5,2x),8x,2(d12.5,2x),\
2(/\002 \002,6x,2(d12.5,2x),8x,2(d12.5,2x)))";
    static char fmt_9940[] = "(\002 \002,i4,2x,2(d12.5,2x),8x,2(d12.5,2x),\
3(/\002 \002,6x,2(d12.5,2x),8x,2(d12.5,2x)))";
    static char fmt_9980[] = "(\002 \002,i4,2x,2(d12.5,2x),8x,2(d12.5,2x),\
8(/\002 \002,6x,2(d12.5,2x),8x,2(d12.5,2x)))";
    static char fmt_9930[] = "(\002 \002,i4,2x,2(d12.5,2x),8x,2(d12.5,2x),\
4(/\002 \002,6x,2(d12.5,2x),8x,2(d12.5,2x)))";

    /* System generated locals */
    integer nf_dim1, nf_offset, i__1, i__2;

    /* Builtin functions */
    integer s_wsfe(), e_wsfe(), do_fio();

    /* Local variables */
    static integer node[2];
    static doublereal work[20]	/* was [2][10] */;
    static integer i, j, k, l, n;
    extern integer errmes_();
    static integer ierror;

    /* Fortran I/O blocks */
    static cilist io___3 = { 0, 0, 0, fmt_9920, 0 };
    static cilist io___9 = { 0, 0, 0, fmt_9910, 0 };
    static cilist io___11 = { 0, 0, 0, fmt_9910, 0 };
    static cilist io___12 = { 0, 0, 0, fmt_9990, 0 };
    static cilist io___13 = { 0, 0, 0, fmt_9970, 0 };
    static cilist io___15 = { 0, 0, 0, fmt_9960, 0 };
    static cilist io___16 = { 0, 0, 0, fmt_9950, 0 };
    static cilist io___17 = { 0, 0, 0, fmt_9940, 0 };
    static cilist io___18 = { 0, 0, 0, fmt_9980, 0 };
    static cilist io___19 = { 0, 0, 0, fmt_9930, 0 };


/* -----------------------------------------------------------------------
 */
/* PURPOSE */
/*      CPRTVL prints out the nodal values of the solution in a */
/*      standard format */

/* HISTORY */

/*      Copyright (C) 1996 : CCLRC, Rutherford Appleton Laboratory */
/*                           Chilton, Didcot, Oxfordshire OX11 0QX */

/*      Release 3.0  20 Jan 1984 (CRIE) */
/*      Commented     1 Nov 1985 (CG) */
/*      Release 4.0   2 Oct 1996 (CG) */

/* ARGUMENTS in */
/*      VAL     complex vector of dimension (2*IVAL) containing */
/*              solution values,and prescribed boundary values */
/*      IVAL    dimension of VAL (.GE. TOTAL number of freedoms */
/*              in system) */
/*      NF      integer array of dimension (INF, JNF) containing */
/*              freedom numbers associated with each NODE */
/*      INF     first dimension of NF (.GE. TOTNOD) */
/*      JNF     second dimension of NF (.GE. DOFNOD) */
/*      DOFNOD  number of degrees of freedom at each NODE */
/*      TOTNOD  total number of nodes in mesh */
/*      NOUT    fortran unit number */
/*      ITEST   error checking option */

/* ROUTINES called */
/*      ERRMES */

/*      SUBROUTINE CPRTVL(VAL,IVAL,NF,INF,JNF,DOFNOD,TOTNOD,NOUT,ITEST) */

/* ***********************************************************************
 */

    /* Parameter adjustments */
    nf_dim1 = *inf;
    nf_offset = nf_dim1 + 1;
    nf -= nf_offset;
    val -= 3;

    /* Function Body */

/*     Parameter checking */

    if (*itest != -1) {
	ierror = 0;
	if (*inf < *totnod || *jnf < *dofnod) {
	    ierror = 2;
	}
	if (*dofnod <= 0 || *totnod <= 0 || *nout <= 0) {
	    ierror = 1;
	}
	if (*dofnod > 10) {
	    ierror = 1;
	}
	*itest = errmes_(&ierror, itest, srname, 6L);
	if (*itest != 0) {
	    return 0;
	}
    }

/*     Main loops */

    if (*dofnod == 1) {

/*     For DOFNOD=1 */

	io___3.ciunit = *nout;
	s_wsfe(&io___3);
	e_wsfe();
	n = 0;
	i__1 = *totnod;
	for (i = 1; i <= i__1; ++i) {
	    ++n;
	    node[n - 1] = i;
	    k = nf[i + nf_dim1];
	    work[(n << 1) - 2] = 0.;
	    work[(n << 1) - 1] = 0.;
	    if (k != 0) {
		work[(n << 1) - 2] = val[(k << 1) + 1];
		work[(n << 1) - 1] = val[(k << 1) + 2];
	    }
	    if (n == 2) {
		io___9.ciunit = *nout;
		s_wsfe(&io___9);
		for (j = 1; j <= 2; ++j) {
		    do_fio(&c__1, (char *)&node[j - 1], (ftnlen)sizeof(
			    integer));
		    do_fio(&c__1, (char *)&work[(j << 1) - 2], (ftnlen)sizeof(
			    doublereal));
		    do_fio(&c__1, (char *)&work[(j << 1) - 1], (ftnlen)sizeof(
			    doublereal));
		}
		e_wsfe();
		n = 0;
	    }
/* L1000: */
	}
	if (n != 0) {
	    io___11.ciunit = *nout;
	    s_wsfe(&io___11);
	    do_fio(&c__1, (char *)&node[0], (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&work[0], (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&work[1], (ftnlen)sizeof(doublereal));
	    e_wsfe();
	}
    } else {
	io___12.ciunit = *nout;
	s_wsfe(&io___12);
	e_wsfe();
	i__1 = *totnod;
	for (i = 1; i <= i__1; ++i) {
	    i__2 = *dofnod;
	    for (j = 1; j <= i__2; ++j) {
		k = nf[i + j * nf_dim1];
		work[(j << 1) - 2] = 0.;
		work[(j << 1) - 1] = 0.;
		if (k != 0) {
		    work[(j << 1) - 2] = val[(k << 1) + 1];
		    work[(j << 1) - 1] = val[(k << 1) + 2];
		}
/* L1010: */
	    }

/*     SELECT number of degrees of freedom */

	    switch ((int)*dofnod) {
		case 1:  goto L1080;
		case 2:  goto L1020;
		case 3:  goto L1060;
		case 4:  goto L1030;
		case 5:  goto L1060;
		case 6:  goto L1040;
		case 7:  goto L1060;
		case 8:  goto L1050;
		case 9:  goto L1060;
		case 10:  goto L1070;
	    }
	    goto L1060;
L1020:
	    io___13.ciunit = *nout;
	    s_wsfe(&io___13);
	    do_fio(&c__1, (char *)&i, (ftnlen)sizeof(integer));
	    i__2 = *dofnod;
	    for (l = 1; l <= i__2; ++l) {
		do_fio(&c__1, (char *)&work[(l << 1) - 2], (ftnlen)sizeof(
			doublereal));
		do_fio(&c__1, (char *)&work[(l << 1) - 1], (ftnlen)sizeof(
			doublereal));
	    }
	    e_wsfe();
	    goto L1080;
L1030:
	    io___15.ciunit = *nout;
	    s_wsfe(&io___15);
	    do_fio(&c__1, (char *)&i, (ftnlen)sizeof(integer));
	    i__2 = *dofnod;
	    for (l = 1; l <= i__2; ++l) {
		do_fio(&c__1, (char *)&work[(l << 1) - 2], (ftnlen)sizeof(
			doublereal));
		do_fio(&c__1, (char *)&work[(l << 1) - 1], (ftnlen)sizeof(
			doublereal));
	    }
	    e_wsfe();
	    goto L1080;
L1040:
	    io___16.ciunit = *nout;
	    s_wsfe(&io___16);
	    do_fio(&c__1, (char *)&i, (ftnlen)sizeof(integer));
	    i__2 = *dofnod;
	    for (l = 1; l <= i__2; ++l) {
		do_fio(&c__1, (char *)&work[(l << 1) - 2], (ftnlen)sizeof(
			doublereal));
		do_fio(&c__1, (char *)&work[(l << 1) - 1], (ftnlen)sizeof(
			doublereal));
	    }
	    e_wsfe();
	    goto L1080;
L1050:
	    io___17.ciunit = *nout;
	    s_wsfe(&io___17);
	    do_fio(&c__1, (char *)&i, (ftnlen)sizeof(integer));
	    i__2 = *dofnod;
	    for (l = 1; l <= i__2; ++l) {
		do_fio(&c__1, (char *)&work[(l << 1) - 2], (ftnlen)sizeof(
			doublereal));
		do_fio(&c__1, (char *)&work[(l << 1) - 1], (ftnlen)sizeof(
			doublereal));
	    }
	    e_wsfe();
	    goto L1080;
L1060:
	    io___18.ciunit = *nout;
	    s_wsfe(&io___18);
	    do_fio(&c__1, (char *)&i, (ftnlen)sizeof(integer));
	    i__2 = *dofnod;
	    for (l = 1; l <= i__2; ++l) {
		do_fio(&c__1, (char *)&work[(l << 1) - 2], (ftnlen)sizeof(
			doublereal));
		do_fio(&c__1, (char *)&work[(l << 1) - 1], (ftnlen)sizeof(
			doublereal));
	    }
	    e_wsfe();
	    goto L1080;
L1070:
	    io___19.ciunit = *nout;
	    s_wsfe(&io___19);
	    do_fio(&c__1, (char *)&i, (ftnlen)sizeof(integer));
	    i__2 = *dofnod;
	    for (l = 1; l <= i__2; ++l) {
		do_fio(&c__1, (char *)&work[(l << 1) - 2], (ftnlen)sizeof(
			doublereal));
		do_fio(&c__1, (char *)&work[(l << 1) - 1], (ftnlen)sizeof(
			doublereal));
	    }
	    e_wsfe();
L1080:
	    ;
	}
    }

} /* cprtvl_ */

