/* prtval.f -- translated by f2c (version of 28 August 1991  0:07:02).
   You must link the resulting object file with the libraries:
	-lF77 -lI77 -lm -lc   (in that order)
*/

#include "f2c.h"

/* Table of constant values */

static integer c__1 = 1;


/* Subroutine */ int prtval_(val, ival, nf, inf, jnf, dofnod, totnod, nout, 
	itest)
doublereal *val;
integer *ival, *nf, *inf, *jnf, *dofnod, *totnod, *nout, *itest;
{
    /* Initialized data */

    static char srname[6+1] = "PRTVAL";

    /* Format strings */
    static char fmt_9950[] = "(/\002 \002,4(\002NODE\002,5x,\002VALUE\002,5x\
),/\002 \002)";
    static char fmt_9930[] = "(\002 \002,4(i4,1x,d12.5,2x))";
    static char fmt_9990[] = "(/\002 \002,2(\002NODE\002,6x,2(\002VALUE\002,\
9x)),/\002 \002)";
    static char fmt_9980[] = "(/\002 NODE\002,6x,3(\002VALUE\002,9x),/\002\
 \002)";
    static char fmt_9970[] = "(/\002 NODE\002,6x,4(\002VALUE\002,9x),/\002\
 \002)";
    static char fmt_9960[] = "(/\002 NODE\002,6x,5(\002VALUE\002,9x),/\002\
 \002)";
    static char fmt_9920[] = "(\002 \002,2(i4,1x,2(d12.5,2x),5x))";
    static char fmt_9940[] = "(\002 \002,i4,1x,5(d12.5,2x),/\002 \002,5x,5(d\
12.5,2x),/\002 \002,5x,5(d12.5,2x))";

    /* System generated locals */
    integer nf_dim1, nf_offset, i__1, i__2, i__3;

    /* Builtin functions */
    integer s_wsfe(), e_wsfe(), do_fio();

    /* Local variables */
    static integer node[5];
    static doublereal work[5];
    static integer i, j, k, n, jtest;
    extern integer errmes_();
    static integer ierror;

    /* Fortran I/O blocks */
    static cilist io___4 = { 0, 0, 0, fmt_9950, 0 };
    static cilist io___10 = { 0, 0, 0, fmt_9930, 0 };
    static cilist io___12 = { 0, 0, 0, fmt_9930, 0 };
    static cilist io___13 = { 0, 0, 0, fmt_9990, 0 };
    static cilist io___14 = { 0, 0, 0, fmt_9980, 0 };
    static cilist io___15 = { 0, 0, 0, fmt_9970, 0 };
    static cilist io___16 = { 0, 0, 0, fmt_9960, 0 };
    static cilist io___17 = { 0, 0, 0, fmt_9920, 0 };
    static cilist io___18 = { 0, 0, 0, fmt_9940, 0 };
    static cilist io___19 = { 0, 0, 0, fmt_9920, 0 };
    static cilist io___20 = { 0, 0, 0, fmt_9940, 0 };


/* -----------------------------------------------------------------------
 */
/* PURPOSE */
/*      PRTVAL prints out the nodal values of the solution in a */
/*      standard format */

/* HISTORY */

/*      Copyright (C) 1996 : CCLRC, Rutherford Appleton Laboratory */
/*                           Chilton, Didcot, Oxfordshire OX11 0QX */

/*      Release 1.1  29 Oct 1979 (CG) */
/*      Commented    14 Oct 1980 (KR) */
/*      Release 4.0   2 Oct 1996 (CG) */

/* ARGUMENTS in */
/*      VAL     vector of dimension IVAL containing solution */
/*              values,and prescribed boundary values */
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

/*     SUBROUTINE PRTVAL(VAL,IVAL,NF,INF,JNF,DOFNOD,TOTNOD,NOUT,ITEST) */
/* ***********************************************************************
 */

    /* Parameter adjustments */
    nf_dim1 = *inf;
    nf_offset = nf_dim1 + 1;
    nf -= nf_offset;
    --val;

    /* Function Body */

/*     Parameter checking */

    jtest = *itest;
    if (jtest != -1) {
	ierror = 0;
	if (*inf < *totnod || *jnf < *dofnod) {
	    ierror = 2;
	}
	if (*dofnod <= 0 || *totnod <= 0) {
	    ierror = 1;
	}
	*itest = errmes_(&jtest, &ierror, srname, 6L);
	if (*itest != 0) {
	    return 0;
	}
    }

/*     Main body */

    if (*dofnod == 1) {

	io___4.ciunit = *nout;
	s_wsfe(&io___4);
	e_wsfe();
	n = 0;
	i__1 = *totnod;
	for (i = 1; i <= i__1; ++i) {
	    ++n;
	    node[n - 1] = i;
	    k = nf[i + nf_dim1];
	    work[n - 1] = 0.;
	    if (k != 0) {
		if (k <= 0) {
		    k = abs(k);
		}
		work[n - 1] = val[k];
	    }
	    if (n == 4) {
		io___10.ciunit = *nout;
		s_wsfe(&io___10);
		for (j = 1; j <= 4; ++j) {
		    do_fio(&c__1, (char *)&node[j - 1], (ftnlen)sizeof(
			    integer));
		    do_fio(&c__1, (char *)&work[j - 1], (ftnlen)sizeof(
			    doublereal));
		}
		e_wsfe();
		n = 0;
	    }
/* L1000: */
	}

	if (n != 0) {
	    io___12.ciunit = *nout;
	    s_wsfe(&io___12);
	    i__1 = n;
	    for (j = 1; j <= i__1; ++j) {
		do_fio(&c__1, (char *)&node[j - 1], (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&work[j - 1], (ftnlen)sizeof(doublereal)
			);
	    }
	    e_wsfe();
	}
    } else {
	if (*dofnod == 2) {
	    io___13.ciunit = *nout;
	    s_wsfe(&io___13);
	    e_wsfe();
	}
	if (*dofnod == 3) {
	    io___14.ciunit = *nout;
	    s_wsfe(&io___14);
	    e_wsfe();
	}
	if (*dofnod == 4) {
	    io___15.ciunit = *nout;
	    s_wsfe(&io___15);
	    e_wsfe();
	}
	if (*dofnod >= 5) {
	    io___16.ciunit = *nout;
	    s_wsfe(&io___16);
	    e_wsfe();
	}
	n = 0;
	i__1 = *totnod;
	for (i = 1; i <= i__1; ++i) {
	    i__2 = *dofnod;
	    for (j = 1; j <= i__2; ++j) {
		++n;
		node[n - 1] = i;
		k = nf[i + j * nf_dim1];

/*     Range checking on K */

		if (jtest != -1) {
		    ierror = 0;
		    if (*ival < k) {
			ierror = 3;
		    }
		    *itest = errmes_(&jtest, &ierror, srname, 6L);
		    if (*itest != 0) {
			return 0;
		    }
		}

		work[n - 1] = 0.;
		if (k != 0) {
		    if (k <= 0) {
			k = abs(k);
		    }
		    work[n - 1] = val[k];
		}
		if ((*dofnod != 2 || n == 4) && *dofnod == 2) {
		    io___17.ciunit = *nout;
		    s_wsfe(&io___17);
		    do_fio(&c__1, (char *)&node[0], (ftnlen)sizeof(integer));
		    do_fio(&c__1, (char *)&work[0], (ftnlen)sizeof(doublereal)
			    );
		    do_fio(&c__1, (char *)&work[1], (ftnlen)sizeof(doublereal)
			    );
		    do_fio(&c__1, (char *)&node[2], (ftnlen)sizeof(integer));
		    do_fio(&c__1, (char *)&work[2], (ftnlen)sizeof(doublereal)
			    );
		    do_fio(&c__1, (char *)&work[3], (ftnlen)sizeof(doublereal)
			    );
		    e_wsfe();
		    n = 0;
		} else if (*dofnod != 2 && n == *dofnod) {
		    io___18.ciunit = *nout;
		    s_wsfe(&io___18);
		    do_fio(&c__1, (char *)&i, (ftnlen)sizeof(integer));
		    i__3 = *dofnod;
		    for (k = 1; k <= i__3; ++k) {
			do_fio(&c__1, (char *)&work[k - 1], (ftnlen)sizeof(
				doublereal));
		    }
		    e_wsfe();
		    n = 0;
		}
/* L1010: */
	    }
/* L1020: */
	}

	if (n != 0) {
	    if (*dofnod == 2) {
		io___19.ciunit = *nout;
		s_wsfe(&io___19);
		do_fio(&c__1, (char *)&node[0], (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&work[0], (ftnlen)sizeof(doublereal));
		do_fio(&c__1, (char *)&work[1], (ftnlen)sizeof(doublereal));
		e_wsfe();
	    }
	    if (*dofnod != 2) {
		io___20.ciunit = *nout;
		s_wsfe(&io___20);
		do_fio(&c__1, (char *)&i, (ftnlen)sizeof(integer));
		i__1 = n;
		for (k = 1; k <= i__1; ++k) {
		    do_fio(&c__1, (char *)&work[k - 1], (ftnlen)sizeof(
			    doublereal));
		}
		e_wsfe();
	    }
	}
    }

} /* prtval_ */

