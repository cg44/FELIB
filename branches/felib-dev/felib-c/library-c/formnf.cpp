/* formnf.f -- translated by f2c (version of 28 August 1991  0:07:02).
   You must link the resulting object file with the libraries:
	-lF77 -lI77 -lm -lc   (in that order)
*/

#include "f2c.h"


/* Subroutine */ int formnf_(rest, irest, jrest, resnod, totnod, dofnod, nf, 
	inf, jnf, totdof, itest)
integer *rest, *irest, *jrest, *resnod, *totnod, *dofnod, *nf, *inf, *jnf, *
	totdof, *itest;
{
    /* Initialized data */

    static char srname[6+1] = "FORMNF";

    /* System generated locals */
    integer nf_dim1, nf_offset, rest_dim1, rest_offset, i__1, i__2;

    /* Local variables */
    static integer i, j, k, l, m, jtest;
    extern integer errmes_();
    static logical switch_;
    static integer ierror;

/* -----------------------------------------------------------------------
 */
/* PURPOSE */
/*      FORMNF constructs the nodal freedom array from the restrained */
/*      freedom data */

/* HISTORY */

/*      Copyright (C) 1996 : CCLRC, Rutherford Appleton Laboratory */
/*                           Chilton, Didcot, Oxfordshire OX11 0QX */

/*      Release 1.1  29 Oct 1979 (IMS) */
/*      Commented    18 Feb 1980 (KR) */
/*      Recoded      01 Nov 1981 (NB) */
/*      Release 4.0   2 Oct 1996 (CG) */

/* ARGUMENTS in */
/*      REST    integer array; REST(I ,J) contains the I'TH set */
/*              of restraint information - REST(I, 1) contains */
/*              the node number, REST(I, J+1) for J=1(1)DOFNOD */
/*              contains the local freedom numbers of the */
/*              freedoms that are restrained */
/*      IREST   first dimension of array REST (.GE. RESNOD) */
/*      JREST   second dimension of REST (.GE. DOFNOD) */
/*      RESNOD  number of nodes at which freedoms are restrained */
/*      TOTNOD  total number of nodes in mesh */
/*      INF     first dimension of array INF (.GE. TOTNOD) */
/*      JNF     second dimension of INF (.GE. DOFNOD) */
/*      ITEST   error checking option */

/* ARGUMENTS out */
/*      NF      NF(I, J), J=1(1)DOFNOD, contains the freedom */
/*              numbers associated with the I'TH node */
/*      TOTDOF  total number of freedoms in problem under */
/*              consideration */

/* ROUTINES called */
/*      ERRMES */

/*      SUBROUTINE FORMNF(REST,IREST,JREST,RESNOD,TOTNOD,DOFNOD,NF,INF, */

/*     *                  JNF,TOTDOF,ITEST) */
/* ***********************************************************************
 */



    /* Parameter adjustments */
    nf_dim1 = *inf;
    nf_offset = nf_dim1 + 1;
    nf -= nf_offset;
    rest_dim1 = *irest;
    rest_offset = rest_dim1 + 1;
    rest -= rest_offset;

    /* Function Body */

/*     Parameter checking */

    jtest = *itest;
    if (jtest != -1) {
	ierror = 0;
	if (*inf < *totnod || *jnf < *dofnod) {
	    ierror = 3;
	}
	if (*irest < *resnod || *jrest < *dofnod + 1) {
	    ierror = 2;
	}
	if (*resnod < 0 || *totnod <= 0 || *dofnod <= 0) {
	    ierror = 1;
	}
	*itest = errmes_(itest, &ierror, srname, 6L);
	if (*itest != 0) {
	    return 0;
	}
    }

/*     Main body */

    switch_ = TRUE_;

/*     Initialise nodel freedom array */

    i__1 = *totnod;
    for (i = 1; i <= i__1; ++i) {
	i__2 = *dofnod;
	for (j = 1; j <= i__2; ++j) {
	    nf[i + j * nf_dim1] = 1;
/* L1000: */
	}
/* L1010: */
    }

/*     If no restrained nodes branch */

    if (*resnod != 0) {
	i__1 = *resnod;
	for (i = 1; i <= i__1; ++i) {
	    k = rest[i + rest_dim1];
	    i__2 = *dofnod;
	    for (j = 1; j <= i__2; ++j) {
		l = rest[i + (j + 1) * rest_dim1];
		m = abs(l);

/*     Range checking on K and M */

		if (jtest != -1) {
		    ierror = 0;
		    if (k > *totnod || m > *dofnod) {
			ierror = 4;
		    }
		    *itest = errmes_(&jtest, &ierror, srname, 6L);
		    if (*itest != 0) {
			return 0;
		    }
		}

		if (l > 0) {
		    nf[k + l * nf_dim1] = 0;
		}
		if (l < 0) {
		    nf[k + m * nf_dim1] = l;
		    switch_ = FALSE_;
		}
/* L1020: */
	    }
/* L1030: */
	}
    }

/*     Renumber nodal freedom array */

    k = 1;
    i__1 = *totnod;
    for (i = 1; i <= i__1; ++i) {
	i__2 = *dofnod;
	for (j = 1; j <= i__2; ++j) {
	    if (nf[i + j * nf_dim1] > 0) {
		nf[i + j * nf_dim1] = k;
		++k;
	    }
/* L1040: */
	}
/* L1050: */
    }

/*     Set up for prescribed values */

    *totdof = k - 1;
    if (! switch_) {
	i__1 = *totnod;
	for (i = 1; i <= i__1; ++i) {
	    i__2 = *dofnod;
	    for (j = 1; j <= i__2; ++j) {
		if (nf[i + j * nf_dim1] < 0) {
		    nf[i + j * nf_dim1] = -k;
		    ++k;
		}
/* L1060: */
	    }
/* L1070: */
	}
    }
} /* formnf_ */

