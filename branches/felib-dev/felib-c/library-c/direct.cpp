/* direct.f -- translated by f2c (version of 28 August 1991  0:07:02).
   You must link the resulting object file with the libraries:
	-lF77 -lI77 -lm -lc   (in that order)
*/

#include "f2c.h"


/* Subroutine */ int direct_(nele, eltop, ieltop, jeltop, nf, inf, jnf, 
	dofnod, steer, isteer, itest)
integer *nele, *eltop, *ieltop, *jeltop, *nf, *inf, *jnf, *dofnod, *steer, *
	isteer, *itest;
{
    /* Initialized data */

    static char srname[6+1] = "DIRECT";

    /* System generated locals */
    integer eltop_dim1, eltop_offset, nf_dim1, nf_offset, i__1, i__2;

    /* Local variables */
    static integer ideg, inod, jnod, k, nodel, jtest;
    extern integer errmes_();
    static integer ierror;

/* -----------------------------------------------------------------------
 */
/* PURPOSE */
/*      DIRECT constructs the steering vector to direct the assembly of a 
*/
/*      system matrix */

/* HISTORY */

/*      Copyright (C) 1996 : CCLRC, Rutherford Appleton Laboratory */
/*                           Chilton, Didcot, Oxfordshire OX11 0QX */

/*      Release 1.1  29 Oct 1979 (CG) */
/*      Commented    13 Feb 1980 (KR) */
/*      Release 4.0   2 Oct 1996 (CG) */

/* ARGUMENTS in */
/*      NELE    element number */
/*      ELTOP   2D array containing element type, number of */
/*              nodes on the element, and the element topologies */
/*      IELTOP  first dimension of ELTOP (.GE. number of */
/*              elements in problem) */
/*      JELTOP  second dimension of ELTOP (.GE. number of nodes */
/*              on element + 2) */
/*      NF      contains freedom numbers associated with each */
/*              node */
/*      INF     first dimension of NF (.GE. total number of */
/*              nodes in problem) */
/*      JNF     second dimension of NF (.GE. DOFNOD) */
/*      DOFNOD  number of degrees of freedom at each node */
/*      ISTEER  dimension of vector STEER (.GE. total number of */
/*              degrees of freedom on element) */
/*      ITEST   error checking option */

/* ARGUMENTS out */
/*      STEER   vector containing freedom numbers associated */
/*              with element NELE, arranged in local node order */

/* ROUTINES called */
/*      ERRMES */

/*      SUBROUTINE DIRECT(NELE,ELTOP,IELTOP,JELTOP,NF,INF,JNF,DOFNOD, */
/*     *                  STEER,ISTEER,ITEST) */
/* ***********************************************************************
 */

    /* Parameter adjustments */
    --steer;
    nf_dim1 = *inf;
    nf_offset = nf_dim1 + 1;
    nf -= nf_offset;
    eltop_dim1 = *ieltop;
    eltop_offset = eltop_dim1 + 1;
    eltop -= eltop_offset;

    /* Function Body */

/*     Parameter checking */

    jtest = *itest;
    if (jtest != -1) {
	ierror = 0;
	if (*jnf < *dofnod) {
	    ierror = 3;
	}
	if (*ieltop < *nele) {
	    ierror = 2;
	}
	if (*nele <= 0 || *dofnod <= 0) {
	    ierror = 1;
	}
	*itest = errmes_(&jtest, &ierror, srname, 6L);
	if (*itest != 0) {
	    return 0;
	}
    }

/*     Main loops */

    nodel = eltop[*nele + (eltop_dim1 << 1)];

/*     Range check on NODEL : number of node in element */

    if (jtest != -1) {
	ierror = 0;
	if (*jeltop < nodel + 2) {
	    ierror = 4;
	}
	if (*isteer < *dofnod * nodel) {
	    ierror = 6;
	}
	*itest = errmes_(&jtest, &ierror, srname, 6L);
	if (*itest != 0) {
	    return 0;
	}
    }

    k = 1;
    i__1 = nodel;
    for (inod = 1; inod <= i__1; ++inod) {
	jnod = eltop[*nele + (inod + 2) * eltop_dim1];

/*     Range check on JNOD : node number */

	if (jtest != -1) {
	    ierror = 0;
	    if (*inf < jnod) {
		ierror = 5;
	    }
	    *itest = errmes_(&jtest, &ierror, srname, 6L);
	    if (*itest != 0) {
		return 0;
	    }
	}

/*     Construct steering vector */

	i__2 = *dofnod;
	for (ideg = 1; ideg <= i__2; ++ideg) {
	    steer[k] = nf[jnod + ideg * nf_dim1];
	    ++k;
/* L1000: */
	}
/* L1010: */
    }

} /* direct_ */

