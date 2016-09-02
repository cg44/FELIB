/* fredif.f -- translated by f2c (version of 28 August 1991  0:07:02).
   You must link the resulting object file with the libraries:
	-lF77 -lI77 -lm -lc   (in that order)
*/

#include "f2c.h"


/* Subroutine */ int fredif_(iele, eltop, ieltop, jeltop, nf, inf, jnf, 
	dofnod, first, dif, itest)
integer *iele, *eltop, *ieltop, *jeltop, *nf, *inf, *jnf, *dofnod;
logical *first;
integer *dif, *itest;
{
    /* Initialized data */

    static char srname[6+1] = "FREDIF";

    /* System generated locals */
    integer eltop_dim1, eltop_offset, nf_dim1, nf_offset, i__1, i__2;

    /* Local variables */
    static integer ideg, inod, i, j, nodel, jtest;
    extern integer errmes_(), maxint_();
    static integer ierror, minmum, maxmum;

/* -----------------------------------------------------------------------
 */
/* PURPOSE */
/*      FREDIF calculates the maximum freedom number difference for an */
/*      element */

/* HISTORY */

/*      Copyright (C) 1996 : CCLRC, Rutherford Appleton Laboratory */
/*                           Chilton, Didcot, Oxfordshire OX11 0QX */

/*      Release 1.1  29 Oct 1979 (CG) */
/*      Commented    18 Feb 1980 (KR) */
/*      Release 4.0   2 Oct 1996 (CG) */

/* ARGUMENTS in */
/*      IELE    element number */
/*      ELTOP   ELTOP(I, 1) = element type of element I */
/*              ELTOP(I, 2) = number of nodes on element I */
/*              ELTOP(I, J+2), J=1(1)NUMBER of nodes on element, */
/*              contains the nodes associated with element I */
/*      IELTOP  FIRST dimension of array ELTOP (.GE. IELE) */
/*      JELTOP  second dimension of ELTOP (.GE. number of nodes */
/*              on element) */
/*      NF      NF(I, J) contains the freedom numbers associated */
/*              with node I */
/*      INF     FIRST dimension of NF (.GE. maximum node number */
/*              on element) */
/*      JNF     second dimension of NF (.GE. DOFNOD) */
/*      DOFNOD  number of degrees of freedom per node on the */
/*              element */
/*      FIRST   must be set to .true. for the FIRST call to */
/*              FREDIF and .false. for subsequent calls */
/*      DIF     must be zero for FIRST call to FREDIF */
/*              subsequently contains the maximum freedom */
/*              difference prior to the current call */
/*      ITEST   error checking option */

/* ARGUMENTS out */
/*      FIRST   set to .FALSE. */
/*      DIF     maximum freedom difference for all elements up */
/*              to and including element nele */

/* ROUTINES called */
/*      ERRMES MAXINT */

/*      SUBROUTINE FREDIF(IELE,ELTOP,IELTOP,JELTOP,NF,INF,JNF,DOFNOD, */
/*     *                  FIRST,DIF,ITEST) */
/* ***********************************************************************
 */



    /* Parameter adjustments */
    nf_dim1 = *inf;
    nf_offset = nf_dim1 + 1;
    nf -= nf_offset;
    eltop_dim1 = *ieltop;
    eltop_offset = eltop_dim1 + 1;
    eltop -= eltop_offset;

    /* Function Body */

/*     Initialisation */

    if (*first) {
	*first = FALSE_;
	*dif = 0;
    }

/*     Parameter checking */

    jtest = *itest;
    if (jtest != -1) {
	ierror = 0;
	if (*ieltop < *iele) {
	    ierror = 3;
	}
	if (*jnf < *dofnod) {
	    ierror = 2;
	}
	if (*iele <= 0 || *dofnod <= 0) {
	    ierror = 1;
	}
	*itest = errmes_(&jtest, &ierror, srname, 6L);
	if (*itest != 0) {
	    return 0;
	}
    }

/*     Main body */

    nodel = eltop[*iele + (eltop_dim1 << 1)];

/*     Range checking on NODEL */

    if (jtest != -1) {
	ierror = 0;
	if (*jeltop < nodel + 2) {
	    ierror = 4;
	}
	*itest = errmes_(&jtest, &ierror, srname, 6L);
	if (*itest != 0) {
	    return 0;
	}
    }

    maxmum = 0;
    minmum = maxint_(&maxmum);
    i__1 = nodel;
    for (i = 1; i <= i__1; ++i) {
	inod = eltop[*iele + (i + 2) * eltop_dim1];

/*     Range checking on INOD */

	if (jtest != -1) {
	    ierror = 0;
	    if (*inf < inod) {
		ierror = 5;
	    }
	    *itest = errmes_(&jtest, &ierror, srname, 6L);
	    if (*itest != 0) {
		return 0;
	    }
	}

	i__2 = *dofnod;
	for (j = 1; j <= i__2; ++j) {
	    ideg = nf[inod + j * nf_dim1];
	    if (ideg > 0) {
		maxmum = max(ideg,maxmum);
		minmum = min(ideg,minmum);
	    }
/* L1000: */
	}
/* L1010: */
    }

/*     Calculate maximum freedom number difference */

/* Computing MAX */
    i__1 = *dif, i__2 = maxmum - minmum;
    *dif = max(i__1,i__2);

} /* fredif_ */

