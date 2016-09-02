/* bndwth.f -- translated by f2c (version of 28 August 1991  0:07:02).
   You must link the resulting object file with the libraries:
	-lF77 -lI77 -lm -lc   (in that order)
*/

#include "f2c.h"


/* Subroutine */ int bndwth_(eltop, ieltop, jeltop, nf, inf, jnf, dofnod, 
	totels, hband, itest)
integer *eltop, *ieltop, *jeltop, *nf, *inf, *jnf, *dofnod, *totels, *hband, *
	itest;
{
    /* Initialized data */

    static char srname[6+1] = "BNDWTH";

    /* System generated locals */
    integer eltop_dim1, eltop_offset, nf_dim1, nf_offset, i__1, i__2, i__3;

    /* Local variables */
    static integer ideg, iele, inod, i, j, nodel, jtest;
    extern integer errmes_(), maxint_();
    static integer ierror, min_, max_;

/* -----------------------------------------------------------------------
 */
/* PURPOSE */
/*      BNDWTH calculates the maximum freedom number difference over */
/*      all the elements in the mesh. */

/* HISTORY */

/*      Copyright (C) 1996 : CCLRC, Rutherford Appleton Laboratory */
/*                           Chilton, Didcot, Oxfordshire OX11 0QX */

/*      Release 3.0  29 Oct 1984 (CG) */
/*      Release 4.0   2 Oct 1996 (CG) */

/* ARGUMENTS in */
/*      ELTOP   ELTOP(I,1) = element type of element I */
/*              ELTOP(I,2) = number of nodes on element I */
/*              ELTOP(I,J+2), J=1(1)NUMBER of nodes on element, */
/*              contains the nodes associated with element I */
/*      IELTOP  first dimension of array ELTOP (.GE. TOTELS) */
/*      JELTOP  second dimension of ELTOP (.GE. number of nodes */
/*              on element) */
/*      NF      NF(I,J) contains the freedom numbers associated */
/*              with node I */
/*      INF     first dimension of NF (.GE. maximum node number */
/*              on element) */
/*      JNF     second dimension of NF (.GE. DOFNOD) */
/*      DOFNOD  number of degrees of freedom per node on the */
/*              element */
/*      TOTELS  number of elements in the mesh */
/*      ITEST   error checking option */

/* ARGUMENTS out */
/*      HBAND   semi bandwidth */

/* ROUTINES called */
/*      ERRMES    MAXINT */

/*      SUBROUTINE BNDWTH(ELTOP,IELTOP,JELTOP,NF,INF,JNF,DOFNOD,TOTELS, */

/*     *                  HBAND,ITEST) */
/* ********************************************************* */



    /* Parameter adjustments */
    nf_dim1 = *inf;
    nf_offset = nf_dim1 + 1;
    nf -= nf_offset;
    eltop_dim1 = *ieltop;
    eltop_offset = eltop_dim1 + 1;
    eltop -= eltop_offset;

    /* Function Body */

/*     Initialisation */

    *hband = 0;
    jtest = *itest;

/*     Parameter checking */

    if (jtest != -1) {
	ierror = 0;
	if (*dofnod <= 0) {
	    ierror = 1;
	}
	if (*jnf < *dofnod) {
	    ierror = 2;
	}
	if (*totels > *ieltop) {
	    ierror = 3;
	}
	*itest = errmes_(&jtest, &ierror, srname, 6L);
	if (*itest != 0) {
	    return 0;
	}
    }

/*     Main scanning loops */

    i__1 = *totels;
    for (iele = 1; iele <= i__1; ++iele) {
	max_ = 0;
	min_ = maxint_(&max_);
	nodel = eltop[iele + (eltop_dim1 << 1)];

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

/*     Loop over nodes */

	i__2 = nodel;
	for (i = 1; i <= i__2; ++i) {
	    inod = eltop[iele + (i + 2) * eltop_dim1];

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

/*     Loop over freedoms */

	    i__3 = *dofnod;
	    for (j = 1; j <= i__3; ++j) {
		ideg = nf[inod + j * nf_dim1];
		if (ideg != 0) {
		    max_ = max(ideg,max_);
		    min_ = min(ideg,min_);
		}
/* L1000: */
	    }
/* L1010: */
	}

/*     Maximum freedom number difference */

/* Computing MAX */
	i__2 = *hband, i__3 = max_ - min_;
	*hband = max(i__2,i__3);
/* L1020: */
    }

/*     Semi band width */

    ++(*hband);

} /* bndwth_ */

