/* sideno.f -- translated by f2c (version of 28 August 1991  0:07:02).
   You must link the resulting object file with the libraries:
	-lF77 -lI77 -lm -lc   (in that order)
*/

#include "f2c.h"


/* Subroutine */ int sideno_(totels, eltop, ieltop, jeltop, m, bdcnd, ibdcnd, 
	jbdcnd, numsid, blist, iblist, jblist, itest)
integer *totels, *eltop, *ieltop, *jeltop, *m, *bdcnd, *ibdcnd, *jbdcnd, *
	numsid, *blist, *iblist, *jblist, *itest;
{
    /* Initialized data */

    static char srname[6+1] = "SIDENO";

    /* System generated locals */
    integer bdcnd_dim1, bdcnd_offset, blist_dim1, blist_offset, eltop_dim1, 
	    eltop_offset, i__1, i__2, i__3, i__4;

    /* Local variables */
    static integer nele, bnum1, i, j, k, l, n, nodel, jtest, n1, nside1, nb, 
	    ik;
    extern integer errmes_();
    static integer ierror;

/* -----------------------------------------------------------------------
 */
/* PURPOSE */
/*      SIDENO the routine generates a list of element and side */
/*      numbers from a list of boundary nodes. the node list */
/*      must be a continuous sequence of node in the same */
/*      sense as the local node ordering. */

/* HISTORY */

/*      Copyright (C) 1996 : CCLRC, Rutherford Appleton Laboratory */
/*                           Chilton, Didcot, Oxfordshire OX11 0QX */

/*      Release 2.0  Jan 1982 (CG) */
/*      Release 4.0   2 Oct 1996 (CG) */

/* ARGUMENTS in */
/*      TOTELS  the total number of elements in the mesh */
/*      ELTOP   the element topology array containing element */
/*              type number, number of nodes on element and */
/*              topology list */
/*      IELTOP  first dimension of ELTOP */
/*      JELTOP  second dimension of ELTOP */
/*      M       boundary condition list number */
/*      BDCND   boundary condition array containing boundary */
/*              condition type, number of nodes on element side */
/*              number of nodes in list and node list. */
/*      IBDCND  first dimension of BDCND array */
/*      JBDCND  second dimension of BDCND array */
/*      IBLIST  first dimension BLIST */
/*      JBLIST  second dimension of BLIST */
/*      ITEST   error checking option */

/* ARGUMENTS out */
/*      NUMSID  number of element sides generated from the */
/*              current boundary node list */
/*      BLIST   array containing element number and side number */

/* ROUTINES called */
/*      ERRMES */

/*     SUBROUTINE SIDENO(TOTELS,ELTOP,IELTOP,JELTOP,M,BDCND,IBDCND, */
/*    *                  JBDCND,NUMSID,BLIST,IBLIST,JBLIST,ITEST) */
/* ***********************************************************************
 */

    /* Parameter adjustments */
    blist_dim1 = *iblist;
    blist_offset = blist_dim1 + 1;
    blist -= blist_offset;
    bdcnd_dim1 = *ibdcnd;
    bdcnd_offset = bdcnd_dim1 + 1;
    bdcnd -= bdcnd_offset;
    eltop_dim1 = *ieltop;
    eltop_offset = eltop_dim1 + 1;
    eltop -= eltop_offset;

    /* Function Body */

/*     Parameter checking */

    jtest = *itest;
    if (jtest != -1) {
	ierror = 0;
	if (*jblist < 2) {
	    ierror = 4;
	}
	if (*ibdcnd < *m) {
	    ierror = 3;
	}
	if (*ieltop < *totels) {
	    ierror = 2;
	}
	if (*totels <= 0 || *m <= 0) {
	    ierror = 1;
	}
	*itest = errmes_(&jtest, &ierror, srname, 6L);
	if (*itest != 0) {
	    return 0;
	}
    }

/*     Main body */

    nside1 = bdcnd[*m + bdcnd_dim1 * 3] - 1;
    bnum1 = bdcnd[*m + (bdcnd_dim1 << 1)] - nside1;
    n = 0;
    i__1 = bnum1;
    i__2 = nside1;
    for (i = 1; i__2 < 0 ? i >= i__1 : i <= i__1; i += i__2) {

/*     Range checking on I */

	if (jtest != -1) {
	    ierror = 0;
	    if (i > *jbdcnd) {
		ierror = 6;
	    }
	    *itest = errmes_(&jtest, &ierror, srname, 6L);
	    if (*itest != 0) {
		return 0;
	    }
	}

	nb = bdcnd[*m + (i + 3) * bdcnd_dim1];
	i__3 = *totels;
	for (j = 1; j <= i__3; ++j) {
	    nodel = eltop[j + (eltop_dim1 << 1)];
	    i__4 = nodel;
	    for (k = 1; k <= i__4; ++k) {

/*     Range checking on K */

		if (jtest != -1) {
		    ierror = 0;
		    if (k + 2 > *jeltop) {
			ierror = 5;
		    }
		    *itest = errmes_(&jtest, &ierror, srname, 6L);
		    if (*itest != 0) {
			return 0;
		    }
		}

		if (nb == eltop[j + (k + 2) * eltop_dim1]) {
		    goto L1010;
		}
/* L1000: */
	    }
	    goto L1030;

L1010:
	    nele = j;
	    if (nodel == 10) {
		--nodel;
	    }
	    n1 = k;
	    i__4 = nside1;
	    for (k = 1; k <= i__4; ++k) {
		l = (n1 + k - 1) % nodel + 3;
		ik = i + k + 3;

/*     Range checking on L */

		if (jtest != -1) {
		    ierror = 0;
		    if (l > *jeltop) {
			ierror = 5;
		    }
		    *itest = errmes_(&jtest, &ierror, srname, 6L);
		    if (*itest != 0) {
			return 0;
		    }
		}

		if (bdcnd[*m + ik * bdcnd_dim1] != eltop[j + l * eltop_dim1]) 
			{
		    goto L1030;
		}
/* L1020: */
	    }
	    goto L1040;
L1030:
	    ;
	}
	goto L1050;

L1040:
	++n;

/*     Range checking on N */

	if (jtest != -1) {
	    ierror = 0;
	    if (n > *iblist) {
		ierror = 7;
	    }
	    *itest = errmes_(&jtest, &ierror, srname, 6L);
	    if (*itest != 0) {
		return 0;
	    }
	}

	blist[n + blist_dim1] = nele;
	blist[n + (blist_dim1 << 1)] = (n1 - 1) / nside1 + 1;
L1050:
	;
    }

    *numsid = n;

} /* sideno_ */

