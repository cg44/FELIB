/* bqqua.f -- translated by f2c (version of 28 August 1991  0:07:02).
   You must link the resulting object file with the libraries:
	-lF77 -lI77 -lm -lc   (in that order)
*/

#include "f2c.h"


/* Subroutine */ int bqqua_(abss, iabss, jabss, work, iwork, nqp, sidnum, 
	coef, itest)
doublereal *abss;
integer *iabss, *jabss;
doublereal *work;
integer *iwork, *nqp, *sidnum;
doublereal *coef;
integer *itest;
{
    /* Initialized data */

    static char srname[6+1] = "BQQUA ";

    /* System generated locals */
    integer abss_dim1, abss_offset, i__1;

    /* Local variables */
    static integer i, j, k;
    extern integer errmes_();
    static integer ierror;
    static doublereal val;

/* -----------------------------------------------------------------------
 */
/* PURPOSE */
/*      BQQUA forms an equivalent one-dimensional quadrature rule of */
/*      a given two-dimensional rule for integration along the */
/*      specified side of a rectangular element. */

/* HISTORY */

/*      Copyright (C) 1996 : CCLRC, Rutherford Appleton Laboratory */
/*                           Chilton, Didcot, Oxfordshire OX11 0QX */

/*      Release 2.0   1 Feb 1981 (CG) */
/*      Commented    24 Jun 1981 (CG) */
/*      Release 4.0   2 Oct 1996 (CG) */

/* ARGUMENTS in */
/*      ABSS    array holding the abscissae of the two- */
/*              dimensional quadrature rule to be used */
/*      IABSS   first dimension of array ABSS */
/*      JABSS   second dimension of array ABSS */
/*      IWORK   dimension of array WORK */
/*      NQP     number of quadrature points in the rule */
/*      SIDNUM  the side number for which the one- */
/*              dimensional rule is required */
/*      ITEST   error checking option */

/* ARGUMENTS out */
/*      WORK    array containing the abscissae of the */
/*              equivalent one-dimensional rule */
/*      COEF    a multiplier of the rule */

/* ROUTINES called */
/*      ERRMES */

/*      SUBROUTINE BQQUA(ABSS,IABSS,JABSS,WORK,IWORK,NQP,SIDNUM,COEF, */
/*     *                 ITEST) */
/* ***********************************************************************
 */



    /* Parameter adjustments */
    --work;
    abss_dim1 = *iabss;
    abss_offset = abss_dim1 + 1;
    abss -= abss_offset;

    /* Function Body */

/*     Parameter checking */

    if (*itest != -1) {
	ierror = 0;
	if (*iabss < 2 || *jabss < *nqp) {
	    ierror = 4;
	}
	if (*iwork < *nqp) {
	    ierror = 3;
	}
	if (*sidnum <= 0 || *sidnum > 4) {
	    ierror = 2;
	}
	if (*nqp <= 0) {
	    ierror = 1;
	}
	*itest = errmes_(itest, &ierror, srname, 6L);
	if (*itest != 0) {
	    return 0;
	}
    }

/*     Main loops */

    *coef = 1.;
    switch ((int)*sidnum) {
	case 1:  goto L1000;
	case 2:  goto L1010;
	case 3:  goto L1020;
	case 4:  goto L1030;
    }

/*     Side number 1 : xi = -1 */

L1000:
    i = 1;
    j = 2;
    val = -1.;
    goto L1040;

/*     Side number 2 : xi = +1 */

L1010:
    i = 2;
    j = 1;
    val = 1.;
    goto L1040;

/*     Side number 3 : eta = -1 */

L1020:
    i = 1;
    j = 2;
    val = 1.;
    goto L1040;

/*     Side number 4 : eta = +1 */

L1030:
    i = 2;
    j = 1;
    val = -1.;

L1040:
    i__1 = *nqp;
    for (k = 1; k <= i__1; ++k) {
	abss[i + k * abss_dim1] = val;
	abss[j + k * abss_dim1] = work[k];
/* L1050: */
    }

} /* bqqua_ */

