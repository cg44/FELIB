/* bqtri.f -- translated by f2c (version of 28 August 1991  0:07:02).
   You must link the resulting object file with the libraries:
	-lF77 -lI77 -lm -lc   (in that order)
*/

#include "f2c.h"


/* Subroutine */ int bqtri_(abss, iabss, jabss, work, iwork, nqp, sidnum, 
	coef, itest)
doublereal *abss;
integer *iabss, *jabss;
doublereal *work;
integer *iwork, *nqp, *sidnum;
doublereal *coef;
integer *itest;
{
    /* Initialized data */

    static char srname[6+1] = "BQTRI ";

    /* System generated locals */
    integer abss_dim1, abss_offset, i__1;

    /* Builtin functions */
    double sqrt();

    /* Local variables */
    static integer i;
    static doublereal l1, l2, l3;
    extern integer errmes_();
    static integer ierror;

/* -----------------------------------------------------------------------
 */
/* PURPOSE */
/*      BQTRI forms an equivalent one-dimensional quadrature rule of */
/*      a given two-dimensional rule for integration along the */
/*      specified side of a triangular element. */

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

/*      SUBROUTINE BQTRI(ABSS,IABSS,JABSS,WORK,IWORK,NQP,SIDNUM,COEF, */
/*     *                 ITEST) */
/* ***********************************************************************
 */



    /* Parameter adjustments */
    --work;
    abss_dim1 = *iabss;
    abss_offset = abss_dim1 + 1;
    abss -= abss_offset;

    /* Function Body */

/*     Statement functions */


/*     Parameter checking */

    if (*itest != -1) {
	ierror = 0;
	if (*iabss <= 2 || *jabss <= *nqp) {
	    ierror = 4;
	}
	if (*iwork <= *nqp) {
	    ierror = 3;
	}
	if (*sidnum <= 0 || *sidnum > 3) {
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

    *coef = .5;
    switch ((int)*sidnum) {
	case 1:  goto L1000;
	case 2:  goto L1020;
	case 3:  goto L1040;
    }

/*     Side number 1 : L3 = 0 */

L1000:
    l3 = 0.;
    i__1 = *nqp;
    for (i = 1; i <= i__1; ++i) {
	l1 = (work[i] + 1.) * .5;
	l2 = 1. - l1;
	abss[i * abss_dim1 + 1] = 1. - (l2 + l3) * 1.5;
	abss[i * abss_dim1 + 2] = sqrt(3.) / 2. * (l3 - l2);
/* L1010: */
    }
    return 0;

/*     Side number 2 : L1 = 0 */

L1020:
    l1 = 0.;
    i__1 = *nqp;
    for (i = 1; i <= i__1; ++i) {
	l2 = (work[i] + 1.) * .5;
	l3 = 1. - l2;
	abss[i * abss_dim1 + 1] = 1. - (l2 + l3) * 1.5;
	abss[i * abss_dim1 + 2] = sqrt(3.) / 2. * (l3 - l2);
/* L1030: */
    }
    return 0;

/*     Side number 3 : L2 = 0 */

L1040:
    l2 = 0.;
    i__1 = *nqp;
    for (i = 1; i <= i__1; ++i) {
	l3 = (work[i] + 1.) * .5;
	l1 = 1. - l3;
	abss[i * abss_dim1 + 1] = 1. - (l2 + l3) * 1.5;
	abss[i * abss_dim1 + 2] = sqrt(3.) / 2. * (l3 - l2);
/* L1050: */
    }
} /* bqtri_ */

