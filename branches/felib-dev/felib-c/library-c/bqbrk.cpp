/* bqbrk.f -- translated by f2c (version of 28 August 1991  0:07:02).
   You must link the resulting object file with the libraries:
	-lF77 -lI77 -lm -lc   (in that order)
*/

#include "f2c.h"


/* Subroutine */ int bqbrk_(abss, iabss, jabss, work, iwork, jwork, nqp, 
	facnum, coef, itest)
doublereal *abss;
integer *iabss, *jabss;
doublereal *work;
integer *iwork, *jwork, *nqp, *facnum;
doublereal *coef;
integer *itest;
{
    /* Initialized data */

    static char srname[6+1] = "BQBRK ";

    /* System generated locals */
    integer abss_dim1, abss_offset, work_dim1, work_offset, i__1;

    /* Local variables */
    static integer i, j, k, l;
    extern integer errmes_();
    static integer ierror;
    static doublereal val;

/* -----------------------------------------------------------------------
 */
/* PURPOSE */
/*      BQBRK forms an equivalent two-dimensional quadrature rule of */
/*      a given three-dimensional rule for integration over the */
/*      specified face of a brick element. */

/* HISTORY */

/*      Copyright (C) 1996 : CCLRC, Rutherford Appleton Laboratory */
/*                           Chilton, Didcot, Oxfordshire OX11 0QX */

/*      Release 2.0   1 Feb 1981 (CG) */
/*      Commented    24 Jun 1981 (CG) */
/*      Release 4.0   2 Oct 1996 (CG) */

/* ARGUMENTS in */
/*      ABSS    array holding the abscissae of the three- */
/*              dimensional quadrature rule to be used */
/*      IABSS   first dimension of array ABSS */
/*      JABSS   second dimension of array ABSS */
/*      IWORK   first dimension of array WORK */
/*      JWORK   second dimension of array WORK */
/*      NQP     number of quadrature points in the rule */
/*      FACNUM  the face number for which the two- */
/*              dimensional rule is required */
/*      ITEST   error checking option */

/* ARGUMENTS out */
/*      WORK    array containing the abscissae of the */
/*              equivalent one-dimensional rule */
/*      COEF    a multiplier of the rule */

/* ROUTINES called */
/*      ERRMES */

/*      SUBROUTINE BQBRK(ABSS,IABSS,JABSS,WORK,IWORK,JWORK,NQP,FACNUM, */
/*     *                 COEF,ITEST) */
/* ***********************************************************************
 */



    /* Parameter adjustments */
    work_dim1 = *iwork;
    work_offset = work_dim1 + 1;
    work -= work_offset;
    abss_dim1 = *iabss;
    abss_offset = abss_dim1 + 1;
    abss -= abss_offset;

    /* Function Body */

/*     Parameter checking */

    if (*itest != -1) {
	ierror = 0;
	if (*nqp == 0) {
	    ierror = 1;
	}
	if (*facnum <= 0 || *facnum > 6) {
	    ierror = 2;
	}
	if (*iwork < 2 || *jwork < *nqp) {
	    ierror = 3;
	}
	if (*iabss < 3 || *jabss < *nqp) {
	    ierror = 4;
	}
	*itest = errmes_(itest, &ierror, srname, 6L);
	if (*itest != 0) {
	    return 0;
	}
    }

/*     Main loops */

    switch ((int)*facnum) {
	case 1:  goto L1000;
	case 2:  goto L1010;
	case 3:  goto L1020;
	case 4:  goto L1030;
	case 5:  goto L1040;
	case 6:  goto L1050;
    }
    goto L1020;

/*     Face number 1 : zeta = -1 */

L1000:
    i = 3;
    j = 1;
    k = 2;
    val = -1.;
    *coef = 1.;
    goto L1060;

/*     Face number 2 : zeta = +1 */

L1010:
    i = 3;
    j = 1;
    k = 2;
    val = 1.;
    *coef = 1.;
    goto L1060;

/*     Face number 3 : xi = -1 */

L1020:
    i = 1;
    j = 2;
    k = 3;
    val = -1.;
    *coef = 1.;
    goto L1060;

/*     Face number 4 : eta = +1 */

L1030:
    i = 2;
    j = 1;
    k = 3;
    val = 1.;
    *coef = 1.;
    goto L1060;

/*     Face number 5 : xi = +1 */

L1040:
    i = 1;
    j = 2;
    k = 3;
    val = 1.;
    *coef = 1.;
    goto L1060;

/*     Face number 6 : eta = -1 */

L1050:
    i = 2;
    j = 1;
    k = 3;
    val = -1.;
    *coef = 1.;

L1060:
    i__1 = *nqp;
    for (l = 1; l <= i__1; ++l) {
	abss[i + l * abss_dim1] = val;
	abss[j + l * abss_dim1] = work[l * work_dim1 + 1];
	abss[k + l * abss_dim1] = work[l * work_dim1 + 2];
/* L1070: */
    }

} /* bqbrk_ */

