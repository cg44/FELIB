/* qtet4.f -- translated by f2c (version of 28 August 1991  0:07:02).
   You must link the resulting object file with the libraries:
	-lF77 -lI77 -lm -lc   (in that order)
*/

#include "f2c.h"


/* Subroutine */ int qtet4_(wght, iwght, abss, iabss, jabss, nqp, itest)
doublereal *wght;
integer *iwght;
doublereal *abss;
integer *iabss, *jabss, *nqp, *itest;
{
    /* Initialized data */

    static char srname[6+1] = "QTET4 ";

    /* System generated locals */
    integer abss_dim1, abss_offset;

    /* Builtin functions */
    double sqrt();

    /* Local variables */
    static doublereal gama, beta;
    static integer i;
    static doublereal alpha;
    extern integer errmes_();
    static integer ierror;

/* -----------------------------------------------------------------------
 */
/* PURPOSE */
/*      QTET4 returns weights and abscissae of a 4-point rule for */
/*      evaluating the integral of a 3D function over a */
/*      tetrahedral region */

/* HISTORY */

/*      Copyright (C) 1996 : CCLRC, Rutherford Appleton Laboratory */
/*                           Chilton, Didcot, Oxfordshire OX11 0QX */

/*      Release 2.0   6 Oct 1980 (CG) */
/*      Commented    16 Oct 1980 (KR) */
/*      Release 4.0   2 Oct 1996 (CG) */

/* ARGUMENTS in */
/*      IWGHT   dimension of vector WGHT (.GE. NQP) */
/*      IABSS   first dimension of array ABSS (.GE. 3) */
/*      JABSS   second dimension of array ABSS (.GE. NQP) */
/*      ITEST   error checking option */

/* ARGUMENTS out */
/*      WGHT    vector of dimension IWGHT.  contains weights for */
/*              4-point quadrature rule */
/*      ABSS    array of dimension (IABSS,JABSS).  contains */
/*              abscissae of points to be used in quadrature */
/*              rule */
/*      NQP     number of quadrature points (=4) */

/* ROUTINES called */
/*      ERRMES */

/*      SUBROUTINE QTET4(WGHT,IWGHT,ABSS,IABSS,JABSS,NQP,ITEST) */
/* ***********************************************************************
 */

    /* Parameter adjustments */
    abss_dim1 = *iabss;
    abss_offset = abss_dim1 + 1;
    abss -= abss_offset;
    --wght;

    /* Function Body */

    *nqp = 4;

/*     Parameter checking */

    if (*itest != -1) {
	ierror = 0;
	if (*iwght < *nqp) {
	    ierror = 1;
	}
	if (*iabss < 3 || *jabss < *nqp) {
	    ierror = 2;
	}
	*itest = errmes_(itest, &ierror, srname, 6L);
	if (*itest != 0) {
	    return 0;
	}
    }

/*     Main body */

    alpha = .223807;
    beta = .387298;
    gama = .158114;

    for (i = 1; i <= 4; ++i) {
	wght[i] = sqrt(6.) * .0625;
/* L1000: */
    }

    abss[abss_dim1 + 1] = alpha * 2.;
    abss[abss_dim1 + 2] = 0.;
    abss[abss_dim1 + 3] = -gama;
    abss[(abss_dim1 << 1) + 1] = -alpha;
    abss[(abss_dim1 << 1) + 2] = -beta;
    abss[(abss_dim1 << 1) + 3] = -gama;
    abss[abss_dim1 * 3 + 1] = -alpha;
    abss[abss_dim1 * 3 + 2] = beta;
    abss[abss_dim1 * 3 + 3] = -gama;
    abss[(abss_dim1 << 2) + 1] = 0.;
    abss[(abss_dim1 << 2) + 2] = 0.;
    abss[(abss_dim1 << 2) + 3] = gama * 3.;

} /* qtet4_ */

