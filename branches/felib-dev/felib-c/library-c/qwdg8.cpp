/* qwdg8.f -- translated by f2c (version of 28 August 1991  0:07:02).
   You must link the resulting object file with the libraries:
	-lF77 -lI77 -lm -lc   (in that order)
*/

#include "f2c.h"


/* Subroutine */ int qwdg8_(wght, iwght, abss, iabss, jabss, nqp, itest)
doublereal *wght;
integer *iwght;
doublereal *abss;
integer *iabss, *jabss, *nqp, *itest;
{
    /* Initialized data */

    static char srname[6+1] = "QWDG8 ";

    /* System generated locals */
    integer abss_dim1, abss_offset;

    /* Builtin functions */
    double sqrt();

    /* Local variables */
    static integer i;
    static doublereal w;
    extern integer errmes_();
    static integer ierror;

/* -----------------------------------------------------------------------
 */
/* PURPOSE */
/*      QWDG8 returns weights and abscissae for a 8-point quadrature */
/*      rule for evaluating the integral of a 3D function over */
/*      a pentahedral region */

/* HISTORY */

/*      Copyright (C) 1996 : CCLRC, Rutherford Appleton Laboratory */
/*                           Chilton, Didcot, Oxfordshire OX11 0QX */

/*      Release 2.0   6 Oct 1980 (CG) */
/*      Commented    21 Oct 1980 (KR) */
/*      Release 4.0   2 Oct 1996 (CG) */

/* ARGUMENTS in */
/*      IWGHT   dimension of vector WGHT (.GE. NQP(=8)) */
/*      IABSS   first dimension array ABSS (.GE. 3) */
/*      JABSS   second dimension of array ABSS (.GE. NQP(=8)) */
/*      ITEST   error checking option */

/* ARGUMENTS out */
/*      WGHT    vector of length IWGHT.  contains weights of */
/*              quadrature rule */
/*      ABSS    array of dimension (IABSS,JABSS).  contains */
/*              abscissae of quadrature points */

/* ROUTINES called */
/*      ERRMES */

/*     SUBROUTINE QWDG8(WGHT,IWGHT,ABSS,IABSS,JABSS,NQP,ITEST) */
/* ***********************************************************************
 */

    /* Parameter adjustments */
    abss_dim1 = *iabss;
    abss_offset = abss_dim1 + 1;
    abss -= abss_offset;
    --wght;

    /* Function Body */

    *nqp = 8;

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

    w = sqrt(3.) * 3. / 2.;
    wght[1] = w * .375;
    wght[2] = w * .375;
    for (i = 3; i <= 8; ++i) {
	wght[i] = w * .041666666666666664;
/* L1000: */
    }

    abss[abss_dim1 + 1] = 0.;
    abss[abss_dim1 + 2] = 0.;
    abss[abss_dim1 + 3] = 1.;
    abss[(abss_dim1 << 1) + 1] = 0.;
    abss[(abss_dim1 << 1) + 2] = 0.;
    abss[(abss_dim1 << 1) + 3] = -1.;
    abss[abss_dim1 * 3 + 1] = 1.;
    abss[abss_dim1 * 3 + 2] = 0.;
    abss[abss_dim1 * 3 + 3] = -1.;
    abss[(abss_dim1 << 2) + 1] = -.5;
    abss[(abss_dim1 << 2) + 2] = sqrt(3.) / 2.;
    abss[(abss_dim1 << 2) + 3] = -1.;
    abss[abss_dim1 * 5 + 1] = -.5;
    abss[abss_dim1 * 5 + 2] = -abss[(abss_dim1 << 2) + 2];
    abss[abss_dim1 * 5 + 3] = -1.;
    abss[abss_dim1 * 6 + 1] = 1.;
    abss[abss_dim1 * 6 + 2] = 0.;
    abss[abss_dim1 * 6 + 3] = 1.;
    abss[abss_dim1 * 7 + 1] = -.5;
    abss[abss_dim1 * 7 + 2] = abss[(abss_dim1 << 2) + 2];
    abss[abss_dim1 * 7 + 3] = 1.;
    abss[(abss_dim1 << 3) + 1] = -.5;
    abss[(abss_dim1 << 3) + 2] = abss[abss_dim1 * 5 + 2];
    abss[(abss_dim1 << 3) + 3] = 1.;

} /* qwdg8_ */

