/* qlin2.f -- translated by f2c (version of 28 August 1991  0:07:02).
   You must link the resulting object file with the libraries:
	-lF77 -lI77 -lm -lc   (in that order)
*/

#include "f2c.h"


/* Subroutine */ int qlin2_(wght, iwght, abss, iabss, nqp, itest)
doublereal *wght;
integer *iwght;
doublereal *abss;
integer *iabss, *nqp, *itest;
{
    /* Initialized data */

    static char srname[6+1] = "QLIN2 ";

    /* Builtin functions */
    double sqrt();

    /* Local variables */
    extern integer errmes_();
    static integer ierror;

/* -----------------------------------------------------------------------
 */
/* PURPOSE */
/*      QLIN2 returns the weights and abscissae of a 2-point gauss- */
/*      legendre quadrature formula for use in evaluating a 1D */
/*      integral over a finite range */

/* HISTORY */

/*      Copyright (C) 1996 : CCLRC, Rutherford Appleton Laboratory */
/*                           Chilton, Didcot, Oxfordshire OX11 0QX */

/*      Release 1.1  29 Oct 1979 (CG) */
/*      Commented    15 Oct 1980 (KR) */
/*      Release 4.0   2 Oct 1996 (CG) */

/* ARGUMENTS in */
/*      IWGHT   dimension of vector WGHT (.GE. NQP(=2)) */
/*      IABSS   dimension of vector ABSS (.GE. NQP(=2)) */
/*      ITEST   error checking option */

/* ARGUMENTS out */
/*      WGHT    vector of dimension IWGHT.  contains weights to */
/*              be used in the 2-point quadrature formula */
/*      ABSS    vector of dimension IABSS.  contains abscissae */
/*              of points to be used in 2-point quadrature */
/*              formula */
/*      NQP     number of quadrature points (=2) */

/* ROUTINES called */
/*      ERRMES */

/*      SUBROUTINE QLIN2(WGHT,IWGHT,ABSS,IABSS,NQP,ITEST) */
/* ***********************************************************************
 */

    /* Parameter adjustments */
    --abss;
    --wght;

    /* Function Body */

    *nqp = 2;

/*     Parameter checking */

    if (*itest != -1) {
	ierror = 0;
	if (*iwght < *nqp) {
	    ierror = 1;
	}
	if (*iabss < *nqp) {
	    ierror = 2;
	}
	*itest = errmes_(itest, &ierror, srname, 6L);
	if (*itest != 0) {
	    return 0;
	}
    }

/*     Main body */

    abss[1] = 1. / sqrt(3.);
    abss[2] = -abss[1];
    wght[1] = 1.;
    wght[2] = 1.;

} /* qlin2_ */

