/* qlin3.f -- translated by f2c (version of 28 August 1991  0:07:02).
   You must link the resulting object file with the libraries:
	-lF77 -lI77 -lm -lc   (in that order)
*/

#include "f2c.h"


/* Subroutine */ int qlin3_(wght, iwght, abss, iabss, nqp, itest)
doublereal *wght;
integer *iwght;
doublereal *abss;
integer *iabss, *nqp, *itest;
{
    /* Initialized data */

    static char srname[6+1] = "QLIN3 ";

    /* Builtin functions */
    double sqrt();

    /* Local variables */
    extern integer errmes_();
    static integer ierror;

/* -----------------------------------------------------------------------
 */
/* PURPOSE */
/*      QLIN3 returns weights and abscissae of a 3-point gauss- */
/*      legendre quadrature formula for use in evaluating a 1D */
/*      integral over a finite range */

/* HISTORY */

/*      Copyright (C) 1996 : CCLRC, Rutherford Appleton Laboratory */
/*                           Chilton, Didcot, Oxfordshire OX11 0QX */

/*      Release 1.1  29 Oct 1979 (CG) */
/*      Commented    15 Oct 1980 (KR) */
/*      Release 4.0   2 Oct 1996 (CG) */

/* ARGUMENTS in */
/*      IWGHT   dimension of vector WGHT (.GE. NQP(=3)) */
/*      IABSS   dimension of vector ABSS (.GE. NQP(=3)) */
/*      ITEST   error checking option */

/* ARGUMENTS out */
/*      WGHT    vector of dimension IWGHT. Contains weights for */
/*              3-point formula */
/*      ABSS    vector of dimension IABSS.  contains abscissae */
/*              of points for use in 3-point formula */
/*      NQP     number of quadrature points (=3) */

/* ROUTINES called */
/*      ERRMES */

/*      SUBROUTINE QLIN3(WGHT,IWGHT,ABSS,IABSS,NQP,ITEST) */
/* ***********************************************************************
 */

    /* Parameter adjustments */
    --abss;
    --wght;

    /* Function Body */

    *nqp = 3;

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

    abss[1] = sqrt(15.) * .2;
    abss[2] = 0.;
    abss[3] = -abss[1];
    wght[1] = .55555555555555558;
    wght[3] = wght[1];
    wght[2] = .88888888888888884;

} /* qlin3_ */

