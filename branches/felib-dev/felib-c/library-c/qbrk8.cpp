/* qbrk8.f -- translated by f2c (version of 28 August 1991  0:07:02).
   You must link the resulting object file with the libraries:
	-lF77 -lI77 -lm -lc   (in that order)
*/

#include "f2c.h"


/* Subroutine */ int qbrk8_(wght, iwght, abss, iabss, jabss, nqp, itest)
doublereal *wght;
integer *iwght;
doublereal *abss;
integer *iabss, *jabss, *nqp, *itest;
{
    /* Initialized data */

    static char srname[6+1] = "QBRK8 ";

    /* System generated locals */
    integer abss_dim1, abss_offset;

    /* Builtin functions */
    double sqrt();

    /* Local variables */
    static integer i;
    static doublereal w, xy;
    extern integer errmes_();
    static integer ierror;

/* -----------------------------------------------------------------------
 */
/* PURPOSE */
/*      QBRK8 returns weights and abscissae of a 8-point gauss type */
/*      quadrature rule for use in evaluating the integral of a */
/*      3D function over a cube */

/* HISTORY */

/*      Copyright (C) 1996 : CCLRC, Rutherford Appleton Laboratory */
/*                           Chilton, Didcot, Oxfordshire OX11 0QX */

/*      Release 1.1  29 Oct 1979 (CG) */
/*      Commented    15 Oct 1980 (KR) */
/*      Release 4.0   2 Oct 1996 (CG) */

/* ARGUMENTS in */
/*      IWGHT   dimension of vector WGHT (.GE. NQP(=8)) */
/*      IABSS   first dimension of array ABSS (.GE. 3) */
/*      JABSS   second dimension of array ABSS (.GE. NQP(=8)) */
/*      ITEST   error checking option */

/* ARGUMENTS out */
/*      WGHT    vector of dimension IWGHT.  contains weights to */
/*              be used in the 8-point quadrature formula */
/*      ABSS    array of dimension (IABSS, JABSS).  contains */
/*              abscissae of points to be used in quadrature */
/*              formula */
/*      NQP     number of quadrature points to be used (=8) */

/* ROUTINES called */
/*      ERRMES */

/*      SUBROUTINE QBRK8(WGHT,IWGHT,ABSS,IABSS,JABSS,NQP,ITEST) */
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

    w = 1.;
    xy = sqrt(.33333333333333331);
    for (i = 1; i <= 8; ++i) {
	wght[i] = w;
	abss[i * abss_dim1 + 1] = xy;
	abss[i * abss_dim1 + 2] = xy;
	abss[i * abss_dim1 + 3] = xy;
/* L1000: */
    }

    for (i = 1; i <= 4; ++i) {
	abss[i * abss_dim1 + 3] = -abss[abss_dim1 + 1];
/* L1010: */
    }

    for (i = 1; i <= 2; ++i) {
	abss[(i + 1) * abss_dim1 + 1] = -abss[abss_dim1 + 1];
	abss[(i + 5) * abss_dim1 + 1] = -abss[abss_dim1 + 1];
	abss[i * abss_dim1 + 2] = -abss[abss_dim1 + 1];
	abss[(i + 4) * abss_dim1 + 2] = -abss[abss_dim1 + 1];
/* L1020: */
    }

} /* qbrk8_ */

