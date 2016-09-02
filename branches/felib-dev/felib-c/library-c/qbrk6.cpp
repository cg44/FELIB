/* qbrk6.f -- translated by f2c (version of 28 August 1991  0:07:02).
   You must link the resulting object file with the libraries:
	-lF77 -lI77 -lm -lc   (in that order)
*/

#include "f2c.h"


/* Subroutine */ int qbrk6_(wght, iwght, abss, iabss, jabss, nqp, itest)
doublereal *wght;
integer *iwght;
doublereal *abss;
integer *iabss, *jabss, *nqp, *itest;
{
    /* Initialized data */

    static char srname[6+1] = "QBRK6 ";

    /* System generated locals */
    integer abss_dim1, abss_offset;

    /* Local variables */
    static integer i;
    static doublereal w;
    extern integer errmes_();
    static integer ierror;

/* -----------------------------------------------------------------------
 */
/* PURPOSE */
/*      QBRK6 returns weights and abscissae of a six-point gauss type */
/*      quadrature rule for use in evaluating the integral of a */
/*      3D function over a cube */

/* HISTORY */

/*      Copyright (C) 1996 : CCLRC, Rutherford Appleton Laboratory */
/*                           Chilton, Didcot, Oxfordshire OX11 0QX */

/*      Release 1.1  29 Oct 1979 (CG) */
/*      Commented    15 Oct 1980 (KR) */
/*      Release 4.0   2 Oct 1996 (CG) */

/* ARGUMENTS in */
/*      IWGHT   dimension of vector WGHT (.GE. NQP(=6)) */
/*      IABSS   first dimension of array ABSS (.GE. 3) */
/*      JABSS   second dimension of array ABSS (.GE. NQP(=6)) */
/*      ITEST   error checking option */

/* ARGUMENTS out */
/*      WGHT    vector of dimension IWGHT.  contains weights to */
/*              be used in the 6-point quadrature formula */
/*      ABSS    array of dimension (IABSS, JABSS).  contains */
/*              abscissae of points to be used in quadrature */
/*              formula */
/*      NQP     number of quadrature points to be used (=6) */

/* ROUTINES called */
/*      ERRMES */

/*      SUBROUTINE QBRK6(WGHT,IWGHT,ABSS,IABSS,JABSS,NQP,ITEST) */
/* ***********************************************************************
 */

    /* Parameter adjustments */
    abss_dim1 = *iabss;
    abss_offset = abss_dim1 + 1;
    abss -= abss_offset;
    --wght;

    /* Function Body */

    *nqp = 6;

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

    w = 1.3333333333333333;

    for (i = 1; i <= 6; ++i) {
	wght[i] = w;
	abss[i * abss_dim1 + 1] = 0.;
	abss[i * abss_dim1 + 2] = 0.;
	abss[i * abss_dim1 + 3] = 0.;
/* L1000: */
    }

    abss[abss_dim1 * 5 + 1] = -1.;
    abss[abss_dim1 * 6 + 1] = 1.;
    abss[abss_dim1 * 3 + 2] = -1.;
    abss[(abss_dim1 << 2) + 2] = 1.;
    abss[abss_dim1 + 3] = -1.;
    abss[(abss_dim1 << 1) + 3] = 1.;

} /* qbrk6_ */

