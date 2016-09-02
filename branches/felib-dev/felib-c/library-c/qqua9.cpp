/* qqua9.f -- translated by f2c (version of 28 August 1991  0:07:02).
   You must link the resulting object file with the libraries:
	-lF77 -lI77 -lm -lc   (in that order)
*/

#include "f2c.h"


/* Subroutine */ int qqua9_(wght, iwght, abss, iabss, jabss, nqp, itest)
doublereal *wght;
integer *iwght;
doublereal *abss;
integer *iabss, *jabss, *nqp, *itest;
{
    /* Initialized data */

    static char srname[6+1] = "QQUA9 ";

    /* System generated locals */
    integer abss_dim1, abss_offset;

    /* Builtin functions */
    double sqrt();

    /* Local variables */
    static doublereal area, xy;
    extern integer errmes_();
    static integer ierror;

/* -----------------------------------------------------------------------
 */
/* PURPOSE */
/*      QQUA9 returns weights and abscissae of a 9-point gaussian */
/*      product quadrature rule for evaluating the integral of a */
/*      2D function over a rectangular region */

/* HISTORY */

/*      Copyright (C) 1996 : CCLRC, Rutherford Appleton Laboratory */
/*                           Chilton, Didcot, Oxfordshire OX11 0QX */

/*      Release 1.1  29 Oct 1979 (CG) */
/*      Commented    16 Oct 1980 (KR) */
/*      Release 4.0   2 Oct 1996 (CG) */

/* ARGUMENTS in */
/*      IWGHT   dimension of vector WGHT (.GE. NQP) */
/*      IABSS   first dimension of array ABSS (.GE. 2) */
/*      JABSS   second dimension of array ABSS (.GE. NQP) */
/*      ITEST   error checking option */

/* ARGUMENTS out */
/*      WGHT    vector of length IWGHT.  contains weights to be */
/*              used in quadrature formula */
/*      ABSS    contains abscissae of points to be used in */
/*              formula */
/*      NQP     number of quadrature points (=9) */

/* ROUTINES called */
/*      ERRMES */

/*      SUBROUTINE QQUA9(WGHT,IWGHT,ABSS,IABSS,JABSS,NQP,ITEST) */
/* ***********************************************************************
 */

    /* Parameter adjustments */
    abss_dim1 = *iabss;
    abss_offset = abss_dim1 + 1;
    abss -= abss_offset;
    --wght;

    /* Function Body */

    *nqp = 9;

/*     Parameter checking */

    if (*itest != -1) {
	ierror = 0;
	if (*iwght < *nqp) {
	    ierror = 1;
	}
	if (*iabss < 2 || *jabss < *nqp) {
	    ierror = 2;
	}
	*itest = errmes_(itest, &ierror, srname, 6L);
	if (*itest != 0) {
	    return 0;
	}
    }

/*     Main body */

    area = 4.;
    wght[1] = area * 25. / 324.;
    wght[3] = wght[1];
    wght[5] = wght[1];
    wght[7] = wght[1];
    wght[2] = area * 10. / 81.;
    wght[4] = wght[2];
    wght[6] = wght[2];
    wght[8] = wght[2];
    wght[9] = area * 16. / 81.;

    xy = sqrt(.59999999999999998);
    abss[abss_dim1 + 1] = -xy;
    abss[abss_dim1 + 2] = -xy;
    abss[(abss_dim1 << 1) + 1] = -xy;
    abss[(abss_dim1 << 1) + 2] = 0.;
    abss[abss_dim1 * 3 + 1] = -xy;
    abss[abss_dim1 * 3 + 2] = xy;
    abss[(abss_dim1 << 2) + 1] = 0.;
    abss[(abss_dim1 << 2) + 2] = xy;
    abss[abss_dim1 * 5 + 1] = xy;
    abss[abss_dim1 * 5 + 2] = xy;
    abss[abss_dim1 * 6 + 1] = xy;
    abss[abss_dim1 * 6 + 2] = 0.;
    abss[abss_dim1 * 7 + 1] = xy;
    abss[abss_dim1 * 7 + 2] = -xy;
    abss[(abss_dim1 << 3) + 1] = 0.;
    abss[(abss_dim1 << 3) + 2] = -xy;
    abss[abss_dim1 * 9 + 1] = 0.;
    abss[abss_dim1 * 9 + 2] = 0.;

} /* qqua9_ */

