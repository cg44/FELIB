/* qtri7.f -- translated by f2c (version of 28 August 1991  0:07:02).
   You must link the resulting object file with the libraries:
	-lF77 -lI77 -lm -lc   (in that order)
*/

#include "f2c.h"


/* Subroutine */ int qtri7_(wght, iwght, abss, iabss, jabss, nqp, itest)
doublereal *wght;
integer *iwght;
doublereal *abss;
integer *iabss, *jabss, *nqp, *itest;
{
    /* Initialized data */

    static char srname[6+1] = "QTRI7 ";

    /* System generated locals */
    integer abss_dim1, abss_offset;

    /* Builtin functions */
    double sqrt();

    /* Local variables */
    static doublereal area;
    extern integer errmes_();
    static integer ierror;

/* -----------------------------------------------------------------------
 */
/* PURPOSE */
/*      QTRI7 returns weights and abscissae of a 7-point gauss-type */
/*      quadrature rule for evaluating the integral of a 2D */
/*      function over a triangular region */

/* HISTORY */

/*      Copyright (C) 1996 : CCLRC, Rutherford Appleton Laboratory */
/*                           Chilton, Didcot, Oxfordshire OX11 0QX */

/*      Release 1.1  29 Oct 1979 (CG) */
/*      Commented    21 Oct 1980 (KR) */
/*      Release 4.0   2 Oct 1996 (CG) */

/* ARGUMENTS in */
/*      IWGHT   dimension of vector WGHT (.GE. NQP) */
/*      IABSS   first dimension of array ABSS (.GE. 2) */
/*      JABSS   second dimension of array ABSS (.GE. NQP) */
/*      ITEST   error checking option */

/* ARGUMENTS out */
/*      WGHT    vector of length IWGHT, containing NQP weights */
/*              of quadrature formula */
/*      ABSS    array of dimension (IABSS,JABSS).  contains */
/*              abscissae of points to be used in quadrature */
/*              rule */
/*      NQP     number of quadrature points (=7) */

/* ROUTINES called */
/*      ERRMES */

/*      SUBROUTINE QTRI7(WGHT,IWGHT,ABSS,IABSS,JABSS,NQP,ITEST) */
/* ***********************************************************************
 */

    /* Parameter adjustments */
    abss_dim1 = *iabss;
    abss_offset = abss_dim1 + 1;
    abss -= abss_offset;
    --wght;

    /* Function Body */

    *nqp = 7;

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

    area = sqrt(3.) * 3. / 4.;
    wght[1] = area * 3. / 60.;
    wght[2] = area * 8. / 60.;
    wght[3] = wght[1];
    wght[4] = wght[2];
    wght[5] = wght[1];
    wght[6] = wght[2];
    wght[7] = area * 27. / 60.;

    abss[abss_dim1 + 1] = 1.;
    abss[(abss_dim1 << 1) + 1] = .25;
    abss[abss_dim1 * 3 + 1] = -.5;
    abss[(abss_dim1 << 2) + 1] = -.5;
    abss[abss_dim1 * 5 + 1] = -.5;
    abss[abss_dim1 * 6 + 1] = .25;
    abss[abss_dim1 * 7 + 1] = 0.;
    abss[abss_dim1 + 2] = 0.;
    abss[(abss_dim1 << 1) + 2] = -sqrt(3.) / 4.;
    abss[abss_dim1 * 3 + 2] = -sqrt(3.) / 2.;
    abss[(abss_dim1 << 2) + 2] = 0.;
    abss[abss_dim1 * 5 + 2] = sqrt(3.) / 2.;
    abss[abss_dim1 * 6 + 2] = sqrt(3.) / 4.;
    abss[abss_dim1 * 7 + 2] = 0.;

} /* qtri7_ */

