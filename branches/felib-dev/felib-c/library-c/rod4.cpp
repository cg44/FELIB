/* rod4.f -- translated by f2c (version of 28 August 1991  0:07:02).
   You must link the resulting object file with the libraries:
	-lF77 -lI77 -lm -lc   (in that order)
*/

#include "f2c.h"


/* Subroutine */ int rod4_(fun, ifun, der, ider, xi, itest)
doublereal *fun;
integer *ifun;
doublereal *der;
integer *ider;
doublereal *xi;
integer *itest;
{
    /* Initialized data */

    static char srname[6+1] = "ROD4  ";

    extern doublereal veps_();
    static doublereal xmax, dummy;
    extern integer errmes_();
    static integer ierror;

/* -----------------------------------------------------------------------
 */
/* PURPOSE */
/*      ROD4 returns the values of the shape functions and their */
/*      derivatives at a specified point for a 4-noded c0 */
/*      continuous element. The function continuous across element */
/*      boundaries */

/* HISTORY */

/*      Copyright (C) 1996 : CCLRC, Rutherford Appleton Laboratory */
/*                           Chilton, Didcot, Oxfordshire OX11 0QX */

/*      Release 1.1  29 Oct 1979 (CG) */
/*      Commented    11 Oct 1980 (KR) */
/*      Release 4.0   2 Oct 1996 (CG) */

/* ARGUMENTS in */
/*      IFUN    dimension of vector FUN (.GE. 4) */
/*      IDER    dimension of vector DER (.GE. 4) */
/*      XI      value of local coordinate */
/*      ITEST   error checking option */

/* ARGUMENTS out */
/*      FUN     vector of dimension IFUN. FUN(I) contains the */
/*              value of the i'th shape function at XI */
/*      DER     vector of dimension IDER. DER(I) contains the */
/*              value of the derivative of the i'th shape */
/*              function at XI */

/* ROUTINES called */
/*      ERRMES */

/*     SUBROUTINE ROD4(FUN,IFUN,DER,IDER,XI,ITEST) */
/* ***********************************************************************
 */

    /* Parameter adjustments */
    --der;
    --fun;

    /* Function Body */

/*     Parameter checking */

    if (*itest != -1) {
	ierror = 0;
	if (*ifun < 4 || *ider < 4) {
	    ierror = 1;
	}
	xmax = veps_(&dummy) + 1.;
	if (abs(*xi) > xmax) {
	    ierror = 2;
	}
	*itest = errmes_(itest, &ierror, srname, 6L);
	if (*itest != 0) {
	    return 0;
	}
    }

/*     Main body */

    fun[1] = (1. - *xi) * .0625 * (*xi * 9. * *xi - 1.);
    fun[2] = (*xi * *xi - 1.) * .5625 * (*xi * 3. - 1.);
    fun[3] = (1. - *xi * *xi) * .5625 * (*xi * 3. + 1.);
    fun[4] = (*xi + 1.) * .0625 * (*xi * 9. * *xi - 1.);

    der[1] = (-1. - *xi * 18. + *xi * 27. * *xi) * -.0625;
    der[2] = (*xi * 9. * *xi - *xi * 2. - 3.) * .5625;
    der[3] = (3. - *xi * 2. - *xi * 9. * *xi) * .5625;
    der[4] = (*xi * 27. * *xi + *xi * 18. - 1.) * .0625;

} /* rod4_ */

