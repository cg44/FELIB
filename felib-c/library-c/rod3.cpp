/* rod3.f -- translated by f2c (version of 28 August 1991  0:07:02).
   You must link the resulting object file with the libraries:
	-lF77 -lI77 -lm -lc   (in that order)
*/

#include "f2c.h"


/* Subroutine */ int rod3_(fun, ifun, der, ider, xi, itest)
doublereal *fun;
integer *ifun;
doublereal *der;
integer *ider;
doublereal *xi;
integer *itest;
{
    /* Initialized data */

    static char srname[6+1] = "ROD3  ";

    extern doublereal veps_();
    static doublereal xmax, dummy;
    extern integer errmes_();
    static integer ierror;

/* -----------------------------------------------------------------------
 */
/* PURPOSE */
/*      ROD3 returns the values of the shape function and its */
/*      derivative at a specified point for a 3-noded c0 */
/*      continuous line element. The function continuous across */
/*      element boundaries */

/* HISTORY */

/*      Copyright (C) 1996 : CCLRC, Rutherford Appleton Laboratory */
/*                           Chilton, Didcot, Oxfordshire OX11 0QX */

/*      Release 1.1  29 Oct 1979 (CG) */
/*      Commented    11 Oct 1980 (KR) */
/*      Release 4.0   2 Oct 1996 (CG) */

/* ARGUMENTS in */
/*      IFUN    dimension of vector FUN (.GE. 3) */
/*      IDER    dimension of vector DER (.GE. 3) */
/*      XI      value of local coordinate */
/*      ITEST   error checking option */

/* ARGUMENTS out */
/*      FUN     vector of length IFUN.  FUN(I) contains the */
/*              value of the i'th shape function at the point XI */
/*      DER     vector of length IDER.  DER(I) contains the */
/*              value of the derivative of the i'th shape */
/*              function at the point XI */

/* ROUTINES called */
/*      ERRMES */

/*     SUBROUTINE ROD3(FUN,IFUN,DER,IDER,XI,ITEST) */
/* ***********************************************************************
 */

    /* Parameter adjustments */
    --der;
    --fun;

    /* Function Body */

/*     Parameter checking */

    if (*itest != -1) {
	ierror = 0;
	if (*ifun < 3 || *ider < 3) {
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

    fun[1] = *xi * .5 * (*xi - 1.);
    fun[2] = 1. - *xi * *xi;
    fun[3] = *xi * .5 * (*xi + 1.);

    der[1] = (*xi * 2. - 1.) * .5;
    der[2] = *xi * -2.;
    der[3] = (*xi * 2. + 1.) * .5;

} /* rod3_ */

