/* rod2.f -- translated by f2c (version of 28 August 1991  0:07:02).
   You must link the resulting object file with the libraries:
	-lF77 -lI77 -lm -lc   (in that order)
*/

#include "f2c.h"


/* Subroutine */ int rod2_(fun, ifun, der, ider, xi, itest)
doublereal *fun;
integer *ifun;
doublereal *der;
integer *ider;
doublereal *xi;
integer *itest;
{
    /* Initialized data */

    static char srname[6+1] = "ROD2  ";

    extern doublereal veps_();
    static doublereal xmax, dummy;
    extern integer errmes_();
    static integer ierror;

/* -----------------------------------------------------------------------
 */
/* PURPOSE */
/*      ROD2 returns the values of the shape functions and their */
/*      derivatives at a specified point for a 2-noded c0 */
/*      continuous line element.  Only the function will be */
/*      continuous across element boundaries. */

/* HISTORY */

/*      Copyright (C) 1996 : CCLRC, Rutherford Appleton Laboratory */
/*                           Chilton, Didcot, Oxfordshire OX11 0QX */

/*      Release 1.1  29 Oct 1979 (CG) */
/*      Commented    26 Feb 1980 (KR) */
/*      Release 4.0   2 Oct 1996 (CG) */

/* ARGUMENTS in */
/*      IFUN    dimension of the vector FUN (.GE. 2) */
/*      IDER    dimension of the vector DER (.GE. 2) */
/*      XI      specifies the value of the local coordinate at */
/*              which the function and its derivative are */
/*              required */
/*      ITEST   error checking option */

/* ARGUMENTS out */
/*      FUN     FUN(I) contains the value of the i'th shape */
/*              function at XI */
/*      DER     DER(I) contains the value of the derivative of */
/*              the i'th shape function at XI */

/* ROUTINES called */
/*      ERRMES */

/*     SUBROUTINE ROD2(FUN,IFUN,DER,IDER,XI,ITEST) */
/* ***********************************************************************
 */

    /* Parameter adjustments */
    --der;
    --fun;

    /* Function Body */

/*     Parameter checking */

    if (*itest != -1) {
	ierror = 0;
	if (*ifun < 2 || *ider < 2) {
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

    fun[1] = (1. - *xi) * .5;
    fun[2] = (*xi + 1.) * .5;
    der[1] = -.5;
    der[2] = .5;

} /* rod2_ */

