/* quam4.f -- translated by f2c (version of 28 August 1991  0:07:02).
   You must link the resulting object file with the libraries:
	-lF77 -lI77 -lm -lc   (in that order)
*/

#include "f2c.h"


/* Subroutine */ int quam4_(fun, ifun, der, ider, jder, xi, eta, itest)
doublereal *fun;
integer *ifun;
doublereal *der;
integer *ider, *jder;
doublereal *xi, *eta;
integer *itest;
{
    /* Initialized data */

    static char srname[6+1] = "QUAM4 ";

    /* System generated locals */
    integer der_dim1, der_offset;

    /* Local variables */
    static doublereal etam, etap;
    extern doublereal veps_();
    static doublereal dummy;
    extern integer errmes_();
    static integer ierror;
    static doublereal val, xim, xip;

/* -----------------------------------------------------------------------
 */
/* PURPOSE */
/*      QUAM4 returns the values of shape functions and their */
/*      derivatives at a specified point for a 4-noded c0 */
/*      continuous quadrilateral element.  The approximated */
/*      function will be continuous across element boundaries */

/* HISTORY */

/*      Copyright (C) 1996 : CCLRC, Rutherford Appleton Laboratory */
/*                           Chilton, Didcot, Oxfordshire OX11 0QX */

/*      Release 1.1  29 Oct 1979 (CG) */
/*      Commented    21 Oct 1980 (KR) */
/*      Release 4.0   2 Oct 1996 (CG) */

/* ARGUMENTS in */
/*      IFUN    length of vector FUN (.GE. 4) */
/*      IDER    first dimension of array DER (.GE. 2) */
/*      JDER    second dimension of array DER (.GE. 4) */
/*      XI      first local coordinate value */
/*      ETA     second local coordinate value */
/*      ITEST   error checking option */

/* ARGUMENTS out */
/*      FUN     vector of length IFUN. FUN(I) contains the */
/*              value of the i'th shape function at the point */
/*              (XI, ETA) */
/*      DER     array of dimension (IDER, JDER). DER(I, J) */
/*              contains the value of the derivative of the j'th */
/*              shape function with respect to the i'th */
/*              coordinate at the point (XI, ETA) */

/* ROUTINES called */
/*      ERRMES */

/*     SUBROUTINE QUAM4(FUN,IFUN,DER,IDER,JDER,XI,ETA,ITEST) */
/* ***********************************************************************
 */

    /* Parameter adjustments */
    der_dim1 = *ider;
    der_offset = der_dim1 + 1;
    der -= der_offset;
    --fun;

    /* Function Body */

/*     Parameter checking */

    if (*itest != -1) {
	ierror = 0;
	if (*ifun < 4) {
	    ierror = 1;
	}
	if (*ider < 2 || *jder < 4) {
	    ierror = 2;
	}
	val = veps_(&dummy) + 1.;
	if (abs(*xi) > val || abs(*eta) > val) {
	    ierror = 3;
	}
	*itest = errmes_(itest, &ierror, srname, 6L);
	if (*itest != 0) {
	    return 0;
	}
    }

/*     Main body */

    etam = (1. - *eta) * .25;
    etap = (*eta + 1.) * .25;
    xim = (1. - *xi) * .25;
    xip = (*xi + 1.) * .25;

    fun[1] = xim * 4. * etam;
    fun[2] = xim * 4. * etap;
    fun[3] = xip * 4. * etap;
    fun[4] = xip * 4. * etam;

    der[der_dim1 + 1] = -etam;
    der[der_dim1 + 2] = -xim;
    der[(der_dim1 << 1) + 1] = -etap;
    der[(der_dim1 << 1) + 2] = xim;
    der[der_dim1 * 3 + 1] = etap;
    der[der_dim1 * 3 + 2] = xip;
    der[(der_dim1 << 2) + 1] = etam;
    der[(der_dim1 << 2) + 2] = -xip;

} /* quam4_ */

