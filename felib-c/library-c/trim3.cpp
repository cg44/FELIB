/* trim3.f -- translated by f2c (version of 28 August 1991  0:07:02).
   You must link the resulting object file with the libraries:
	-lF77 -lI77 -lm -lc   (in that order)
*/

#include "f2c.h"


/* Subroutine */ int trim3_(fun, ifun, der, ider, jder, xi, eta, itest)
doublereal *fun;
integer *ifun;
doublereal *der;
integer *ider, *jder;
doublereal *xi, *eta;
integer *itest;
{
    /* Initialized data */

    static char srname[6+1] = "TRIM3 ";

    /* System generated locals */
    integer der_dim1, der_offset;

    /* Builtin functions */
    double sqrt();

    /* Local variables */
    static doublereal xmin, ymin;
    extern doublereal veps_();
    static doublereal xmax, ymax, dummy;
    extern integer errmes_();
    static integer ierror;

/* -----------------------------------------------------------------------
 */
/* PURPOSE */
/*      TRIM3 returns the values of the shape functions and their */
/*      derivatives at a specified point for a three-noded c0 */
/*      continuous triangular element. The approximation function */
/*      continuous across element boundaries. */

/* HISTORY */

/*      Copyright (C) 1996 : CCLRC, Rutherford Appleton Laboratory */
/*                           Chilton, Didcot, Oxfordshire OX11 0QX */

/*      Release 1.1  29 Oct 1979 (CG) */
/*      Commented    22 Oct 1980 (KR) */
/*      Release 4.0   2 Oct 1996 (CG) */

/* ARGUMENTS in */
/*      IFUN    length of vector FUN (.GE. 3) */
/*      IDER    first dimension of array DER (.GE. 2) */
/*      JDER    second dimension of array DER (.GE. 3) */
/*      XI      first local coordinate */
/*      ETA     second local coordinate */
/*      ITEST   error checking option */

/* ARGUMENTS out */
/*      FUN     vector of length IFUN. FUN(I) contains the */
/*              value of the i'th shape function at the point */
/*              (XI,ETA), for i=1(1)3 */
/*      DER     array of dimension (IDER, JDER).  DER(I, J) */
/*              contains the derivative of the j'th shape */
/*              function with respect to the i'th coordinate for */
/*              i=1(1)2 and j=1(1)3 */

/* ROUTINES called */
/*      ERRMES */

/*     SUBROUTINE TRIM3(FUN,IFUN,DER,IDER,JDER,XI,ETA,ITEST) */
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
	if (*ifun < 3) {
	    ierror = 1;
	}
	if (*ider < 2 || *jder < 3) {
	    ierror = 2;
	}
	ymin = 1. / sqrt(3.) * (*xi - 1.) - veps_(&dummy);
	ymax = 1. / sqrt(3.) * (1. - *xi) + veps_(&dummy);
	xmin = -(veps_(&dummy) + .5);
	xmax = veps_(&dummy) + 1.;
	if (*xi < xmin || *xi > xmax || (*eta < ymin || *eta > ymax)) {
	    ierror = 3;
	}
	*itest = errmes_(itest, &ierror, srname, 6L);
	if (*itest != 0) {
	    return 0;
	}
    }

/*     Main body */

    fun[1] = (*xi * 2. + 1.) * .33333333333333331;
    fun[2] = (1. - *xi - sqrt(3.) * *eta) * .33333333333333331;
    fun[3] = (1. - *xi + sqrt(3.) * *eta) * .33333333333333331;

/*     Derivatives */

    der[der_dim1 + 1] = .66666666666666663;
    der[(der_dim1 << 1) + 1] = -.33333333333333331;
    der[der_dim1 * 3 + 1] = -.33333333333333331;
    der[der_dim1 + 2] = 0.;
    der[(der_dim1 << 1) + 2] = -1. / sqrt(3.);
    der[der_dim1 * 3 + 2] = 1. / sqrt(3.);

} /* trim3_ */

