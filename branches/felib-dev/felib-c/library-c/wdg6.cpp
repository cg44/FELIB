/* wdg6.f -- translated by f2c (version of 28 August 1991  0:07:02).
   You must link the resulting object file with the libraries:
	-lF77 -lI77 -lm -lc   (in that order)
*/

#include "f2c.h"


/* Subroutine */ int wdg6_(fun, ifun, der, ider, jder, xi, eta, zeta, itest)
doublereal *fun;
integer *ifun;
doublereal *der;
integer *ider, *jder;
doublereal *xi, *eta, *zeta;
integer *itest;
{
    /* Initialized data */

    static char srname[6+1] = "WDG6  ";

    /* System generated locals */
    integer der_dim1, der_offset;

    /* Builtin functions */
    double sqrt();

    /* Local variables */
    static doublereal xmin, ymin;
    extern doublereal veps_();
    static doublereal xmax, ymax, zval, zetam, zetap, dummy, l1, l2, l3;
    extern integer errmes_();
    static integer ierror;

/* -----------------------------------------------------------------------
 */
/* PURPOSE */
/*      WDG6 returns the values of shape functions and derivatives */
/*      at a specified point for an 6-noded pentahedral element. */
/*      the function is continuous across element boundaries. */

/* HISTORY */

/*      Copyright (C) 1996 : CCLRC, Rutherford Appleton Laboratory */
/*                           Chilton, Didcot, Oxfordshire OX11 0QX */

/*      Release 1.1  29 Oct 1979 (CG) */
/*      Commented    10 Oct 1980 (KR) */
/*      Release 4.0   2 Oct 1996 (CG) */

/* ARGUMENTS in */
/*      IFUN    dimension of vector FUN (.GE. 6) */
/*      IDER    first dimension of array DER (.GE. 3) */
/*      JDER    second dimension of array DER (.GE. 6) */
/*      XI      value of local coordinate at which function and */
/*              derivative values required */
/*      ETA     value of second local coordinate */
/*      ZETA    value of third local coordinate */
/*      ITEST   error checking option */

/* ARGUMENTS out */
/*      FUN     real vector of dimension IFUN. FUN(I) contains */
/*              value of i'th shape function at (XI, ETA, ZETA) */
/*      DER     real array of dimensions (IDER, JDER).  DER(I, J) */
/*              contains the derivative of the j'th shape */
/*              function with respect to the i'th coordinate at */
/*              the point (XI, ETA, ZETA) */

/* ROUTINES called */
/*      ERRMES */

/*     SUBROUTINE WDG6(FUN,IFUN,DER,IDER,JDER,XI,ETA,ZETA,ITEST) */
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
	if (*ifun < 6) {
	    ierror = 1;
	}
	if (*ider < 3 || *jder < 6) {
	    ierror = 2;
	}
	ymin = 1. / sqrt(3.) * (*xi - 1.) - veps_(&dummy);
	ymax = 1. / sqrt(3.) * (1. - *xi) + veps_(&dummy);
	xmin = -(veps_(&dummy) + .5);
	xmax = veps_(&dummy) + 1.;
	zval = veps_(&dummy) + 1.;
	if (*xi < xmin || *xi > xmax || (*eta < ymin || *eta > ymax) || abs(*
		zeta) > zval) {
	    ierror = 3;
	}
	*itest = errmes_(itest, &ierror, srname, 6L);
	if (*itest != 0) {
	    return 0;
	}
    }

/*     Main body */

    l1 = (*xi * 2. + 1.) * .16666666666666666;
    l2 = (1. - *xi - sqrt(3.) * *eta) * .16666666666666666;
    l3 = (1. - *xi + sqrt(3.) * *eta) * .16666666666666666;
    zetam = 1. - *zeta;
    zetap = *zeta + 1.;

/*     Shape functions */

    fun[1] = l1 * zetam;
    fun[2] = l2 * zetam;
    fun[3] = l3 * zetam;
    fun[4] = l1 * zetap;
    fun[5] = l2 * zetap;
    fun[6] = l3 * zetap;

/*     Derivatives */

    der[der_dim1 + 1] = zetam * .33333333333333331;
    der[der_dim1 + 2] = 0.;
    der[der_dim1 + 3] = -l1;
    der[(der_dim1 << 1) + 1] = zetam * -.16666666666666666;
    der[(der_dim1 << 1) + 2] = -1. / (sqrt(3.) * 2.) * zetam;
    der[(der_dim1 << 1) + 3] = -l2;
    der[der_dim1 * 3 + 1] = zetam * -.16666666666666666;
    der[der_dim1 * 3 + 2] = 1. / (sqrt(3.) * 2.) * zetam;
    der[der_dim1 * 3 + 3] = -l3;
    der[(der_dim1 << 2) + 1] = zetap * .33333333333333331;
    der[(der_dim1 << 2) + 2] = 0.;
    der[(der_dim1 << 2) + 3] = l1;
    der[der_dim1 * 5 + 1] = zetap * -.16666666666666666;
    der[der_dim1 * 5 + 2] = -1. / (sqrt(3.) * 2.) * zetap;
    der[der_dim1 * 5 + 3] = l2;
    der[der_dim1 * 6 + 1] = zetap * -.16666666666666666;
    der[der_dim1 * 6 + 2] = 1. / (sqrt(3.) * 2.) * zetap;
    der[der_dim1 * 6 + 3] = l3;

} /* wdg6_ */

