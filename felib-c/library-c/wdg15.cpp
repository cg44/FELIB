/* wdg15.f -- translated by f2c (version of 28 August 1991  0:07:02).
   You must link the resulting object file with the libraries:
	-lF77 -lI77 -lm -lc   (in that order)
*/

#include "f2c.h"


/* Subroutine */ int wdg15_(fun, ifun, der, ider, jder, xi, eta, zeta, itest)
doublereal *fun;
integer *ifun;
doublereal *der;
integer *ider, *jder;
doublereal *xi, *eta, *zeta;
integer *itest;
{
    /* Initialized data */

    static char srname[6+1] = "WDG15 ";

    /* System generated locals */
    integer der_dim1, der_offset;

    /* Builtin functions */
    double sqrt();

    /* Local variables */
    static doublereal xmin, ymin;
    extern doublereal veps_();
    static doublereal xmax, ymax, zval, dummy, l1, l2, l3, zm, zp;
    extern integer errmes_();
    static integer ierror;
    static doublereal dl1x, dl1y, dl2x, dl2y, dl3x, dl3y;

/* -----------------------------------------------------------------------
 */
/* PURPOSE */
/*      WDG15 returns the values of shape functions and derivatives */
/*      at a specified point for an 6-noded pentahedral element. */
/*      the function is continuous across element boundaries. */

/* HISTORY */

/*      Copyright (C) 1996 : CCLRC, Rutherford Appleton Laboratory */
/*                           Chilton, Didcot, Oxfordshire OX11 0QX */

/*      Release 1.1  29 Oct 1979 (CG) */
/*      Commented    10 Oct 1980 (KR) */
/*      Release 4.0   2 Oct 1996 (CG) */

/* ARGUMENTS in */
/*      IFUN    dimension of vector FUN (.GE. 15) */
/*      IDER    first dimension of array DER (.GE. 3) */
/*      JDER    second dimension of array DER (.GE. 15) */
/*      XI      value of local coordinate at which function and */
/*              derivative values required */
/*      ETA     value of second local coordinate */
/*      ZETA    value of third local coordinate */
/*      ITEST   error checking option */

/* ARGUMENTS out */
/*      FUN     real vector of dimension IFUN. FUN(I) contains */
/*              value of i'th shape function at (XI, ETA, ZETA) */
/*      DER     real array of dimensions (IDER, JDER). DER(I, J) */
/*              contains the derivative of the j'th shape */
/*              function with respect to the i'th coordinate at */
/*              the point (XI, ETA, ZETA) */

/* ROUTINES called */
/*      ERRMES */

/*     SUBROUTINE WDG15(FUN,IFUN,DER,IDER,JDER,XI,ETA,ZETA,ITEST) */
/* ***********************************************************************
 */


    /* Parameter adjustments */
    der_dim1 = *ider;
    der_offset = der_dim1 + 1;
    der -= der_offset;
    --fun;

    /* Function Body */

/*     Statement functions */


/*     Parameter checking */

    if (*itest != -1) {
	ierror = 0;
	ymin = 1. / sqrt(3.) * (*xi - 1.) - veps_(&dummy);
	ymax = 1. / sqrt(3.) * (1. - *xi) + veps_(&dummy);
	xmin = -(veps_(&dummy) + .5);
	xmax = veps_(&dummy) + 1.;
	zval = veps_(&dummy) + 1.;
	if (*xi < xmin || *xi > xmax || (*eta < ymin || *eta > ymax) || abs(*
		zeta) > zval) {
	    ierror = 3;
	}
	if (*ider < 3 || *jder < 15) {
	    ierror = 2;
	}
	if (*ifun < 15) {
	    ierror = 1;
	}
	*itest = errmes_(itest, &ierror, srname, 6L);
	if (*itest != 0) {
	    return 0;
	}
    }

/*     Main body */

    l1 = (*xi * 2. + 1.) * .33333333333333331;
    l2 = (1. - *xi - sqrt(3.) * *eta) * .33333333333333331;
    l3 = (1. - *xi + sqrt(3.) * *eta) * .33333333333333331;
    dl1x = .66666666666666663;
    dl1y = 0.;
    dl2x = -.33333333333333331;
    dl2y = -1. / sqrt(3.);
    dl3x = -.33333333333333331;
    dl3y = 1. / sqrt(3.);
    zp = 1.;
    zm = -1.;

/*     Shape functions */

    fun[1] = l1 * .5 * (l1 * 2. - 1.) * (zm * *zeta + 1.) - l1 * .5 * (1. - *
	    zeta * *zeta);
    fun[2] = l1 * 2. * l2 * (zm * *zeta + 1.);
    fun[3] = l2 * .5 * (l2 * 2. - 1.) * (zm * *zeta + 1.) - l2 * .5 * (1. - *
	    zeta * *zeta);
    fun[4] = l2 * 2. * l3 * (zm * *zeta + 1.);
    fun[5] = l3 * .5 * (l3 * 2. - 1.) * (zm * *zeta + 1.) - l3 * .5 * (1. - *
	    zeta * *zeta);
    fun[6] = l3 * 2. * l1 * (zm * *zeta + 1.);
    fun[7] = l1 * (1. - *zeta * *zeta);
    fun[8] = l2 * (1. - *zeta * *zeta);
    fun[9] = l3 * (1. - *zeta * *zeta);
    fun[10] = l1 * .5 * (l1 * 2. - 1.) * (zp * *zeta + 1.) - l1 * .5 * (1. - *
	    zeta * *zeta);
    fun[11] = l1 * 2. * l2 * (zp * *zeta + 1.);
    fun[12] = l2 * .5 * (l2 * 2. - 1.) * (zp * *zeta + 1.) - l2 * .5 * (1. - *
	    zeta * *zeta);
    fun[13] = l2 * 2. * l3 * (zp * *zeta + 1.);
    fun[14] = l3 * .5 * (l3 * 2. - 1.) * (zp * *zeta + 1.) - l3 * .5 * (1. - *
	    zeta * *zeta);
    fun[15] = l3 * 2. * l1 * (zp * *zeta + 1.);

/*     Derivatives */

    der[der_dim1 + 1] = ((l1 * 4. - 1.) * (zm * *zeta + 1.) + *zeta * *zeta - 
	    1.) * .5 * dl1x;
    der[der_dim1 + 2] = ((l1 * 4. - 1.) * (zm * *zeta + 1.) + *zeta * *zeta - 
	    1.) * .5 * dl1y;
    der[der_dim1 + 3] = l1 * .5 * (*zeta * 2. + zm * (l1 * 2. - 1.));
    der[(der_dim1 << 1) + 1] = l2 * 2. * (zm * *zeta + 1.) * dl1x + l1 * 2. * 
	    (zm * *zeta + 1.) * dl2x;
    der[(der_dim1 << 1) + 2] = l2 * 2. * (zm * *zeta + 1.) * dl1y + l1 * 2. * 
	    (zm * *zeta + 1.) * dl2y;
    der[(der_dim1 << 1) + 3] = l1 * 2. * l2 * zm;
    der[der_dim1 * 3 + 1] = ((l2 * 4. - 1.) * (zm * *zeta + 1.) + *zeta * *
	    zeta - 1.) * .5 * dl2x;
    der[der_dim1 * 3 + 2] = ((l2 * 4. - 1.) * (zm * *zeta + 1.) + *zeta * *
	    zeta - 1.) * .5 * dl2y;
    der[der_dim1 * 3 + 3] = l2 * .5 * (*zeta * 2. + zm * (l2 * 2. - 1.));
    der[(der_dim1 << 2) + 1] = l3 * 2. * (zm * *zeta + 1.) * dl2x + l2 * 2. * 
	    (zm * *zeta + 1.) * dl3x;
    der[(der_dim1 << 2) + 2] = l3 * 2. * (zm * *zeta + 1.) * dl2y + l2 * 2. * 
	    (zm * *zeta + 1.) * dl3y;
    der[(der_dim1 << 2) + 3] = l2 * 2. * l3 * zm;
    der[der_dim1 * 5 + 1] = ((l3 * 4. - 1.) * (zm * *zeta + 1.) + *zeta * *
	    zeta - 1.) * .5 * dl3x;
    der[der_dim1 * 5 + 2] = ((l3 * 4. - 1.) * (zm * *zeta + 1.) + *zeta * *
	    zeta - 1.) * .5 * dl3y;
    der[der_dim1 * 5 + 3] = l3 * .5 * (*zeta * 2. + zm * (l3 * 2. - 1.));
    der[der_dim1 * 6 + 1] = l1 * 2. * (zm * *zeta + 1.) * dl3x + l3 * 2. * (
	    zm * *zeta + 1.) * dl1x;
    der[der_dim1 * 6 + 2] = l1 * 2. * (zm * *zeta + 1.) * dl3y + l3 * 2. * (
	    zm * *zeta + 1.) * dl1y;
    der[der_dim1 * 6 + 3] = l3 * 2. * l1 * zm;
    der[der_dim1 * 7 + 1] = (1. - *zeta * *zeta) * dl1x;
    der[der_dim1 * 7 + 2] = (1. - *zeta * *zeta) * dl1y;
    der[der_dim1 * 7 + 3] = l1 * -2. * *zeta;
    der[(der_dim1 << 3) + 1] = (1. - *zeta * *zeta) * dl2x;
    der[(der_dim1 << 3) + 2] = (1. - *zeta * *zeta) * dl2y;
    der[(der_dim1 << 3) + 3] = l2 * -2. * *zeta;
    der[der_dim1 * 9 + 1] = (1. - *zeta * *zeta) * dl3x;
    der[der_dim1 * 9 + 2] = (1. - *zeta * *zeta) * dl3y;
    der[der_dim1 * 9 + 3] = l3 * -2. * *zeta;
    der[der_dim1 * 10 + 1] = ((l1 * 4. - 1.) * (zp * *zeta + 1.) + *zeta * *
	    zeta - 1.) * .5 * dl1x;
    der[der_dim1 * 10 + 2] = ((l1 * 4. - 1.) * (zp * *zeta + 1.) + *zeta * *
	    zeta - 1.) * .5 * dl1y;
    der[der_dim1 * 10 + 3] = l1 * .5 * (*zeta * 2. + zp * (l1 * 2. - 1.));
    der[der_dim1 * 11 + 1] = l2 * 2. * (zp * *zeta + 1.) * dl1x + l1 * 2. * (
	    zp * *zeta + 1.) * dl2x;
    der[der_dim1 * 11 + 2] = l2 * 2. * (zp * *zeta + 1.) * dl1y + l1 * 2. * (
	    zp * *zeta + 1.) * dl2y;
    der[der_dim1 * 11 + 3] = l1 * 2. * l2 * zp;
    der[der_dim1 * 12 + 1] = ((l2 * 4. - 1.) * (zp * *zeta + 1.) + *zeta * *
	    zeta - 1.) * .5 * dl2x;
    der[der_dim1 * 12 + 2] = ((l2 * 4. - 1.) * (zp * *zeta + 1.) + *zeta * *
	    zeta - 1.) * .5 * dl2y;
    der[der_dim1 * 12 + 3] = l2 * .5 * (*zeta * 2. + zp * (l2 * 2. - 1.));
    der[der_dim1 * 13 + 1] = l3 * 2. * (zp * *zeta + 1.) * dl2x + l2 * 2. * (
	    zp * *zeta + 1.) * dl3x;
    der[der_dim1 * 13 + 2] = l3 * 2. * (zp * *zeta + 1.) * dl2y + l2 * 2. * (
	    zp * *zeta + 1.) * dl3y;
    der[der_dim1 * 13 + 3] = l2 * 2. * l3 * zp;
    der[der_dim1 * 14 + 1] = ((l3 * 4. - 1.) * (zp * *zeta + 1.) + *zeta * *
	    zeta - 1.) * .5 * dl3x;
    der[der_dim1 * 14 + 2] = ((l3 * 4. - 1.) * (zp * *zeta + 1.) + *zeta * *
	    zeta - 1.) * .5 * dl3y;
    der[der_dim1 * 14 + 3] = l3 * .5 * (*zeta * 2. + zp * (l3 * 2. - 1.));
    der[der_dim1 * 15 + 1] = l1 * 2. * (zp * *zeta + 1.) * dl3x + l3 * 2. * (
	    zp * *zeta + 1.) * dl1x;
    der[der_dim1 * 15 + 2] = l1 * 2. * (zp * *zeta + 1.) * dl3y + l3 * 2. * (
	    zp * *zeta + 1.) * dl1y;
    der[der_dim1 * 15 + 3] = l3 * 2. * l1 * zp;

} /* wdg15_ */

