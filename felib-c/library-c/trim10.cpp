/* trim10.f -- translated by f2c (version of 28 August 1991  0:07:02).
   You must link the resulting object file with the libraries:
	-lF77 -lI77 -lm -lc   (in that order)
*/

#include "f2c.h"


/* Subroutine */ int trim10_(fun, ifun, der, ider, jder, xi, eta, itest)
doublereal *fun;
integer *ifun;
doublereal *der;
integer *ider, *jder;
doublereal *xi, *eta;
integer *itest;
{
    /* Initialized data */

    static char srname[6+1] = "TRIM10";

    /* System generated locals */
    integer der_dim1, der_offset;

    /* Builtin functions */
    double sqrt();

    /* Local variables */
    static doublereal xmin, ymin;
    extern doublereal veps_();
    static doublereal xmax, ymax, dummy, l1, l2, l3;
    extern integer errmes_();
    static integer ierror;
    static doublereal dlx[3], dly[3];

/* -----------------------------------------------------------------------
 */
/* PURPOSE */
/*      TRIM10 returns the values of the shape functions and their */
/*      derivatives at a specified point for a 10-noded c0 */
/*      continuous triangular element. The function is */
/*      continuous across element boundaries. */

/* HISTORY */

/*      Copyright (C) 1996 : CCLRC, Rutherford Appleton Laboratory */
/*                           Chilton, Didcot, Oxfordshire OX11 0QX */

/*      Release 1.1  29 Oct 1979 (CG) */
/*      Commented    22 Oct 1980 (KR) */
/*      Release 4.0   2 Oct 1996 (CG) */

/* ARGUMENTS in */
/*      IFUN    length of vector FUN (.GE. 10) */
/*      IDER    first dimension of array DER (.GE. 2) */
/*      JDER    second dimension of array DER (.GE. 10) */
/*      XI      first local coordinate */
/*      ETA     second local coordinate */
/*      ITEST   error checking option */

/* ARGUMENTS out */
/*      FUN     vector of length IFUN. FUN(I) contains the */
/*              value of the i'th shape function at the point */
/*              (XI,ETA), for i=1(1)10 */
/*      DER     array of dimension (IDER, JDER). DER(I, J) */
/*              contains the value of the derivative of the j'th */
/*              shape function with respect to the i'th local */
/*              coordinate, for i=1(1)2 and j=1(1)10 */

/* ROUTINES called */
/*      ERRMES */

/*     SUBROUTINE TRIM10(FUN,IFUN,DER,IDER,JDER,XI,ETA,ITEST) */
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
	if (*ifun < 10) {
	    ierror = 1;
	}
	if (*ider < 2 || *jder < 10) {
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

    l1 = (*xi * 2. + 1.) * .33333333333333331;
    l2 = (1. - *xi - sqrt(3.) * *eta) * .33333333333333331;
    l3 = (1. - *xi + sqrt(3.) * *eta) * .33333333333333331;
    dlx[0] = .66666666666666663;
    dly[0] = 0.;
    dlx[1] = -.33333333333333331;
    dly[1] = -1. / sqrt(3.);
    dlx[2] = -.33333333333333331;
    dly[2] = 1. / sqrt(3.);

/*     Shape functions */

    fun[1] = (l1 * 3. - 1.) * .5 * (l1 * 3. - 2.) * l1;
    fun[2] = l1 * 4.5 * l2 * (l1 * 3. - 1.);
    fun[3] = l2 * 4.5 * l1 * (l2 * 3. - 1.);
    fun[4] = (l2 * 3. - 1.) * .5 * (l2 * 3. - 2.) * l2;
    fun[5] = l2 * 4.5 * l3 * (l2 * 3. - 1.);
    fun[6] = l3 * 4.5 * l2 * (l3 * 3. - 1.);
    fun[7] = (l3 * 3. - 1.) * .5 * (l3 * 3. - 2.) * l3;
    fun[8] = l3 * 4.5 * l1 * (l3 * 3. - 1.);
    fun[9] = l1 * 4.5 * l3 * (l1 * 3. - 1.);
    fun[10] = l1 * 27. * l2 * l3;

/*     Derivatives */

    der[der_dim1 + 1] = (l1 * 27. * l1 - l1 * 18. + 2.) * .5 * dlx[0];
    der[der_dim1 + 2] = (l1 * 27. * l1 - l1 * 18. + 2.) * .5 * dly[0];
    der[(der_dim1 << 1) + 1] = l2 * 4.5 * (l1 * 6. - 1.) * dlx[0] + l1 * 4.5 *
	     (l1 * 3. - 1.) * dlx[1];
    der[(der_dim1 << 1) + 2] = l2 * 4.5 * (l1 * 6. - 1.) * dly[0] + l1 * 4.5 *
	     (l1 * 3. - 1.) * dly[1];
    der[der_dim1 * 3 + 1] = l1 * 4.5 * (l2 * 6. - 1.) * dlx[1] + l2 * 4.5 * (
	    l2 * 3. - 1.) * dlx[0];
    der[der_dim1 * 3 + 2] = l1 * 4.5 * (l2 * 6. - 1.) * dly[1] + l2 * 4.5 * (
	    l2 * 3. - 1.) * dly[0];
    der[(der_dim1 << 2) + 1] = (l2 * 27. * l2 - l2 * 18. + 2.) * .5 * dlx[1];
    der[(der_dim1 << 2) + 2] = (l2 * 27. * l2 - l2 * 18. + 2.) * .5 * dly[1];
    der[der_dim1 * 5 + 1] = l3 * 4.5 * (l2 * 6. - 1.) * dlx[1] + l2 * 4.5 * (
	    l2 * 3. - 1.) * dlx[2];
    der[der_dim1 * 5 + 2] = l3 * 4.5 * (l2 * 6. - 1.) * dly[1] + l2 * 4.5 * (
	    l2 * 3. - 1.) * dly[2];
    der[der_dim1 * 6 + 1] = l2 * 4.5 * (l3 * 6. - 1.) * dlx[2] + l3 * 4.5 * (
	    l3 * 3. - 1.) * dlx[1];
    der[der_dim1 * 6 + 2] = l2 * 4.5 * (l3 * 6. - 1.) * dly[2] + l3 * 4.5 * (
	    l3 * 3. - 1.) * dly[1];
    der[der_dim1 * 7 + 1] = (l3 * 27. * l3 - l3 * 18. + 2.) * .5 * dlx[2];
    der[der_dim1 * 7 + 2] = (l3 * 27. * l3 - l3 * 18. + 2.) * .5 * dly[2];
    der[(der_dim1 << 3) + 1] = l1 * 4.5 * (l3 * 6. - 1.) * dlx[2] + l3 * 4.5 *
	     (l3 * 3. - 1.) * dlx[0];
    der[(der_dim1 << 3) + 2] = l1 * 4.5 * (l3 * 6. - 1.) * dly[2] + l3 * 4.5 *
	     (l3 * 3. - 1.) * dly[0];
    der[der_dim1 * 9 + 1] = l3 * 4.5 * (l1 * 6. - 1.) * dlx[0] + l1 * 4.5 * (
	    l1 * 3. - 1.) * dlx[2];
    der[der_dim1 * 9 + 2] = l3 * 4.5 * (l1 * 6. - 1.) * dly[0] + l1 * 4.5 * (
	    l1 * 3. - 1.) * dly[2];
    der[der_dim1 * 10 + 1] = (dlx[0] * l2 * l3 + l1 * dlx[1] * l3 + l1 * l2 * 
	    dlx[2]) * 27.;
    der[der_dim1 * 10 + 2] = (dly[0] * l2 * l3 + l1 * dly[1] * l3 + l1 * l2 * 
	    dly[2]) * 27.;

} /* trim10_ */

