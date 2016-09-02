/* brk20.f -- translated by f2c (version of 28 August 1991  0:07:02).
   You must link the resulting object file with the libraries:
	-lF77 -lI77 -lm -lc   (in that order)
*/

#include "f2c.h"


/* Subroutine */ int brk20_(fun, ifun, der, ider, jder, xi, eta, zeta, itest)
doublereal *fun;
integer *ifun;
doublereal *der;
integer *ider, *jder;
doublereal *xi, *eta, *zeta;
integer *itest;
{
    /* Initialized data */

    static doublereal m1 = -1.;
    static doublereal p1 = 1.;
    static char srname[6+1] = "BRK20 ";

    /* System generated locals */
    integer der_dim1, der_offset;

    /* Local variables */
    extern doublereal veps_();
    static doublereal dummy;
    extern integer errmes_();
    static integer ierror;
    static doublereal val;

/* -----------------------------------------------------------------------
 */
/* PURPOSE */
/*      BRK20 calculates shape functions and derivatives at a */
/*      specified point for 20-noded brick element. The shape function */
/*      is continuous across element boundaries. */

/* HISTORY */

/*      Copyright (C) 1996 : CCLRC, Rutherford Appleton Laboratory */
/*                           Chilton, Didcot, Oxfordshire OX11 0QX */

/*      Release 1.1  29 Oct 1979 (CG) */
/*      Commented    10 Oct 1980 (KR) */
/*      Release 4.0   2 Oct 1996 (CG) */

/* ARGUMENTS in */
/*      IFUN    dimension of vector FUN (.GE. 20) */
/*      IDER    first dimension of real array DER (.GE. 3) */
/*      JDER    second dimension of real array DER (.GE. 20) */
/*      XI      value of first local coordinate */
/*      ETA     value of second local coordinate */
/*      ZETA    value of third local coordinate */
/*      ITEST   error checking option */

/* ARGUMENTS out */
/*      FUN     vector of shape function values.  FUN(I) */
/*              contains the value of the i'th shape function at */
/*              (XI, ETA, ZETA) */
/*      DER     array of shape function derivative values. */
/*              DER(I, J) contains the value of the derivative of */
/*              the j'th shape function with respect to the i'th */
/*              coordinate at (XI, ETA, ZETA) */

/* ROUTINES called */
/*      ERRMES    VEPS */

/*      SUBROUTINE BRK20(FUN,IFUN,DER,IDER,JDER,XI,ETA,ZETA,ITEST) */
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
	if (*ifun < 20) {
	    ierror = 1;
	}
	if (*ider < 3 || *jder < 20) {
	    ierror = 2;
	}
	val = veps_(&dummy) + 1.;
	if (abs(*xi) > val || abs(*eta) > val || abs(*zeta) > val) {
	    ierror = 3;
	}
	*itest = errmes_(itest, &ierror, srname, 6L);
	if (*itest != 0) {
	    return 0;
	}
    }

/*     Set shape functions */

    fun[1] = (*xi * m1 + 1.) * .125 * (*eta * m1 + 1.) * (*zeta * m1 + 1.) * (
	    *xi * m1 + *eta * m1 + *zeta * m1 - 2.);
    fun[2] = (*xi * m1 + 1.) * .25 * (1. - *eta * *eta) * (*zeta * m1 + 1.);
    fun[3] = (*xi * m1 + 1.) * .125 * (*eta * p1 + 1.) * (*zeta * m1 + 1.) * (
	    *xi * m1 + *eta * p1 + *zeta * m1 - 2.);
    fun[4] = (1. - *xi * *xi) * .25 * (*eta * p1 + 1.) * (*zeta * m1 + 1.);
    fun[5] = (*xi * p1 + 1.) * .125 * (*eta * p1 + 1.) * (*zeta * m1 + 1.) * (
	    *xi * p1 + *eta * p1 + *zeta * m1 - 2.);
    fun[6] = (*xi * p1 + 1.) * .25 * (1. - *eta * *eta) * (*zeta * m1 + 1.);
    fun[7] = (*xi * p1 + 1.) * .125 * (*eta * m1 + 1.) * (*zeta * m1 + 1.) * (
	    *xi * p1 + *eta * m1 + *zeta * m1 - 2.);
    fun[8] = (1. - *xi * *xi) * .25 * (*eta * m1 + 1.) * (*zeta * m1 + 1.);
    fun[9] = (*xi * m1 + 1.) * .25 * (*eta * m1 + 1.) * (1. - *zeta * *zeta);
    fun[10] = (*xi * m1 + 1.) * .25 * (*eta * p1 + 1.) * (1. - *zeta * *zeta);

    fun[11] = (*xi * p1 + 1.) * .25 * (*eta * p1 + 1.) * (1. - *zeta * *zeta);

    fun[12] = (*xi * p1 + 1.) * .25 * (*eta * m1 + 1.) * (1. - *zeta * *zeta);

    fun[13] = (*xi * m1 + 1.) * .125 * (*eta * m1 + 1.) * (*zeta * p1 + 1.) * 
	    (*xi * m1 + *eta * m1 + *zeta * p1 - 2.);
    fun[14] = (*xi * m1 + 1.) * .25 * (1. - *eta * *eta) * (*zeta * p1 + 1.);
    fun[15] = (*xi * m1 + 1.) * .125 * (*eta * p1 + 1.) * (*zeta * p1 + 1.) * 
	    (*xi * m1 + *eta * p1 + *zeta * p1 - 2.);
    fun[16] = (1. - *xi * *xi) * .25 * (*eta * p1 + 1.) * (*zeta * p1 + 1.);
    fun[17] = (*xi * p1 + 1.) * .125 * (*eta * p1 + 1.) * (*zeta * p1 + 1.) * 
	    (*xi * p1 + *eta * p1 + *zeta * p1 - 2.);
    fun[18] = (*xi * p1 + 1.) * .25 * (1. - *eta * *eta) * (*zeta * p1 + 1.);
    fun[19] = (*xi * p1 + 1.) * .125 * (*eta * m1 + 1.) * (*zeta * p1 + 1.) * 
	    (*xi * p1 + *eta * m1 + *zeta * p1 - 2.);
    fun[20] = (1. - *xi * *xi) * .25 * (*eta * m1 + 1.) * (*zeta * p1 + 1.);

/*     Set derivatives */

    der[der_dim1 + 1] = m1 * .125 * (*eta * m1 + 1.) * (*zeta * m1 + 1.) * (*
	    xi * 2. * m1 + *eta * m1 + *zeta * m1 - 1.);
    der[der_dim1 + 2] = m1 * .125 * (*xi * m1 + 1.) * (*zeta * m1 + 1.) * (*
	    xi * m1 + *eta * 2. * m1 + *zeta * m1 - 1.);
    der[der_dim1 + 3] = m1 * .125 * (*xi * m1 + 1.) * (*eta * m1 + 1.) * (*xi 
	    * m1 + *eta * m1 + *zeta * 2. * m1 - 1.);
    der[(der_dim1 << 1) + 1] = m1 * .25 * (1. - *eta * *eta) * (*zeta * m1 + 
	    1.);
    der[(der_dim1 << 1) + 2] = *eta * -.5 * (*xi * m1 + 1.) * (*zeta * m1 + 
	    1.);
    der[(der_dim1 << 1) + 3] = m1 * .25 * (1. - *eta * *eta) * (*xi * m1 + 1.)
	    ;
    der[der_dim1 * 3 + 1] = m1 * .125 * (*eta * p1 + 1.) * (*zeta * m1 + 1.) *
	     (*xi * 2. * m1 + *eta * p1 + *zeta * m1 - 1.);
    der[der_dim1 * 3 + 2] = p1 * .125 * (*xi * m1 + 1.) * (*zeta * m1 + 1.) * 
	    (*xi * m1 + *eta * 2. * p1 + *zeta * m1 - 1.);
    der[der_dim1 * 3 + 3] = m1 * .125 * (*xi * m1 + 1.) * (*eta * p1 + 1.) * (
	    *xi * m1 + *eta * p1 + *zeta * 2. * m1 - 1.);
    der[(der_dim1 << 2) + 1] = *xi * -.5 * (*eta * p1 + 1.) * (*zeta * m1 + 
	    1.);
    der[(der_dim1 << 2) + 2] = p1 * .25 * (1. - *xi * *xi) * (*zeta * m1 + 1.)
	    ;
    der[(der_dim1 << 2) + 3] = m1 * .25 * (1. - *xi * *xi) * (*eta * p1 + 1.);

    der[der_dim1 * 5 + 1] = p1 * .125 * (*eta * p1 + 1.) * (*zeta * m1 + 1.) *
	     (*xi * 2. * p1 + *eta * p1 + *zeta * m1 - 1.);
    der[der_dim1 * 5 + 2] = p1 * .125 * (*xi * p1 + 1.) * (*zeta * m1 + 1.) * 
	    (*xi * p1 + *eta * 2. * p1 + *zeta * m1 - 1.);
    der[der_dim1 * 5 + 3] = m1 * .125 * (*xi * p1 + 1.) * (*eta * p1 + 1.) * (
	    *xi * p1 + *eta * p1 + *zeta * 2. * m1 - 1.);
    der[der_dim1 * 6 + 1] = p1 * .25 * (1. - *eta * *eta) * (*zeta * m1 + 1.);

    der[der_dim1 * 6 + 2] = *eta * -.5 * (*xi * p1 + 1.) * (*zeta * m1 + 1.);
    der[der_dim1 * 6 + 3] = m1 * .25 * (1. - *eta * *eta) * (*xi * p1 + 1.);
    der[der_dim1 * 7 + 1] = p1 * .125 * (*eta * m1 + 1.) * (*zeta * m1 + 1.) *
	     (*xi * 2. * p1 + *eta * m1 + *zeta * m1 - 1.);
    der[der_dim1 * 7 + 2] = m1 * .125 * (*xi * p1 + 1.) * (*zeta * m1 + 1.) * 
	    (*xi * p1 + *eta * 2. * m1 + *zeta * m1 - 1.);
    der[der_dim1 * 7 + 3] = m1 * .125 * (*xi * p1 + 1.) * (*eta * m1 + 1.) * (
	    *xi * p1 + *eta * m1 + *zeta * 2. * m1 - 1.);
    der[(der_dim1 << 3) + 1] = *xi * -.5 * (*eta * m1 + 1.) * (*zeta * m1 + 
	    1.);
    der[(der_dim1 << 3) + 2] = m1 * .25 * (1. - *xi * *xi) * (*zeta * m1 + 1.)
	    ;
    der[(der_dim1 << 3) + 3] = m1 * .25 * (1. - *xi * *xi) * (*eta * m1 + 1.);

    der[der_dim1 * 9 + 1] = m1 * .25 * (*eta * m1 + 1.) * (1. - *zeta * *zeta)
	    ;
    der[der_dim1 * 9 + 2] = m1 * .25 * (*xi * m1 + 1.) * (1. - *zeta * *zeta);

    der[der_dim1 * 9 + 3] = *zeta * -.5 * (*xi * m1 + 1.) * (*eta * m1 + 1.);
    der[der_dim1 * 10 + 1] = m1 * .25 * (*eta * p1 + 1.) * (1. - *zeta * *
	    zeta);
    der[der_dim1 * 10 + 2] = p1 * .25 * (*xi * m1 + 1.) * (1. - *zeta * *zeta)
	    ;
    der[der_dim1 * 10 + 3] = *zeta * -.5 * (*xi * m1 + 1.) * (*eta * p1 + 1.);

    der[der_dim1 * 11 + 1] = p1 * .25 * (*eta * p1 + 1.) * (1. - *zeta * *
	    zeta);
    der[der_dim1 * 11 + 2] = p1 * .25 * (*xi * p1 + 1.) * (1. - *zeta * *zeta)
	    ;
    der[der_dim1 * 11 + 3] = *zeta * -.5 * (*xi * p1 + 1.) * (*eta * p1 + 1.);

    der[der_dim1 * 12 + 1] = p1 * .25 * (*eta * m1 + 1.) * (1. - *zeta * *
	    zeta);
    der[der_dim1 * 12 + 2] = m1 * .25 * (*xi * p1 + 1.) * (1. - *zeta * *zeta)
	    ;
    der[der_dim1 * 12 + 3] = *zeta * -.5 * (*xi * p1 + 1.) * (*eta * m1 + 1.);

    der[der_dim1 * 13 + 1] = m1 * .125 * (*eta * m1 + 1.) * (*zeta * p1 + 1.) 
	    * (*xi * 2. * m1 + *eta * m1 + *zeta * p1 - 1.);
    der[der_dim1 * 13 + 2] = m1 * .125 * (*xi * m1 + 1.) * (*zeta * p1 + 1.) *
	     (*xi * m1 + *eta * 2. * m1 + *zeta * p1 - 1.);
    der[der_dim1 * 13 + 3] = p1 * .125 * (*xi * m1 + 1.) * (*eta * m1 + 1.) * 
	    (*xi * m1 + *eta * m1 + *zeta * 2. * p1 - 1.);
    der[der_dim1 * 14 + 1] = m1 * .25 * (1. - *eta * *eta) * (*zeta * p1 + 1.)
	    ;
    der[der_dim1 * 14 + 2] = *eta * -.5 * (*xi * m1 + 1.) * (*zeta * p1 + 1.);

    der[der_dim1 * 14 + 3] = p1 * .25 * (1. - *eta * *eta) * (*xi * m1 + 1.);
    der[der_dim1 * 15 + 1] = m1 * .125 * (*eta * p1 + 1.) * (*zeta * p1 + 1.) 
	    * (*xi * 2. * m1 + *eta * p1 + *zeta * p1 - 1.);
    der[der_dim1 * 15 + 2] = p1 * .125 * (*xi * m1 + 1.) * (*zeta * p1 + 1.) *
	     (*xi * m1 + *eta * 2. * p1 + *zeta * p1 - 1.);
    der[der_dim1 * 15 + 3] = p1 * .125 * (*xi * m1 + 1.) * (*eta * p1 + 1.) * 
	    (*xi * m1 + *eta * p1 + *zeta * 2. * p1 - 1.);
    der[(der_dim1 << 4) + 1] = *xi * -.5 * (*eta * p1 + 1.) * (*zeta * p1 + 
	    1.);
    der[(der_dim1 << 4) + 2] = p1 * .25 * (1. - *xi * *xi) * (*zeta * p1 + 1.)
	    ;
    der[(der_dim1 << 4) + 3] = p1 * .25 * (1. - *xi * *xi) * (*eta * p1 + 1.);

    der[der_dim1 * 17 + 1] = p1 * .125 * (*eta * p1 + 1.) * (*zeta * p1 + 1.) 
	    * (*xi * 2. * p1 + *eta * p1 + *zeta * p1 - 1.);
    der[der_dim1 * 17 + 2] = p1 * .125 * (*xi * p1 + 1.) * (*zeta * p1 + 1.) *
	     (*xi * p1 + *eta * 2. * p1 + *zeta * p1 - 1.);
    der[der_dim1 * 17 + 3] = p1 * .125 * (*xi * p1 + 1.) * (*eta * p1 + 1.) * 
	    (*xi * p1 + *eta * p1 + *zeta * 2. * p1 - 1.);
    der[der_dim1 * 18 + 1] = p1 * .25 * (1. - *eta * *eta) * (*zeta * p1 + 1.)
	    ;
    der[der_dim1 * 18 + 2] = *eta * -.5 * (*xi * p1 + 1.) * (*zeta * p1 + 1.);

    der[der_dim1 * 18 + 3] = p1 * .25 * (1. - *eta * *eta) * (*xi * p1 + 1.);
    der[der_dim1 * 19 + 1] = p1 * .125 * (*eta * m1 + 1.) * (*zeta * p1 + 1.) 
	    * (*xi * 2. * p1 + *eta * m1 + *zeta * p1 - 1.);
    der[der_dim1 * 19 + 2] = m1 * .125 * (*xi * p1 + 1.) * (*zeta * p1 + 1.) *
	     (*xi * p1 + *eta * 2. * m1 + *zeta * p1 - 1.);
    der[der_dim1 * 19 + 3] = p1 * .125 * (*xi * p1 + 1.) * (*eta * m1 + 1.) * 
	    (*xi * p1 + *eta * m1 + *zeta * 2. * p1 - 1.);
    der[der_dim1 * 20 + 1] = *xi * -.5 * (*eta * m1 + 1.) * (*zeta * p1 + 1.);

    der[der_dim1 * 20 + 2] = m1 * .25 * (1. - *xi * *xi) * (*zeta * p1 + 1.);
    der[der_dim1 * 20 + 3] = p1 * .25 * (1. - *xi * *xi) * (*eta * m1 + 1.);

} /* brk20_ */

